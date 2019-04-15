import requests
from lookups import mutation_type as mut
import logging
import re
import copy
from normalizers import gene_enricher

from lookups.prot_to_hgvs import to_hgvs_p


def _enrich_ensemble(feature, transcript_id, exon, provenance_rule):
    """ get coordinates from ensembl
        curl -s 'http://grch37.rest.ensembl.org/lookup/id/ENST00000275493?expand=1'
        | jq '{start:.start, end:.end, chromosome:.seq_region_name, strand:.strand}'
    """
    headers = {'Content-type': 'application/json'}
    url = 'http://grch37.rest.ensembl.org/lookup/id/{}?expand=1' \
        .format(transcript_id)
    r = requests.get(url, timeout=60, headers=headers)
    transcript = r.json()

    if 'Exon' in transcript:
        exon_ = transcript['Exon'][exon-1]
        feature['chromosome'] = str(exon_['seq_region_name'])
        feature['start'] = int(exon_['start'])
        feature['end'] = int(exon_['end'])
        feature['referenceName'] = 'GRCh37'
        if 'provenance' not in feature:
            feature['provenance'] = []
        feature['provenance'].append(url)
        feature['provenance_rule'] = provenance_rule
    return feature


def _enrich_gene(feature, gene=None, provenance_rule='default'):
    """ description contains a gene, get its location """
    if not gene:
        gene = feature['description']
    parms = 'fields=genomic_pos_hg19'
    url = "http://mygene.info/v3/query?q={}&{}".format(gene, parms)
    r = requests.get(url, timeout=60)
    hit = None
    hits = r.json()
    if 'hits' in hits:
        for a_hit in hits['hits']:
            if 'genomic_pos_hg19' in a_hit:
                hit = a_hit['genomic_pos_hg19']
                break

    if isinstance(hit, list):
        alternatives = hit
        for alternative in alternatives:
            if alternative['chr'] in ['20', '21', '22', '23', '1', '3', '2',
                                      '5', '4', '7', '6', '9', '8', 'Y', 'X',
                                      '11', '10', '13', '12', '15', '14', '17',
                                      '16', '19', '18']:
                hit = alternative
    if hit:
        if 'chr' in hit and ('chromosome' not in feature or not feature['chromosome']):
            feature['chromosome'] = str(hit['chr'])
        if 'start' in hit:
            feature['start'] = hit['start']
        if 'end' in hit:
            feature['end'] = hit['end']
        feature['referenceName'] = 'GRCh37'
        if 'provenance' not in feature:
            feature['provenance'] = []
        feature['provenance'].append(url)
        feature['provenance_rule'] = provenance_rule

    # FIXME: this function is called for variants like 'exon 18 deletion', but it gets the start, end coords of the entire gene
    # perhaps it should get the start, end coords of exon 18 instead?

    return feature


def fetch_variant_pos(gene_symbol, name, loose_search=True):
    """
    Fetches positional info for a variant from myvariant.info.

    Uses snpeff's genename and hgvs_p fields to perform the matching. If conversion to an hgvs_p string
    fails and loose_search is True, uses the variant name in a loose all-fields search. Otherwise, returns
    None.
    :param gene_symbol: the hugo symbol for the variant, e.g. BRAF
    :param name: the variant's colloquial name, e.g. V600E
    :param loose_search: whether to search for the variant anyway if HGVS protein string conversion fails
    :return: a tuple consisting of the hit, if available (None otherwise), and the URL that was attempted
    """

    # NOTE: the below was the original query, which gets incorrect variants and thus produces incorrect coordinates(!!)
    # url = "http://myvariant.info/v1/query?q={}".format(feature['description'])
    # this is a fixed URL which queries specifically for snpeff entries by genename and protein-level change
    url = "http://myvariant.info/v1/query"
    params = {'fields': 'hg19,vcf,snpeff,chrom'}
    remapped = to_hgvs_p(name)

    if remapped is not None:
        params['q'] = 'snpeff.ann.genename:%s AND snpeff.ann.hgvs_p:%s' % (gene_symbol, remapped)
    elif loose_search:
        params['q'] = 'snpeff.ann.genename:%s AND %s' % (gene_symbol, name)
    else:
        return None

    r = requests.get(url, params=params, timeout=60)
    hits = r.json()

    if 'hits' in hits:
        for a_hit in hits['hits']:
            if 'hg19' in a_hit and 'vcf' in a_hit:
                return a_hit, r.url

    # either no hits were found, or none contained the necessary positional information
    return None, r.url


def _enrich_feature(feature, provenance_rule='default'):
    """ description contains a gene + variant, get its location """
    #  curl -s http://myvariant.info/v1/query?q=FLT3%20N676D |
    # jq '.hits[0] |
    # { referenceName: "GRCh37", chromosome: .chrom,
    # start: .hg19.start, end: .hg19.end, ref: .vcf.ref, alt: .vcf.alt  }'
    # {
    #   "referenceName": "GRCh37",
    #   "chromosome": "13",
    #   "start": 28644637,
    #   "end": 28644637,
    #   "ref": "T",
    #   "alt": "A"
    # }

    hit, url = fetch_variant_pos(feature['geneSymbol'], feature['name'])

    if hit:
        hg19 = hit.get('hg19')
        vcf = hit.get('vcf')
        if 'ref' in vcf:
            feature['ref'] = vcf['ref']
        if 'alt' in vcf:
            feature['alt'] = vcf['alt']
        if 'chrom' in hit:
            feature['chromosome'] = str(hit['chrom'])
        if 'start' in hg19:
            feature['start'] = hg19['start']
        if 'end' in hg19:
            feature['end'] = hg19['end']
        feature['referenceName'] = 'GRCh37'
        if 'provenance' not in feature:
            feature['provenance'] = []
        feature['provenance'].append(url)
        feature['provenance_rule'] = provenance_rule

        if 'biomarker_type' not in feature:
            if 'cadd' in hit and 'type' in hit['cadd']:
                feature['biomarker_type'] = \
                    mut.norm_biomarker(hit['cadd']['type'])

    return feature


def _is_gene(symbols):
    """ return true if all symbols exist"""
    for symbol in symbols:
        try:
            # this returns a value, but we just care if it throws an exception
            gene_enricher.get_gene(symbol)
        except ValueError:
            return False
    return True


def enrich(feature, feature_association):
    """
    given a feature, decorate it with genomic location
    """
    enriched_features = [feature]

    # noinspection PyBroadException
    try:
        # return if already there
        # ^ (FAISAL: i assume this means if we already have a start position, although the other fields might not be filled)
        if feature.get('start', None):
            feature['provenance_rule'] = 'from_source'
            return [feature]

        # make sure it has a name and a description
        if not feature.get('description', None):
            feature['description'] = feature.get('name', None)
        if not feature.get('name', None):
            feature['name'] = feature.get('description', None)

        # we can't normalize things without a description
        if not feature.get('description', None):
            feature['provenance_rule'] = 'missing_description'
            return [feature]

        # apply rules
        description_parts = re.split(' +', feature['description'].strip())
        source = feature_association['source'] if 'source' in feature_association else None
        exonMatch = re.match(r'.* Exon ([0-9]*) .*', feature['name'], re.M|re.I)  # re.M,I: multiline, ignore case
        feat_desc_lc = feature['description'].lower()

        enriched_features = []

        # we know the description has at least [0] == gene symbol

        if not _is_gene([description_parts[1]]) and len(description_parts[1].split('-')) == 2 and _is_gene(description_parts[1].split('-')):
            # it's a fusion, e.g. BRAF-EGFR etc.
            # we determine that by first checking if the thing isn't a gene as a whole ("BRAF-EGFR" is of course not),
            # then by checking if we get exactly two pieces when we split it on the hyphen,
            # then by checking if each piece *is* a gene

            # split out the fusion donor, acceptor and enrich gene information for each one individually
            fusion_donor, fusion_acceptor = description_parts[1].split('-')

            feature_fusion_donor = _enrich_gene(copy.deepcopy(feature), fusion_donor, provenance_rule='is_fusion_donor')
            feature_fusion_donor['geneSymbol'] = fusion_donor
            enriched_features.append(feature_fusion_donor)

            feature_fusion_acceptor = _enrich_gene(copy.deepcopy(feature), fusion_acceptor, provenance_rule='is_fusion_acceptor')
            feature_fusion_acceptor['geneSymbol'] = fusion_acceptor
            enriched_features.append(feature_fusion_acceptor)

        elif len(description_parts) == 1:
            # it's just a gene name, e.g. "EGFR"; we can only report info about the gene
            feature = _enrich_gene(feature, provenance_rule='gene_only')
            enriched_features.append(feature)

        elif ('oncokb' == source and 'clinical' in feature_association['oncokb'] and
              exonMatch and 'Isoform' in feature_association['oncokb']['clinical']):
            # it's an oncokb clinical entry which mentions an exon and contains isoform information
            # e.g., "EGFR Exon 19
            # (which we need to determine exon mappings)

            isoform = feature_association['oncokb']['clinical']['Isoform']
            feature = _enrich_ensemble(feature, transcript_id=isoform, exon=int(exonMatch.group(1)), provenance_rule='is_oncokb_exon')
            enriched_features.append(feature)

        elif 'deletion' in feat_desc_lc or 'del ' in feat_desc_lc:
            # FIXME: some civic variants report the representative transcript and start, end coords of the exon
            # we should be reporting that instead of the extents of the entire gene
            # but, of course, there are examples like EGFR Ex19 del L858R
            feature = _enrich_gene(feature, description_parts[0], provenance_rule='is_deletion')
            enriched_features.append(feature)

        elif 'amplification' in feat_desc_lc or 'amp ' in feat_desc_lc:
            # this is another gene-level event, e.g. "BRAF AMPLIFICATION"
            # we can only report gene-level information
            feature = _enrich_gene(feature, description_parts[0], provenance_rule='is_amplification')
            enriched_features.append(feature)

        elif 'loss' in feat_desc_lc:
            feature = _enrich_gene(feature, description_parts[0], provenance_rule='is_loss')
            enriched_features.append(feature)

        elif 'mutation' in feat_desc_lc:
            feature = _enrich_gene(feature, description_parts[0], provenance_rule='is_mutation')
            enriched_features.append(feature)

        elif 'inact mut' in feat_desc_lc:
            feature = _enrich_gene(feature, description_parts[0], provenance_rule='is_inact_mut')
            enriched_features.append(feature)

        else:
            feature = _enrich_feature(feature, provenance_rule='default_feature')
            enriched_features.append(feature)

    except Exception as e:
        logging.error(feature, exc_info=1)

    return enriched_features
