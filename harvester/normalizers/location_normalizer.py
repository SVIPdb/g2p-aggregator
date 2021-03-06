import traceback

import requests
import re
import logging
import copy

import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.parser
from hgvs.exceptions import HGVSError
from hgvs.sequencevariant import SequenceVariant
import hgvs.dataproviders.uta
import hgvs.assemblymapper

from lookups.hgvs_parsing import hgvsparser, am
from normalizers.feature_enricher import enrich

from normalizers.reference_genome_normalizer import normalize as normalize_referencename

from utils_ex.instrumentation import add_crawl_status

# these shared assembly mappers will allow us to convert HGVS g. variants to c. and p. later on


def _complement(bases, reverse=True):
    """
    return complement of bases string
    """
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if reverse:
        return ''.join([complements.get(base, base) for base in reversed(bases)])
    else:
        return ''.join([complements.get(base, base) for base in bases])


def _get_ref_alt(description):
    """
    parse ref & alt from description
    return none if no match
    """
    p = ".*[c]\.[0-9]*_[0-9]*(.*)"
    m = re.match(p, description)
    if m and len(m.groups()) == 1:
        return m.groups()[0]
    p = ".*[A-Z][0-9]*_[A-Z][0-9]*(.*)"
    m = re.match(p, description)
    if m and len(m.groups()) == 1:
        return m.groups()[0]
    return None


def allele_registry(hgvs):
    """
    call allele registry with hgvs notation, return allele_registry
    """
    url = 'http://reg.genome.network/allele?hgvs={}' \
        .format(requests.utils.quote(hgvs))
    r = requests.get(url, headers={'Accept': 'application/json'})
    if r.status_code not in [200, 400, 404]:
        logging.info('unexpected allele_registry {} {}'.format(url,
                                                               r.status_code))
    rsp = r.json()
    return rsp, url


def genomic_hgvs(feature, complement=False, description=False):
    """
    given a feature, create a hgvs genomic representation
    http://varnomen.hgvs.org/bg-material/refseq/#DNAg
    """
    ac_map = {
        '1': 'NC_000001.10',
        '2': 'NC_000002.11',
        '3': 'NC_000003.11',
        '4': 'NC_000004.11',
        '5': 'NC_000005.9',
        '6': 'NC_000006.11',
        '7': 'NC_000007.13',
        '8': 'NC_000008.10',
        '9': 'NC_000009.11',
        '10': 'NC_000010.10',
        '11': 'NC_000011.9',
        '12': 'NC_000012.11',
        '13': 'NC_000013.10',
        '14': 'NC_000014.8',
        '15': 'NC_000015.9',
        '16': 'NC_000016.9',
        '17': 'NC_000017.10',
        '18': 'NC_000018.9',
        '19': 'NC_000019.9',
        '20': 'NC_000020.10',
        '21': 'NC_000021.8',
        '22': 'NC_000022.10',
        'X': 'NC_000023.10',
        '23': 'NC_000023.10',
        'Y': 'NC_000024.9',
    }
    ac = ac_map[feature['chromosome']]

    # Make an edit object
    ref = feature.get('ref', None)
    if ref == '-' or ref == '':
        ref = None
    alt = feature.get('alt', None)
    if alt == '-' or alt == '':
        alt = None

    if complement:
        ref = _complement(ref) if ref else None
        alt = _complement(alt) if alt else None

    if 'start' in feature and feature['start']:
        start = hgvs.location.SimplePosition(base=int(feature['start']))
    else:
        start = None

    if 'end' in feature and feature['end']:
        end = hgvs.location.SimplePosition(base=int(feature['end']))
    else:
        # if 'end' isn't specified, we can kind of assume it's going to cover the reference, at least
        # (it's -1 because a SNP's start == end)
        end = hgvs.location.SimplePosition(base=int(feature['start']) + max(0, len(ref)-1))

    iv = hgvs.location.Interval(start=start, end=end)

    feature_description = feature.get('description', feature.get('name', None))
    edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
    posedit = hgvs.posedit.PosEdit(pos=iv, edit=edit)

    if description:
        # override with ref_alt from description or name
        # override edit if we have a hint from description that its a dup
        ref_alt = _get_ref_alt(feature_description)
        posedit.edit = ref_alt
    else:
        if 'dup' in feature_description:
            ref_alt = _get_ref_alt(feature_description)
            if ref_alt:
                posedit.edit = ref_alt + alt

    # Make the variant
    var = SequenceVariant(ac=ac, type='g', posedit=posedit)
    # https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37.p13

    return var


def _key_na_in_dict(k, o):
    return k not in o or o[k] is None or o[k] == 'None'


def normalize(feature):
    # if we don't have referenceName (aka NCBI refseq id), the chromosome, or the reference sequence,
    # we can't do anything
    if any(_key_na_in_dict(k, feature) for k in ('referenceName', 'chromosome', 'ref')):
        return None, None

    # if strand info is present (e.g., from COSMIC), use it to skip a likely round-trip w/allele registry
    inferred_complement = (feature.get('strand') == '-')

    # if hgvs_g has already been specified, use that, otherwise attempt to compute it ourselves
    hgvs_rep = hgvsparser.parse_g_variant(feature['hgvs_g']) if 'hgvs_g' in feature and feature['hgvs_g'] \
        else genomic_hgvs(feature, complement=inferred_complement)
    allele = None
    provenance = None

    if hgvs_rep:
        (allele, provenance) = allele_registry(str(hgvs_rep))

        if 'errorType' in allele and allele['errorType'] == 'IncorrectReferenceAllele':
            message = allele['message']
            actualAllele = allele['actualAllele']

            complement_ref = _complement(feature['ref'])
            if complement_ref == actualAllele:
                # print 'reverse strand re-try'
                # try the opposite of what we tried before
                hgvs_rep = genomic_hgvs(feature, complement=not inferred_complement)
                (allele, provenance) = allele_registry(str(hgvs_rep))
            # else:
            #     print 'complement_ref {} m[0] {}'.format(complement_ref,
            #                                              actualAllele)

        if 'errorType' in allele and allele['errorType'] == 'IncorrectHgvsPosition':
            # print 'position error re-try'
            hgvs_rep = genomic_hgvs(feature, description=True)
            (allele, provenance) = allele_registry(str(hgvs_rep))

        if allele:
            # hgvs_rep = hgnorm.normalize(hgvs_rep)
            allele['hgvs_g'] = str(hgvs_rep)

            # FIXME: it's quite possible that the allele registry has returned hgvs strings for us in the genomicAlleles,
            #  and transcriptAlleles keys. we should probably consult those if HGVS remapping fails, e.g. due to a missing
            #  refseq

            # also use the variant mapper to convert hgvs_rep to the c. and p. versions as well, if possible
            # to accomplish that, we'll need the reference sequence. in this case, the hgvs coordinates are more
            # stringently validated than for genomic coords (e.g., delins position coords have to match the edit length)
            normal_reference_name = normalize_referencename(feature['referenceName'])
            if feature.get('refseq'):
                try:
                    # FIXME: we get a few HGVS.g strings that look like NC_000007.13:g.140453135delinsAT,
                    #  which doesn't specify a range for the deletion and thus isn't a valid HGVS string.
                    #  we should track down where these are coming from and fix it.
                    if not feature.get('hgvs_c'):
                        var_c = am[normal_reference_name].g_to_c(hgvs_rep, tx_ac=feature['refseq'])
                        allele['hgvs_c'] = str(var_c)
                    else:
                        var_c = hgvsparser.parse_c_variant(feature.get('hgvs_c'))
                        allele['hgvs_c'] = str(var_c)

                    # now that we have an hgvs.c, let's see if we can get an hgvs.p
                    if var_c and not feature.get('hgvs_p'):
                        var_p = am[normal_reference_name].c_to_p(var_c)
                        if var_p.posedit:
                            # FIXME: verify that we should be flagging this as certain
                            var_p.posedit.uncertain = False
                        allele['hgvs_p'] = str(var_p)

                except HGVSError:
                    message = (
                        "Couldn't produce HGVS c. and p. strings for g. string '%s' with assembly/refseq '%s'/'%s'" %
                        (allele['hgvs_g'], normal_reference_name, feature['refseq'])
                    )
                    edit_msg = "Edit: %s, %s" % (repr(hgvs_rep.posedit.edit), str(hgvs_rep.posedit.edit))
                    traceback_msg = traceback.format_exc()

                    logging.warn(message)
                    logging.warn(edit_msg)
                    logging.warn(traceback_msg)

                    add_crawl_status(feature, __name__, {
                        'message': message,
                        'edit': edit_msg,
                        'traceback': traceback_msg
                    })

    return allele, provenance


def _apply_allele_registry(feature, allele_registry, provenance):
    # there is a lot of info in registry, just get synonyms and links
    links = feature.get('links', [])
    synonyms = feature.get('synonyms', [])
    links.append(allele_registry['@id'])

    if 'externalRecords' in allele_registry:
        externalRecords = allele_registry['externalRecords']
        links = links + [externalRecords[r][0].get('@id')
                         for r in externalRecords
                         if '@id' in externalRecords[r][0]
                         ]
        synonyms = synonyms + [externalRecords[r][0].get('id')
                               for r in externalRecords
                               if 'id' in externalRecords[r][0]
                               ]

        # get the dbsnp IDs and myvariant URL from the allele registry, too, if possible
        if 'dbSNP' in externalRecords:
            feature['dbsnp_ids'] = [str(x['rs']) for x in externalRecords['dbSNP']]
        if 'MyVariantInfo_hg19' in externalRecords and len(externalRecords['MyVariantInfo_hg19']) > 0:
            feature['myvariant_hg19'] = externalRecords['MyVariantInfo_hg19'][0]['@id']

    if 'genomicAlleles' in allele_registry:
        for genomicAllele in allele_registry['genomicAlleles']:
            synonyms = synonyms + genomicAllele['hgvs']
            links.append(genomicAllele['referenceSequence'])

    synonyms = list(set(synonyms))
    links = list(set(links))
    if len(synonyms) > 0:
        feature['synonyms'] = synonyms
    if len(links) > 0:
        feature['links'] = links
    if 'provenance' not in feature:
        feature['provenance'] = []
    feature['provenance'].append(provenance)

    # capture the extremely useful hgvs_g field from the allele registry as well
    if _key_na_in_dict('hgvs_g', feature):
        feature['hgvs_g'] = allele_registry['hgvs_g']

    # also try to get hgvs_c, hgvs_p, but they may not be available if assemblymapper fails
    if _key_na_in_dict('hgvs_c', feature):
        feature['hgvs_c'] = allele_registry.get('hgvs_c')

    if _key_na_in_dict('hgvs_p', feature):
        feature['hgvs_p'] = allele_registry.get('hgvs_p')


def _fix_location_end(feature):
    """ if end not present, set it based on start, ref & alt length"""
    end = feature.get('end', 0)
    start = feature.get('start', 0)
    ref_len = 0
    alt_len = 0
    if feature.get('ref', None):
        ref_len = len(feature.get('ref', ''))
    if feature.get('alt', None):
        alt_len = len(feature.get('alt', ''))
    offset = max(ref_len, alt_len)
    if start > 0 and end == 0:
        end = max(start, start + (offset - 1))
        feature['end'] = end
    return feature


def normalize_feature_association(feature_association):
    """ given the 'final' g2p feature_association,
    update it with genomic location """
    allele_registry_instance = None
    normalized_features = []

    for feature in feature_association['features']:
        # skip features that already have hgvs_g, hgvs_c, and hgvs_p fields
        if all(not _key_na_in_dict(k, feature) for k in ('hgvs_g', 'hgvs_c', 'hgvs_p')):
            continue

        try:
            # ensure we have location, enrich can create new features
            enriched_features = enrich(copy.deepcopy(feature), feature_association)

            for enriched_feature in enriched_features:
                # go get AR info
                (allele_registry_instance, provenance) = normalize(enriched_feature)

                if allele_registry_instance and '@id' in allele_registry_instance:
                    _apply_allele_registry(enriched_feature, allele_registry_instance, provenance)

                enriched_feature = _fix_location_end(enriched_feature)
                normalized_features.append(enriched_feature)

            feature_association['features'] = normalized_features

        except Exception as e:
            message = 'exception {}'.format(e)
            logging.exception(message)
            add_crawl_status(feature, __name__, {
                'message': message
            })


def _test(feature):
    allele_reg = normalize(feature)
    if allele_reg and '@id' not in allele_reg:
        print 'FAIL', allele_reg['message']
        print "\t", allele_reg['hgvs_g']
        print "\t", feature
        print "\t", allele_reg
    elif allele_reg:
        print 'OK', allele_reg['hgvs_g']
    else:
        print 'OK', 'not normalized'


if __name__ == '__main__':
    import yaml
    import logging.config
    path = 'logging.yml'
    with open(path) as f:
        config = yaml.load(f)
    logging.config.dictConfig(config)

    # _test('YAP1-MAMLD1 Fusion')
    # _test('YAP1-MAMLD1')
    _test({
      "entrez_id": 238,
      "end": 29445213,
      "name": "HIP1-ALK I1171N",
      "start": 29445213,
      "biomarker_type": "snp",
      "referenceName": "GRCh37",
      "geneSymbol": "ALK",
      "alt": "T",
      "ref": "A",
      "chromosome": "2"
    })

    _test({"end": "55242483",
           "name": "EGFR c.2237_2253delAATTAAGAGAAGCAACAinsTC",
           "start": "55242467", "biomarker_type": "snp",
           "referenceName": "GRCh37", "alt": "TC", "ref": "-",
           "chromosome": "7",
           "description": "EGFR c.2237_2253delAATTAAGAGAAGCAACAinsTC"})

    _test({"entrez_id": 9968, "end": 70349258, "name": "L1224F",
           "start": 70349258, "biomarker_type": "snp",
           "referenceName": "GRCh37", "geneSymbol": "MED12",
           "alt": "T", "ref": "C", "chromosome": "23"})

    _test({"name": "ABL1:V289I", "start": 133747558,
           "biomarker_type": "mutant", "referenceName": "GRCh37",
           "geneSymbol": "ABL1", "alt": "A", "ref": "G", "chromosome": "9",
           "description": "ABL1:I242T,M244V,K247R,L248V,G250E,G250R,Q252R,Q252H,Y253F,Y253H,E255K,E255V,M237V,E258D,W261L,L273M,E275K,E275Q,D276G,T277A,E279K,V280A,V289A,V289I,E292V,E292Q,I293V,L298V,V299L,F311L,F311I,T315I,F317L,F317V,F317I,F317C,Y320C,L324Q,Y342H,M343T,A344V,A350V,M351T,E355D,E355G,E355A,F359V,F359I,F359C,F359L,D363Y,L364I,A365V,A366G,L370P,V371A,E373K,V379I,A380T,F382L,L384M,L387M,L387F,L387V,M388L,Y393C,H396P,H396R,H396A,A397P,S417F,S417Y,I418S,I418V,A433T,S438C,E450K,E450G,E450A,E450V,E453K,E453G,E453A,E453V,E459K,E459G,E459A,E459V,M472I,P480L,F486S,E507G"})  # NOQA

    _test({"end": 28592642, "name": "FLT3 D835N", "referenceName": "GRCh37",
           "start": 28592642, "biomarker_type": "snp", "geneSymbol": "FLT3",
           "attributes": {"amino_acid_change": {"string_value": "D835N"},
                            "germline": {"string_value": None},
                            "partner_gene": {"string_value": None},
                            "cytoband": {"string_value": None},
                            "exons": {"string_value": "20"},
                            "notes": {"string_value": None},
                            "cosmic": {"string_value": "789"},
                            "effect": {"string_value": None},
                            "cnv_type": {"string_value": None},
                            "id": {"string_value": 139},
                            "variant_type": {"string_value": "missense"},
                            "dna_change": {"string_value": "2503G>A"},
                            "codons": {"string_value": "835"},
                            "chromosome_based_cnv": {"string_value": False},
                            "transcript": {"string_value": "ENST00000241453"},
                            "description_type": {"string_value": "HGVS"},
                            "chromosome": {"string_value": None},
                            "description": {"string_value": None}},
          "alt": "A", "ref": "G", "chromosome": "13"})

    _test({"end": "55592183", "name": "KIT p.S501_A502dup",
           "start": "55592183",
           "biomarker_type": "polymorphism", "referenceName": "GRCh37",
           "alt": "CTGCCT", "ref": "-", "chromosome": "4",
           "description": "KIT p.S501_A502dup"})

    _test({"end": "37880996", "name": "ERBB2 Y772_A775dup",
           "start": "37880995",
           "biomarker_type": "nonsense", "referenceName": "GRCh37",
           "alt": "ATACGTGATGGC", "ref": "-", "chromosome": "17",
           "description": "ERBB2 Y772_A775dup"})

    _test({"end": 140453136, "name": "BRAF V600E", "referenceName": "GRCh37",
           "start": 140453136, "biomarker_type": "snp", "geneSymbol": "BRAF",
           "attributes": {"amino_acid_change": {"string_value": "V600E"},
                          "germline": {"string_value": None},
                          "partner_gene": {"string_value": None},
                          "cytoband": {"string_value": None},
                          "exons": {"string_value": "15"},
                          "notes": {"string_value": None},
                          "cosmic": {"string_value": "476"},
                          "effect": {"string_value": None},
                          "cnv_type": {"string_value": None},
                          "id": {"string_value": 112},
                          "variant_type": {"string_value": "missense"},
                          "dna_change": {"string_value": "1799T>A"},
                          "codons": {"string_value": "600"},
                          "chromosome_based_cnv": {"string_value": False},
                          "transcript": {"string_value": "ENST00000288602"},
                          "description_type": {"string_value": "HGVS"},
                          "chromosome": {"string_value": None},
                          "description": {"string_value": None}},
           "alt": "A", "ref": "T", "chromosome": "7"})

    # _test({"end": 1806125, "name": "FGFR3 F384L", "referenceName": "GRCh37", "start": 1806125, "biomarker_type": "snp", "geneSymbol": "FGFR3", "attributes": {"amino_acid_change": {"string_value": "F384L"}, "germline": {"string_value": False}, "partner_gene": {"string_value": None}, "cytoband": {"string_value": None}, "exons": {"string_value": "9"}, "notes": {"string_value": None}, "cosmic": {"string_value": None}, "effect": {"string_value": None}, "cnv_type": {"string_value": None}, "id": {"string_value": 391}, "variant_type": {"string_value": "missense"}, "dna_change": {"string_value": "1150T>C"}, "codons": {"string_value": "384"}, "chromosome_based_cnv": {"string_value": False}, "transcript": {"string_value": "ENST00000340107"}, "description_type": {"string_value": "HGVS"}, "chromosome": {"string_value": None}, "description": {"string_value": None}}, "alt": "C", "ref": "T", "chromosome": "4"})

    _test( {"entrez_id": 1956, "end": None, "name": "R776C", "start": None, "referenceName": "GRCh37", "geneSymbol": "EGFR", "alt": "None", "ref": "None", "chromosome": "None"})


    _test({"end": "55242478", "description": "EGFR E746_E749delELRE", "links": ["https://api.molecularmatch.com/v2/mutation/get?name=EGFR+E746_E749delELRE"], "start": "55242467", "biomarker_type": "nonsense", "referenceName": "GRCh37", "alt": "-", "ref": "12", "chromosome": "7", "name": "EGFR E746_E749delELRE"})
