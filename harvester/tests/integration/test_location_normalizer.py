import json
import sys

sys.path.append('.')  # NOQA
from normalizers.location_normalizer import normalize_feature_association


# noinspection PyPep8Naming
def test_normalize_feature_association_BRAF_V600E():
    fa = {
        'features': [
            {
                "geneSymbol": "BRAF",
                "name": "V600E",
                "description": "BRAF V600E"
            }
        ]
    }
    normalize_feature_association(fa)
    print(fa)

    assert fa == {
        'features': [{
            'provenance_rule': 'default_feature',
            'end': 140453155,
            'name': 'V600E',
            'links': [u'http://reg.genome.network/refseq/RS000031', u'http://reg.genome.network/allele/CA369543156',
                u'http://reg.genome.network/refseq/RS000055',
                u'http://myvariant.info/v1/variant/chr7:g.140453155C>A?assembly=hg19',
                u'http://reg.genome.network/refseq/RS000682',
                u'http://myvariant.info/v1/variant/chr7:g.140753355C>A?assembly=hg38',
                u'http://reg.genome.network/refseq/RS000007'],
            'provenance': ['http://myvariant.info/v1/query?q=BRAF V600E',
                'http://reg.genome.network/allele?hgvs=NC_000007.13%3Ag.140453155C%3EA'],
            'hgvs_g_suffix': [u'g.176410G>T', u'g.140099624C>A', u'g.140753355C>A', u'g.140453155C>A'],
            'start': 140453155,
            'synonyms': [u'ENST00000497784.1:n.1815G>T', u'XP_005250103.1:p.Asp594Tyr', u'CM000669.2:g.140753355C>A',
                u'chr7:g.140453155C>A', u'XR_927523.1:n.1703-3937G>T', u'XP_011514832.1:p.=', u'NM_004333.5:c.1780G>T',
                u'LRG_299t1:c.1780G>T', u'NM_001354609.1:c.1780G>T', u'XP_016868047.1:p.Asp634Tyr',
                u'ENST00000479537.5:n.64G>T', u'XR_927522.1:n.1703-3937G>T', u'NM_001374258.1:c.1900G>T',
                u'NM_001374244.1:c.1900G>T', u'NM_004333.4:c.1780G>T', u'XR_001744858.1:n.1823-3937G>T',
                u'XR_242190.1:n.1788G>T', u'XM_005250045.1:c.1780G>T', u'NP_001341538.1:p.Asp594Tyr',
                u'XP_011514831.1:p.Asp594Tyr', u'ENSP00000288602.6:p.Asp594Tyr', u'NM_001354609.2:c.1780G>T',
                u'XM_011516530.1:c.1695-3937G>T', u'ENST00000288602.10:c.1780G>T', u'XR_001744857.1:n.1908G>T',
                u'XM_005250046.1:c.1780G>T', u'NG_007873.3:g.176410G>T', u'XM_011516529.1:c.1780G>T',
                u'XR_927520.1:n.1788G>T', u'NP_001361187.1:p.Asp634Tyr', u'ENSP00000420119.1:p.=',
                u'CM000669.1:g.140453155C>A', u'XP_005250102.1:p.Asp594Tyr', u'XP_016868048.1:p.Asp634Tyr',
                u'chr7:g.140753355C>A', u'NM_004333.6:c.1780G>T', u'LRG_299:g.176410G>T',
                u'NC_000007.14:g.140753355C>A', u'NP_001361173.1:p.Asp634Tyr', u'ENSP00000418033.1:p.Asp22Tyr',
                u'NP_004324.2:p.Asp594Tyr', u'NC_000007.12:g.140099624C>A', '\n          ', u'XM_017012559.1:c.1900G>T',
                u'XM_017012558.1:c.1900G>T', u'NC_000007.13:g.140453155C>A', u'NR_148928.1:n.2878G>T',
                u'ENST00000496384.6:n.603G>T', u'XR_927521.1:n.1788G>T'],
            'biomarker_type': 'nonsense',
            'hgvs': [u'XP_005250103.1:p.Asp594Tyr', u'CM000669.2:g.140753355C>A', u'NP_001361173.1:p.Asp634Tyr',
                u'NC_000007.12:g.140099624C>A', u'NG_007873.3:g.176410G>T', u'ENSP00000288602.6:p.D594Y',
                u'XP_011514832.1:p.=', u'NP_001341538.1:p.D594Y', u'NP_001361187.1:p.Asp634Tyr',
                u'NC_000007.13:g.140453155C>A', u'XP_016868047.1:p.Asp634Tyr', u'XP_005250102.1:p.D594Y',
                u'NP_004324.2:p.D594Y', u'CM000669.1:g.140453155C>A', u'XP_005250102.1:p.Asp594Tyr',
                u'XP_016868048.1:p.Asp634Tyr', u'XP_016868047.1:p.D634Y', u'chr7:g.140753355C>A',
                u'XP_011514831.1:p.D594Y', u'XP_016868048.1:p.D634Y', u'LRG_299:g.176410G>T',
                u'ENSP00000418033.1:p.D22Y', u'NC_000007.14:g.140753355C>A', u'chr7:g.140453155C>A',
                u'ENSP00000418033.1:p.Asp22Tyr', u'NP_001341538.1:p.Asp594Tyr', u'XP_011514831.1:p.Asp594Tyr',
                u'ENSP00000288602.6:p.Asp594Tyr', u'NP_004324.2:p.Asp594Tyr', u'ENSP00000420119.1:p.=',
                u'XP_005250103.1:p.D594Y', u'NP_001361187.1:p.D634Y', u'NP_001361173.1:p.D634Y'],
            'referenceName': 'GRCh37',
            'hgvs_p_suffix': [u'p.D22Y', u'p.=', u'p.D594Y', u'p.D634Y', u'p.Asp22Tyr', u'p.Asp594Tyr', u'p.Asp634Tyr'],
            'geneSymbol': 'BRAF',
            'alt': u'A',
            'ref': u'C',
            'chromosome': '7',
            'description': 'BRAF V600E'
        }]
    }

# noinspection PyPep8Naming
def test_normalize_feature_association_ALK_D1203N():
    fa = {
        'features': [
            {
                "geneSymbol": "ALK",
                # should it work with "ALK D1203N "? most of our entries just have a protein change in 'name'
                # let's try it without...
                # "name": "ALK D1203N "
                "name": "D1203N",
                "description": "ALK D1203N"
            }
        ]
    }
    normalize_feature_association(fa)
    assert fa == {
        'features': [{
            'provenance_rule': 'default_feature',
            'end': 29940443,
            'name': 'ALK D1203N ',
            'links': [u'http://reg.genome.network/allele/CA346587606', u'http://reg.genome.network/refseq/RS000050',
                u'http://reg.genome.network/refseq/RS000026', u'http://reg.genome.network/refseq/RS000002',
                u'http://myvariant.info/v1/variant/chr2:g.29940443C>T?assembly=hg19',
                u'http://myvariant.info/v1/variant/chr2:g.29717577C>T?assembly=hg38',
                u'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=309051',
                u'http://reg.genome.network/refseq/RS001597'],
            'provenance': ['http://myvariant.info/v1/query?q=ALK D1203N ',
                'http://reg.genome.network/allele?hgvs=NC_000002.11%3Ag.29940443C%3ET'],
            'start': 29940443,
            'synonyms': [u'LRG_488:g.208990G>A', u'chr2:g.29717577C>T', u'NC_000002.10:g.29793947C>T',
                u'CM000664.1:g.29940443C>T', u'chr2:g.29940443C>T', u'NC_000002.12:g.29717577C>T',
                u'NC_000002.11:g.29940443C>T', u'CM000664.2:g.29717577C>T', u'COSM309051', u'NG_009445.1:g.208990G>A'],
            'referenceName': 'GRCh37',
            'geneSymbol': 'ALK',
            'alt': u'T',
            'ref': u'C',
            'chromosome': '2',
            'description': 'ALK D1203N '
        }]
    }


# noinspection PyPep8Naming
def test_normalize_feature_association_AR_amplification():
    fa = {
        'features': [
            {
                "name": 'AR amplification'
            }
        ]
    }
    normalize_feature_association(fa)
    assert fa == {
        'features': [{
            'alt': None,
            'provenance_rule': 'is_amplification',
            'end': 66950461,
            'name': 'AR amplification',
            'provenance': ['http://mygene.info/v3/query?q=AR&fields=genomic_pos_hg19'],
            'start': 66764465,
            'ref': None,
            'referenceName': 'GRCh37',
            'chromosome': 'X',
            'description': 'AR amplification'
        }]
    }


# noinspection PyPep8Naming
def test_normalize_feature_association_ABL1_BCR_fusion():
    fa = {
        'features': [
            {
                "name": "ABL1-BCR fusion"
            }
        ]
    }
    normalize_feature_association(fa)
    print fa
    assert len(fa['features']) == 2


# FA: march 2nd, 2020 duplicate hgvs.g detected

# noinspection PyPep8Naming
def test_normalize_feature_association_PTEN_L325V():
    fa = {
        'features': [
            {
                "geneSymbol": "PTEN",
                # should it work with "ALK D1203N "? most of our entries just have a protein change in 'name'
                # let's try it without...
                # "name": "ALK D1203N "
                "name": "L325V",
                "description": "PTEN L325V"
            }
        ]
    }
    normalize_feature_association(fa)
    assert fa == {
        'features': [{
            'alt': u'G',
            'biomarker_type': 'nonsense',
            'chromosome': '10',
            'description': 'PTEN L325V',
            'end': 89692970,
            'geneSymbol': 'PTEN',
            'hgvs_c': None,
            'hgvs_g': 'NC_000010.10:g.89692970C>G',
            'hgvs_p': None,
            'links': [u'http://reg.genome.network/refseq/RS000034',
                u'http://reg.genome.network/refseq/RS000010',
                u'http://reg.genome.network/refseq/RS000058',
                u'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=5180',
                u'http://myvariant.info/v1/variant/chr10:g.89692970C>G?assembly=hg19',
                u'http://reg.genome.network/allele/CA377482756',
                u'http://reg.genome.network/refseq/RS000591',
                u'http://myvariant.info/v1/variant/chr10:g.87933213C>G?assembly=hg38'],
            'myvariant_hg19': u'http://myvariant.info/v1/variant/chr10:g.89692970C>G?assembly=hg19',
            'name': 'L325V',
            'provenance': [
                u'http://myvariant.info/v1/query?q=snpeff.ann.genename%3APTEN+AND+snpeff.ann.hgvs_p%3Ap.Leu325Val&fields=hg19%2Cvcf%2Csnpeff%2Cchrom%2Ccadd',
                'http://reg.genome.network/allele?hgvs=NC_000010.10%3Ag.89692970C%3EG'],
            'provenance_rule': 'default_feature',
            'ref': u'C',
            'referenceName': 'GRCh37',
            'start': 89692970,
            'synonyms': [u'LRG_311:g.74775C>G',
                u'chr10:g.87933213C>G',
                u'NC_000010.11:g.87933213C>G',
                u'CM000672.1:g.89692970C>G',
                u'COSM5180',
                u'CM000672.2:g.87933213C>G',
                u'NC_000010.9:g.89682950C>G',
                u'NC_000010.10:g.89692970C>G',
                u'NG_007466.2:g.74775C>G',
                u'chr10:g.89692970C>G']
        }]
    }


# noinspection PyPep8Naming
def test_PTEN_L325V_notequal_PTEN_L152V():
    fa1 = {
        'features': [
            {
                "geneSymbol": "PTEN",
                "name": "L325V",
                "description": "PTEN L325V"
            }
        ]
    }
    fa2 = {
        'features': [
            {
                "geneSymbol": "PTEN",
                "name": "L152V",
                "description": "PTEN L152V"
            }
        ]
    }
    normalize_feature_association(fa1)
    normalize_feature_association(fa2)

    for fa in (fa1, fa2):
        for x in fa['features']:
            del x['name']
            del x['description']
            del x['provenance']

    assert json.dumps(fa1, indent=2) != json.dumps(fa2, indent=2)
