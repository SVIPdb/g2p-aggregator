
# normalized vocabulary for evidence_direction
# 'resistant', 'sensitive', or 'no benefit'


# FIXME: 'evidence direction' may be a misnomer, since this value is always assigned to an association's
#  'response_type' field

def evidence_direction(evidence, na=False):
    res_type = {
        'resistant': ['resistant', 'resistance', 'poor outcome', 'decreased response'],
        'sensitive': ['sensitive', 'predictive of response'],
        'no benefit': ['no benefit']
    }

    for item in res_type:
        for opt in res_type[item]:
            if evidence and opt in evidence.lower():
                return item

    # if we reach here, it didn't match anything
    return 'NA' if na else evidence


def evidence_direction_biological(evidence, na=False):
    patho_mapping = {
        ('oncogenic',):
            'Pathogenic',
        ('likely oncogenic',):
            'Likely Pathogenic',
        ('neutral',):
            'Benign',
        ('likely neutral',):
            'Likely Benign',
        ('inconclusive',):
            'Uncertain Significance',
    }

    for patterns, code in patho_mapping.items():
        if evidence.lower() in patterns:
            return code

    return 'NA' if na else None
