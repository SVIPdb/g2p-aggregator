
# normalized vocabulary for evidence_direction
# 'resistant', 'sensitive', or 'no benefit'


def evidence_direction(evidence, association, na=False):
    resistant = ['resistant', 'resistance', 'poor outcome',
                 'decreased response']
    sensitive = ['sensitive', 'predictive of response']
    nb = ['no benefit']

    res_type = {
        'resistant': resistant,
        'sensitive': sensitive,
        'no benefit': nb
    }

    for item in res_type:
        for opt in res_type[item]:
            if evidence and opt in evidence.lower():
                association['response_type'] = item

    if 'response_type' not in association:
        if na:
            association['response_type'] = 'NA'
        else:
            association['response_type'] = evidence

    return association


def evidence_direction_biological(evidence, association, na=False):
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
            association['response_type'] = code
            break
    else:
        association['response_type'] = 'NA' if na else None

    return association
