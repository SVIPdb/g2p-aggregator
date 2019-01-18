import requests
import urllib
import logging
import re


def reference_names(n):
    return {
        '37': 'GRCh37',
        37: 'GRCh37',
        'Hg37': 'GRCh37',
        'GRCh37/hg19': 'GRCh37',
        'grch37_hg19': 'GRCh37',
        'Hg38': 'GRCh38',
        '38': 'GRCh38',
        38: 'GRCh38'
    }.get(n, n)


def normalize(name):
    try:
        name = name.encode('utf8')
    except Exception as e:
        pass

    return reference_names(name)


def normalize_feature_association(feature_association):
    """ given the 'final' g2p feature_association,
    update it with a standard reference name """
    # nothing to read?, return
    if 'features' not in feature_association:
        return

    for feature in feature_association['features']:
        if 'referenceName' in feature:
            feature['referenceName'] = normalize(feature['referenceName'])
