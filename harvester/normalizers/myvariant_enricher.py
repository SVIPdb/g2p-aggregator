import logging
import requests

SOURCES = (
    'cadd', 'dbnsfp', 'exac', 'gnomad_genome'
)


def get_myvariant_info(hg19_url, sources=SOURCES):
    """
    Fetches myvariant.info data for the sources specified in SOURCES
    None.
    :param sources:
    :param hg19_url: the full URL to myvariant.info's hg19-assembly version of this variant
    """
    url = hg19_url
    params = {'fields': ",".join(SOURCES)}
    r = requests.get(url, params=params, timeout=60)
    data = r.json()

    if any(k in data for k in SOURCES):
        return data

    # either no hits were found, or they didn't have what we were looking for
    return None


# noinspection PyTypeChecker
def normalize_feature_association(feature_association):
    for feature in feature_association['features']:
        try:
            if 'myvariant_hg19' not in feature or not feature['myvariant_hg19']:
                continue

            # ensure we have location, enrich can create new features
            retrieved = get_myvariant_info(feature['myvariant_hg19'])
            if retrieved:
                feature['mv_info'] = {k: retrieved.get(k) for k in SOURCES}

        except Exception as e:
            logging.exception('mv_enricher: exception {} feature {}'.format(e, feature))
