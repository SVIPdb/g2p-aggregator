
def filter_feature_association(feature_association):
    """
    Only retain feature associations that have a HGVS.c string. This filters out
    "weird" variants like fusions, frameshifts, etc. for which we don't do a good job
    of inferring the HGVS string.

    :param feature_association: the feature association to filter
    :return: True if the association has an HGVS.c string, False otherwise
    """
    try:
        if not feature_association['features'][0]['hgvs_g']:
            return False
    except (KeyError, IndexError):
        return False

    return True
