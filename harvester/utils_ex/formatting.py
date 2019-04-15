import json


def unicode_or_none(x):
    """
    Simply returns a string representation of x if it's not None, or a literal None if it is.
    :param x: the value to stringify
    :return: str(x) if x is not None, otherwise None
    """
    return unicode(x) if x is not None else None


def jsonify_or_none(x):
    return json.dumps(x) if x is not None else None


def capitalize_words(x):
    if not x:
        return x
    return " ".join(tok[0].upper() + tok[1:].lower() for tok in x.split(" "))


