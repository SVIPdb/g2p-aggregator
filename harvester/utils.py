def str_or_none(x):
    """
    Simply returns a string representation of x if it's not None, or a literal None if it is.
    :param x: the value to stringify
    :return: str(x) if x is not None, otherwise None
    """
    return str(x) if x is not None else None
