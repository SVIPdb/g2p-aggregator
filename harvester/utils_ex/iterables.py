"""
Contains miscellaneous iterable-manipulation methods.
"""


def try_flatten(seq, scalar_classes=(str,)):
    """
    Performs one level of flatting 'seq', but returns non-interables as-is.

    :param seq: the input sequence to flatten
    :param scalar_classes: anything that is in this list is treated as a scalar and not unpacked

    >>> list(try_flatten([1,2,3,4]))
    [1, 2, 3, 4]
    >>> list(try_flatten([1,2,[3,4],5]))
    [1, 2, 3, 4, 5]
    >>> list(try_flatten([1,2,[3,4,[5]],6]))  # only performs one level of un-nesting
    [1, 2, 3, 4, [5], 6]
    >>> list(try_flatten(["hello", "there"]))  # don't flatten iterables in scalar_classes (i.e., strings by default)
    ['hello', 'there']
    >>> list(try_flatten(["hello", "there", ["friend"]]))  # still flatten lists of scalar_classes, though
    ['hello', 'there', 'friend']
    """
    for item in seq:
        if any(isinstance(item, c) for c in scalar_classes):
            yield item
        else:
            try:
                for subitem in item:
                    yield subitem
            except TypeError:
                yield item
