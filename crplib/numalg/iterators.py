# coding=utf-8

"""
Convenience module holding "custom" iterators for generic data structures
like masked arrays
"""


def iter_consecutive_blocks(positions):
    """ Find consecutive blocks in a series of positions,
    i.e. array indices, and return start:end tuples s.t.
    access to the array returns all values in [start ... end]
    :param positions:
    :return:
    """
    start = positions[0]
    for idx, val in enumerate(positions):
        try:
            if val + 1 == positions[idx+1]:
                continue
            else:
                yield start, positions[idx] + 1
                start = positions[idx + 1]
        except IndexError:
            yield start, positions[idx] + 1
    return
