import itertools


def pairwise(iterable):
    #~ "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def is_within(xa, xb):
    if xa[0] >= xb[0] and xa[1] <= xb[1]:
        return True
    else:
        return False
