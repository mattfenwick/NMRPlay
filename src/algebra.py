import itertools


def product(xs, ys):
    return list(itertools.product(xs, ys))


def inner_join(predicate, xs, ys):
    """
    (a -> b -> Bool) -> [a] -> [b] -> [(a, b)]
    """
    return [pair for pair in product(xs, ys) if predicate(*pair)]


def ffilter(predicate, coll):
    if isinstance(coll, list):
        return filter(predicate, coll)
    elif isinstance(coll, dict):
        out = {}
        for (key, value) in dict.iteritems():
            if predicate(value):
                out[key] = value
        return out
    raise TypeError(('invalid type to ffilter', coll))


def fmap(f, coll):
    if isinstance(coll, list):
        return map(f, coll)
    elif isinstance(coll, dict):
        out = {}
        for (key, value) in dict.iteritems():
            out[key] = f(value)
        return out
    raise TypeError(('invalid type to fmap', coll))


def f_id(x):
    return x


# can't use itertools.groupby, b/c it:
#  1) assumes its input is sorted
#  2) output shares generator or something, making results very delicate to access
def groupBy(fkey, xs, fval=f_id):
    """
    (a -> b) -> [a] -> Maybe (a -> c) -> Map b [c]
    """
    out = {}
    for x in xs:
        v = fkey(x)
        if not out.has_key(v):
            out[v] = []
        out[v].append(fval(x))
    return out
