from unparse import combinators as c
from unparse import conslist
from . import model


item = c.itemPosition
literal, satisfy, not1, string = c.tokenPosition

def _oneOf(cs):
    chars = set(cs)
    return satisfy(lambda x: x in chars)

_dig         = _oneOf('0123456789')
_digit       = c.fmap(int, _dig)
_int         = c.fmap(''.join, c.many1(_dig))
_integer     = c.fmap(int, _int)
_float       = c.app(lambda s, d, e, f: float(''.join([s, d, '.', f])), 
                     c.optional('+', literal('-')),
                     _int, 
                     literal('.'), 
                     _int)

_newline     = _oneOf('\n\r\f')
_space       = _oneOf(' \t')


line1 = c.app(lambda x, y, z: y, 
              string('# Number of dimensions '), 
              _digit, 
              _newline)


# the dim numbers are ignored; they're assumed
# to be the integers 1..n, in order and without
# any skipping
dim = c.app(lambda _1, _2, _3, i, _4: i, 
            string('# INAME '), 
            _digit, 
            _space, 
            item, 
            _newline)

def dims(n):
    return c.all_([dim] * n)


# what's the spec say about leading whitespace -- is it optional or required?
_ws_integer = c.seq2R(c.many0(_space), _integer)
_ws_float   = c.seq2R(c.many1(_space), _float)
_ws_field   = c.seq2R(c.many1(_space), c.many1(not1(c.plus(_newline, _space))))

def peakline(n):
    return c.app(lambda ident, shifts, _fields1, height, _fields2, _newl: (ident, model.Peak(shifts, float(height))), 
                 _ws_integer,
                 c.bind(c.getState, lambda s: c.commit(('bad peakline or something', s), c.all_([_ws_float] * n))),
                 c.all_([_ws_field] * 2),
                 c.fmap(''.join, _ws_field),
                 c.many0(_ws_field), 
                 _newline)


def xeasyAction(dims, peaks):
    return model.PeakFile.fromSimple(dims, peaks)

endCheck = c.bind(c.getState, 
                  lambda p: c.not0(c.seq2L(item, c.error(('unparsed input remaining', p)))))

xeasy = c.seq2L(c.bind(line1, 
                       lambda n: c.app(xeasyAction, 
                                       dims(n), 
                                       c.many0(peakline(n)))),
                endCheck)

def run(parser, inp):
    return parser.parse(conslist.ConsList.fromIterable(inp), (1, 1))
