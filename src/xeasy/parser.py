from unparse import (position, ConsList, fmap, app, optional,
                     many1, seq, many0, bind, getState, commit,
                     seq2L, not0, error, seq2R, alt)
from . import model


item = position.item
literal, satisfy, not1, string = position.literal, position.satisfy, position.not1, position.string

def _oneOf(cs):
    chars = set(cs)
    return satisfy(lambda x: x in chars)

_dig         = _oneOf('0123456789')
_digit       = fmap(int, _dig)
_int         = fmap(''.join, many1(_dig))
_integer     = fmap(int, _int)
_float       = app(lambda s, d, e, f: float(''.join([s, d, '.', f])), 
                   optional('+', literal('-')),
                   _int, 
                   literal('.'), 
                   _int)

_newline     = _oneOf('\n\r\f')
_space       = _oneOf(' \t')


line1 = app(lambda x, y, z: y, 
            string('# Number of dimensions '), 
            _digit, 
            _newline)


# the dim numbers are ignored; they're assumed
# to be the integers 1..n, in order and without
# any skipping
dim = app(lambda _1, _2, _3, i, _4: i, 
          string('# INAME '), 
          _digit, 
          _space, 
          item, 
          _newline)

def dims(n):
    return seq([dim] * n)


# what's the spec say about leading whitespace -- is it optional or required?
_ws_integer = seq2R(many0(_space), _integer)
_ws_float   = seq2R(many1(_space), _float)
_ws_field   = seq2R(many1(_space), many1(not1(alt(_newline, _space))))

def peakline(n):
    return app(lambda ident, shifts, _fields1, height, _fields2, _newl: (ident, model.Peak(shifts, float(height))), 
               _ws_integer,
               bind(getState, lambda s: commit(('bad peakline or something', s), seq([_ws_float] * n))),
               seq([_ws_field] * 2),
               fmap(''.join, _ws_field),
               many0(_ws_field), 
               _newline)


def xeasyAction(dims, peaks):
    return model.PeakFile.fromSimple(dims, peaks)

endCheck = bind(getState, 
                lambda p: not0(seq2L(item, error(('unparsed input remaining', p)))))

xeasy = seq2L(bind(line1, 
                   lambda n: app(xeasyAction, 
                                 dims(n), 
                                 many0(peakline(n)))),
              endCheck)

def run(parser, inp):
    return parser.parse(ConsList.fromIterable(inp), (1, 1))
