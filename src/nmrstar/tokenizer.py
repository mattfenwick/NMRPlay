from parse.standard import Parser
from tokens import Token


# token regexes
#
# newline = "^[\n\r\f]"
# blank = "^[ \t]"
# space = 

# tokens
#   needs error reporting
#   and take care of return values

def lit(c):
    return Parser.satisfy(lambda t: t.char == c)

def str(cs):
    return Parser.all(map(lit, cs))

def extract(cs):
    return ''.join([x.char for x in cs])


newline = Parser.any(map(lit, '\n\r\f'))

blank = lit(' ').plus(lit('\t'))

space = blank.plus(newline)

special = Parser.any(map(lit, '"#_')).plus(space)

stop = str("stop_").fmap(lambda xs: Token('stop', extract(xs), xs[0].meta))

saveclose = str("save_").fmap(lambda xs: Token('saveclose', extract(xs), xs[0].meta))

loop = str("loop_").fmap(lambda xs: Token('loop', extract(xs), xs[0].meta))


def commentAction(pd, chars):
    return Token('comment', extract(chars), pd.meta)
    
comment = Parser.app(commentAction, lit('#'), newline.not1().many0())


def dataAction(open, name):
    return Token('dataopen', extract(name), open[0].meta)

dataopen = Parser.app(dataAction, str("data_"), space.not1().many1())
    
    
def saveAction(open, name):
    return Token('saveopen', extract(name), open[0].meta)

saveopen = Parser.app(saveAction, str("save_"), space.not1().many1())


def idAction(und, name):
    return Token('identifier', extract(name), und.meta)

identifier = Parser.app(idAction, lit('_'), space.not1().many1())


def uqAction(c, cs):
    return Token('value', extract([c] + cs), c.meta)
    
unquoted = Parser.app(uqAction, special.not1(), space.not1().many0())

sc, sq, dq = map(lit, ';\'"')

endsc = newline.seq2R(sc)

def scAction(open, body, close):
    return Token('value', extract(body), open.meta)
    
scstring = Parser.app(scAction, sc, endsc.not1().many0(), endsc)

def sqRest(open):
    def action(cs):
        return Token('value', extract(cs), open.meta)
    return sq.not1().many1().fmap(action).seq2L(sq)

sqstring = sq.bind(sqRest)

def dqRest(open):
    def action(cs):
        return Token('value', extract(cs), open.meta)
    return dq.not1().many1().fmap(action).seq2L(dq)

dqstring = dq.bind(dqRest)

value = Parser.any([sqstring, dqstring, scstring, unquoted])


whitespace = blank.many1().fmap(lambda xs: Token('whitespace', extract(xs), xs[0].meta))

newlines = newline.many1().fmap(lambda xs: Token('newlines', extract(xs), xs[0].meta))

token = Parser.any([dataopen, saveopen, saveclose,
                    loop, stop, value, whitespace,
                    newlines, comment, identifier])
                    
def noParse(ts):
    if len(ts) == 0:
        return Parser.zero
    else:
        return token.commit({'message': 'unable to parse token', 'position': ts[0].meta})
        
scanner = Parser.get.bind(noParse).many0()

# scanner = token.many0().seq2L(Parser.get(lambda xs: if len(xs) == 0 then Parser.pure(None) else 

