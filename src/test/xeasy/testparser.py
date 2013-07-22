'''
Created on Apr 29, 2013

@author: mattf
'''
import xeasy.parser as p
import unparse.maybeerror as me
import unparse.conslist as c
import unittest as u
import xeasy.model as mod


m = me.MaybeError
l = c.ConsList.fromIterable

def good(rest, state, result):
    return m.pure({'rest': rest, 'state': state, 'result': result})

def runParser(parser, inp):
    return parser.parse(l(inp), (1, 1))


class TestXEasy(u.TestCase):
    
    def testLine1(self):
        inp = '# Number of dimensions 2\nhello'
        self.assertEqual(good(l(inp[-5:]), (2,1), 2), runParser(p.line1, inp))
    
    def testDims(self):
        inp = '# INAME 1 h\n# INAME 2 N\noops'
        self.assertEqual(good(l(inp[-4:]), (3,1), ['h', 'N']), runParser(p.dims(2), inp))
        self.assertEqual(me.MaybeError.zero, runParser(p.dims(2), inp[:20]))
    
    def testPeakLine(self):
        inp = ' 12 1.2 3.4 a b 1123\noops'
        self.assertEqual(good(l(inp[-4:]), (2,1), (12, mod.Peak([1.2, 3.4], 1123))), runParser(p.peakline(2), inp))
#        self.assertEqual(me.MaybeError.zero, runParser(p.dims(2), inp1[:20]))

    def testPeakLine3d(self):
        inp = '4436   9.559 129.430 175.123 1 T          1.077e+02  0.00e+00 a   0    0    0    0 0\nhi'
        self.assertEqual(good(l(inp[-2:]), (2,1), (4436, mod.Peak([9.559, 129.430, 175.123], 107.7))), runParser(p.peakline(3), inp))

    def testXEeasy(self):
        inp = '''
# Number of dimensions 1
# INAME 1 H
 1 8.795 1 T          0.000e+02  0.00e+00 a   0    0    0    0 0
 2 7.229 1 T          0.000e+02  0.00e+00 a   0    0    0    0 0
 3 6.243 1 T          0.000e+02  0.00e+00 a   0    0    0    0 0
 4 7.781 1 T          0.000e+02  0.00e+00 a   0    0    0    0 0
'''
        self.assertEqual(good(l([]), (7, 1), mod.PeakFile(['H'], {1: mod.Peak([8.795], 0), 2: mod.Peak([7.229], 0),
                                                                  3: mod.Peak([6.243], 0), 4: mod.Peak([7.781], 0)})), runParser(p.xeasy, inp[1:]))
