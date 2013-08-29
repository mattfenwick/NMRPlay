from .algebra import fmap, concatMap, groupBy
import re


def parseData(path='bmrb/avgshifts.txt'):
    shifts = {}
    with open(path, 'r') as infile:
        lines = infile.read().split('\n')
        keys = lines[0].split()
        for line in lines[1:]:
            row = dict(zip(keys, line.split()))
            aa, name, shift = row['aminoacid'], row['atomname'], row['avg']
            if aa not in shifts:
                shifts[aa] = {}
            shifts[aa][name] = float(shift)
    return ShiftStats(shifts)


class ShiftStats(object):
    
    def __init__(self, data):
        self.data = data
    
    def __repr__(self):
        return repr(self.__dict__)
    
    def getAtoms(self):
        keys = []
        for v in self.data.values():
            keys.extend(v.keys())
        return set(keys)
    
    def getByAtom(self, atomname):
        if not atomname in self.getAtoms():
            raise ValueError('invalid atom name: ' + atomname)
        shifts = fmap(lambda x: x[atomname] if x.has_key(atomname) else None, self.data)
        return sorted(shifts.items(), key=lambda x:x[1])
    
    def getByMatch(self, expr):
        """
        Data -> String -> [(Aatype, Atom, Float)]
        search through all the data, and retrieve the shift statistics
        for atoms with matching names
        """
        regex = re.compile('(' + expr + ')') # make it capturing so I can extract the text later
        shifts = []
        for (aa, atoms) in self.data.items():
            for (a, shift) in atoms.items():
                match = regex.match(a)
                if match:
                    shifts.append((aa, a, shift))
        return shifts
    
    def getClosest(self, atomname, val):
        shifts = self.getByAtom(atomname)
        return [(aatype, avg, abs(val - avg)) for (aatype, avg) in shifts if avg is not None]
    
    def getClosestN(self, atomnames, vals):
        diffs = concatMap(lambda x: self.getClosest(*x), zip(atomnames, vals))
        grouped = groupBy(lambda x: x[0], diffs, lambda x: x[2])
        return sorted(grouped.items(), key=lambda x: sum(x[1]))

def getClosest(val, xs):
    """
    this looks dumb
    Number -> [Number] -> [(Number, Number)]
    """
    pairs = [(x, abs(x - val)) for x in xs]
    return sorted(pairs, key=lambda x: x[1])

#def getClosestN(vals, seqs):
#    """
#    [Number] -> [[Number]] -> [[(Number, Number)]]
#    """
#    for (val, xs) in zip(vals, seqs):
    


data = {
  'ALA': {
    'C': 177.72,  'CA': 53.17,  'CB': 19.07,  'H': 8.19,    'HA': 4.25,   'HB': 1.35,    'N': 123.28
  },
  'ARG': {
    'C': 176.38,  'CA': 56.77,  'CB': 30.73,  'CD': 43.16,  'CG': 27.23,  'CZ': 160.51,  'H': 8.26,    'HA': 4.3,     'HB2': 1.79,
    'HB3': 1.76,  'HD2': 3.11,  'HD3': 3.09,  'HE': 7.49,   'HG2': 1.56,  'HG3': 1.54,   'HH11': 6.9,  'HH12': 6.83,  'HH21': 6.79,
    'HH22': 6.77, 'N': 120.74,  'NE': 91.03,  'NH1': 79.27, 'NH2': 79.32
  },
  'ASN': {
    'C': 175.21,  'CA': 53.54,  'CB': 38.7,   'CG': 176.44, 'H': 8.35,    'HA': 4.66,    'HB2': 2.8,   'HB3': 2.75,
    'HD21': 7.34, 'HD22': 7.14, 'N': 118.89,  'ND2': 112.76
  },
  'ASP': {
    'C': 176.39,  'CA': 54.69,  'CB': 40.9,   'CG': 178.42, 'H': 8.3,  'HA': 4.59,  'HB2': 2.72,  'HB3': 2.67,  'HD2': 6.29,  'N': 120.65
  },
  'CYS': {
    'C': 174.85,  'CA': 58.21,  'CB': 32.83,  'H': 8.38,  'HA': 4.68,  'HB2': 3.22,  'HB3': 3.13,  'HG': 2.09,  'N': 120.65
  },
  'GLN': {
    'C': 176.3,   'CA': 56.6,  'CB': 29.21,  'CD': 179.4,  'CG': 33.78,  'H': 8.22,  'HA': 4.27,  'HB2': 2.04,  'HB3': 2.01,  'HE21': 7.22,
    'HE22': 7.04, 'HG2': 2.31, 'HG3': 2.29,  'N': 119.84,  'NE2': 111.86
  },
  'GLU': {
    'C': 176.82,  'CA': 57.33,  'CB': 30.04,  'CD': 181.8,  'CG': 36.12,  'H': 8.33,  'HA': 4.25,  'HB2': 2.02,  'HB3': 1.99,
    'HE2': 7.49,  'HG2': 2.27,  'HG3': 2.24,  'N': 120.63
  },
  'GLY': {
    'C': 173.85,  'CA': 45.38,  'H': 8.33,  'HA2': 3.97,  'HA3': 3.89,  'N': 109.8
  },
  'HIS': {
    'C': 175.15,  'CA': 56.55,  'CB': 30.29,  'CD2': 119.91, 'CE1': 137.27, 'CG': 131.42,  'H': 8.26,     'HA': 4.62,  'HB2': 3.18,
    'HB3': 3.12,  'HD1': 10.5,  'HD2': 7.23,  'HE1': 7.79,   'HE2': 11.55,  'N': 119.68,   'ND1': 193.63, 'NE2': 180.0
  },
  'ILE': {
    'C': 175.82,  'CA': 61.62,  'CB': 38.61,  'CD1': 13.5,  'CG1': 27.71,  'CG2': 17.57,  'H': 8.27,  'HA': 4.17,  'HB': 1.78,
    'HD1': 0.67,  'HG12': 1.26, 'HG13': 1.19, 'HG2': 0.77,  'N': 121.43
  },
  'LEU': {
    'C': 176.97,  'CA': 55.65,  'CB': 42.3,  'CD1': 24.72, 'CD2': 24.1,  'CG': 26.81,  'H': 8.22,  'HA': 4.32,  'HB2': 1.61,
    'HB3': 1.52,  'HD1': 0.75,  'HD2': 0.72, 'HG': 1.5,    'N': 122.01
  },
  'LYS': {
    'C': 176.59,  'CA': 56.95,  'CB': 32.81, 'CD': 28.97,  'CE': 41.89,  'CG': 24.92,  'H': 8.18,   'HA': 4.26,  'HB2': 1.77,
    'HB3': 1.74,  'HD2': 1.61,  'HD3': 1.6,  'HE2': 2.91,  'HE3': 2.9,   'HG2': 1.36,  'HG3': 1.35, 'HZ': 7.3,   'N': 120.99,  'NZ': 48.41
  },
  'MET': {
    'C': 176.18,  'CA': 56.15,  'CB': 33.0,  'CE': 17.19,  'CG': 32.04,  'H': 8.26,  'HA': 4.43,  'HB2': 2.03,  'HB3': 1.99,
    'HE': 1.7,  'HG2': 2.35,  'HG3': 2.32,  'N': 120.01
  },
  'PHE': {
    'C': 175.42,  'CA': 58.12,  'CB': 39.99,  'CD1': 131.35,  'CD2': 131.38,  'CE1': 130.49,  'CE2': 130.57,  'CG': 136.94,
    'CZ': 129.02,  'H': 8.35,  'HA': 4.62,  'HB2': 2.99,  'HB3': 2.93,  'HD1': 7.03,  'HD2': 7.03,  'HE1': 7.05,  'HE2': 7.05,
    'HZ': 6.99,  'N': 120.4
  },
  'PRO': {
    'C': 176.61,  'CA': 63.32,  'CB': 31.88,  'CD': 50.3,  'CG': 27.24,  'HA': 4.39,  'HB2': 2.07,  'HB3': 1.99,  'HD2': 3.63,
    'HD3': 3.59,  'HG2': 1.92,  'HG3': 1.89,  'N': 138.53
  },
  'SER': {
    'C': 174.61,  'CA': 58.72,  'CB': 63.72,  'H': 8.28,  'HA': 4.48,  'HB2': 3.87,  'HB3': 3.84,  'HG': 5.53,  'N': 116.29
  },
  'THR': {
    'C': 174.52,  'CA': 62.26,  'CB': 69.61,  'CG2': 21.59,  'H': 8.24,  'HA': 4.46,  'HB': 4.17,  'HG1': 5.33,  'HG2': 1.14,  'N': 115.52
  },
  'TRP': {
    'C': 176.15,   'CA': 57.68,   'CB': 30.12,  'CD1': 126.31, 'CD2': 128.17, 'CE2': 138.38, 'CE3': 120.23, 'CG': 110.29,  'CH2': 123.62,
    'CZ2': 114.05, 'CZ3': 121.21, 'H': 8.28,    'HA': 4.69,    'HB2': 3.18,   'HB3': 3.12,   'HD1': 7.13,   'HE1': 10.07,  'HE3': 7.28,  
    'HH2': 6.94,   'HZ2': 7.26,   'HZ3': 6.83,   'N': 121.65,  'NE1': 129.25
  },
  'TYR': {
    'C': 175.37, 'CA': 58.13,  'CB': 39.33,  'CD1': 132.43, 'CD2': 132.43, 'CE1': 117.71, 'CE2': 117.78, 'CG': 128.11,  'CZ': 154.67,
    'H': 8.31,   'HA': 4.62,   'HB2': 2.89,  'HB3': 2.83,   'HD1': 6.91,   'HD2': 6.91,   'HE1': 6.69,   'HE2': 6.69,   'HH': 9.17,   'N': 120.95
  },
  'VAL': {
    'C': 175.63,  'CA': 62.53,  'CB': 32.75,  'CG1': 21.53,  'CG2': 21.33,  'H': 8.28,  'HA': 4.18,  'HB': 1.98,  'HG1': 0.82,
    'HG2': 0.8,   'N': 121.24
  }}


stats = ShiftStats(data)
