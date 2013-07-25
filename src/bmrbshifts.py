from .algebra import fmap


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
