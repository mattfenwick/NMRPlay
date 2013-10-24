from . import queries


header = '''
REMARK Chemical Shift Table for my protein of interest

DATA FIRST_RESID 1

DATA SEQUENCE GGGRDYKDDD DKGTMELELR FFATFREVVG QKSIYWRVDD DATVGDVLRS
DATA SEQUENCE LEAEYDGLAG RLIEDGEVKP HVNVLKNGRE VVHLDGMATA LDDGDAVSVF
DATA SEQUENCE PPVAGG

VARS   RESID RESNAME ATOMNAME SHIFT
FORMAT %4d   %1s     %4s      %8.3f

'''

def line(resid, aatype, atomname, shift):
    print 'atom: ' , resid, aatype, atomname, shift
    return '{:4} {} {:>4} {:8.3f}'.format(resid, aatype, atomname, shift)

def talos(atoms):
    lines = []
    for args in atoms:
        lines.append(line(*args))
    return '\n'.join(lines)

def calculate_shift(shifts):
    num = sum([v for (_, v) in shifts])
    denom = len(shifts)
    return (num * 1.0) / denom

def talos_shifts():
    """
    model -> String
    """
    myShifts = []
    resShifts = queries.getResidueShifts()
    aatypes = dict(enumerate(queries.getData().molecule.residues, start=1))
    for key, val in sorted(resShifts.items(), key=lambda x: x[0]):
        print 'key, val: ', key, val
        aa = aatypes[key]
        for atomtypes, shifts in val.items():
            if atomtypes is None: # a HACK b/c I don't know why there are unassigned shifts in residues
                continue
            for atomtype in atomtypes:
                if not atomtype in set(['HN', 'N', 'CA', 'CB', 'C', 'HA', 'HA2', 'HA3']):
                    continue # QUESTION:  is the carbonyl-oxygen named 'C' or 'CO'?
                myShifts.append([key, aa, atomtype, calculate_shift(shifts)])
#                print 'my shifts: ', myShifts[-1]
    return header + talos(myShifts)
