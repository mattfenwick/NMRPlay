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

def talos_shifts():
    """
    () -> String
    """
    myShifts = []
    shifts = queries.getShifts()
    aatypes = dict(enumerate(queries.getData().molecule.residues, start=1))
    for (resid, atom, shift) in shifts:
        aa = aatypes[resid]
        if not atom in set(['H', 'N', 'CA', 'CB', 'C', 'HA', 'HA2', 'HA3', 'QA']):
            continue # QUESTION:  is the carbonyl-oxygen named 'C' or 'CO'?
        elif atom == 'H':
            myShifts.append([resid, aa, 'HN', shift])
        elif atom == 'QA':
            myShifts.append([resid, aa, 'HA2', shift])
            myShifts.append([resid, aa, 'HA3', shift])
        else:
            myShifts.append([resid, aa, atom, shift])
#                print 'my shifts: ', myShifts[-1]
    return header + talos(myShifts)
