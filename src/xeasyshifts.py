from . import queries



def line(num, shift, atomname, resid):
    print 'resid: ' , resid, atomname, shift, num
    return '{:3} {:7.3f} {:7} {:4} {:3}'.format(num, shift, '0.000', atomname, resid)

def xeasy(atoms):
    count, lines = 1, []
    for (s, atomname, resid) in atoms:
        shift = round(s, 3) if s is not None else '999.000'
        lines.append(line(count, shift, atomname, resid))
        count += 1
    return '\n'.join(lines)

def fill_in_missing_atoms(residue):
    pass

def calculate_shift(shifts):
    num = sum([v for (_, v) in shifts])
    denom = len(shifts)
    return (num * 1.0) / denom

def xeasy_shifts():
    """
    model -> String
    """
    myShifts = []
    resShifts = queries.getResidueShifts()
    for key, val in sorted(resShifts.items(), key=lambda x: x[0]):
        # NEED TO FILL IN MISSING ATOM SHIFTS WITH '999.000', WHICH MEANS NEED TO FIND THE AMINO ACID TYPE AND USE 
        #   IT TO LOOK UP WHICH ATOMTYPES ARE REQUIRED, THEN DO A DIFFERENCE AND FIGURE OUT WHICH ONES AREN'T THERE
        # use something like aatypes = dict(enumerate(queries.getData().molecule.residues, start=1))
        print 'key, val: ', key, val
        for atomtypes, shifts in val.items():
            if atomtypes is None: # a HACK b/c I don't know why there are unassigned shifts in residues
                continue
            for atomtype in atomtypes:
                myShifts.append([calculate_shift(shifts), atomtype, key])
#                print 'my shifts: ', myShifts[-1]
    return xeasy(myShifts)
