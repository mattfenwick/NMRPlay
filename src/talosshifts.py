

def line(num, shift, atomname, resid):
    return ' '.join([
        str(num), 
        str(shift),
        '0.000',
        atomname,
        str(resid)
    ])

def talos(atoms):
    count, lines = 1, []
    for (s, atomname, resid) in atoms:
        shift = s if s is not None else '999.000',
        lines.append(line(count, shift, atomname, resid))
    return '\n'.join(lines)

def fill_in_missing_atoms(residue):
    pass


