


def xez_to_star(paths, star):
    q = xeasy_in('mattsnmrproject', paths)
    star_out(q, star)

def star_to_xez(paths, star):
    q = star_in(star)
    xeasy_out(q, paths)


# some more examples
def star_filter_hsqc_to_xez(nhsqc, star, tag):
    q = star_in(star)
    new_nhsqc = ut.filterSpecPeaks(ut.peakHasTag(tag), q.spectra['nhsqc'])
    q.spectra['nhsqc'] = new_nhsqc
    xeasy_out(q, {'nhsqc': nhsqc})

def star_in(starpath):
    with open(starpath, 'r') as starfile:
        star = nsp.fullParse(starfile.read())
    if star.status != 'success':
        raise ValueError(('unable to parse star project file', star))
    return nc.star2patch(star.value['result'])

def star_out(pmodel, starpath):
    raise ImportError('missing unparsing module for nmr-star serialization')
    text = unp.unparse(nc.patch2star(pmodel))
    with open(starpath, 'w') as outfile:
        outfile.write(text)
    return None

