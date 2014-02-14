from . import inout as pt
from . import model
from . import querymodel
from .bmrbshifts import stats
from .algebra import ffilter, inner_join, groupBy, split, concatMap, fmap
from operator import attrgetter


ROOT = '../PPAbnormal/'

def getData(root=ROOT, simple=True):
    data = pt.json_in(root + 'project.txt')
    if simple:
        return querymodel.fromModel(data)
    return data

def fst(x):
    return x[0]

def snd(x):
    return x[1]

def filterSpecPeaks(pred, spec):
    '''
    Filter peaks from a spectrum based on a predicate,
    returning a new spectrum.
    '''
    peaks = ffilter(pred, spec.peaks)
    return model.Spectrum(spec.axes, peaks)


def getAllPeaks(prj):
    """
    Project -> [Peak]
    """
    peaks = []
    for spectrum in prj.getSpectra():
        peaks.extend(spectrum.getPeaks())
    return peaks


def joinSSToPeaks(prj):
    """
    Project -> [(SpinSystem, [Peak])]
    """
    joined, spins = [], prj.getSpinSystems()
    for ss in spins:
        pks = [prj._spectra[spec]._peaks[pk_id] for (spec, pk_id) in ss.pkids]
        joined.append((ss, pks))
    return joined

def joined2(spins, peaks):
    """
    [SpinSystem] -> [Peak] -> [(SpinSystem, Peak)]
    
    this is nearly identical to 'joinSSToPeaks' except for the lack of grouping
    """
    def p(ss, pk):
        return [pk.spectrum_name, pk.id] in ss.pkids
    joined = inner_join(p, spins, peaks)
    return joined


def findSpinsystemsOfPeak(spins, peak):
    """
    [SpinSystem] -> Peak -> [Int]
    """
    return filter(lambda ss: [peak.spectrum_name, peak.id] in ss.pkids, spins)


def findCloseNHSQCPeaks(proj, HTOL=0.025, NTOL=0.2):
    """
    Project -> [((Peak, [SpinSystem]), (Peak, [SpinSystem]))]
    """
    nhsqc = proj._spectra['nhsqc']
    spins = proj.getSpinSystems()
    
    if nhsqc.axes != ['N', 'H']:
        raise ValueError('oops, unexpected axis order!')
    
    peaks = nhsqc.getPeaks()
    
    for p1 in peaks:
        for p2 in peaks:
            n1, h1 = map(fst, p1.dims)
            n2, h2 = map(fst, p2.dims)
            if abs(n2 - n1) <= NTOL and abs(h2 - h1) <= HTOL and p1.id < p2.id:
                print 'peak1', p1
                print 'peak2', p2
                print 'spin systems1: ', findSpinsystemsOfPeak(spins, p1)
                print 'spin systems2: ', findSpinsystemsOfPeak(spins, p2)
                print '\n'


def analyzeSpinSystems(proj):
    """
    Project -> ([(SpinSystemID, [Float], [(Int, String)])], Map Int Int)
    """
    hncopeaks = proj._spectra['hnco'].getPeaks()
    spins = proj.getSpinSystems()
    
    # 1. find lone NHSQC peaks
    #   - get list of all NHSQC peak ids -- wait, isn't this unnecessary since I have a spin system for each NHSQC peak?
    #   - figure out which spin systems don't have any hnco peaks
    print 'spin systems with more or less than 1 HNCO peak:'
    joined = joined2(spins, hncopeaks)
    grouped = groupBy(fst, joined, snd)
    badss = filter(lambda x: len(x[1]) != 1, grouped)
    for ss in badss:
        print ss
    print '\n'
    
    # 2. find ??? HNCO peaks -- member of < or > 1 spin system
    #   - get list of HNCO peak ids
    #   - look through all spin systems, and total up the hnco peak ids in each
    backbonepks = joined2(spins, filter(lambda pk: 'backbone amide' in pk.tags, hncopeaks))
    grped2 = groupBy(snd, backbonepks, fst)
    badpks = filter(lambda q: len(q[1]) != 1, grped2)
    for pk in badpks:
        print pk


def checkTags(proj):
    """
    Project -> [(String, PeakID)]
    """
    print 'find "problems" with HNCO peak tags'
    for pk in proj._spectra['hnco'].getPeaks():
        if len(pk.tags) != 1:
            print 'problem'
        else:
            print 'good', pk


def findJunkPeaksInSpinSystems(proj):
    """
    Project -> ???
    """
    hnco = proj._spectra['hnco']
    pks = filter(lambda pk: pk.tags != ['backbone amide'], hnco.getPeaks())
    uhohs = joined2(proj.getSpinSystems(), pks)
    for (ss, pk) in uhohs:
        print 'uh-oh, problem with ', pk.id, 'in spin system', ss.id, 'with tags', pk.tags

    # find the spin systems for each hnco peak ???
    #   I don't understand this
    pk_ss = joined2(proj.getSpinSystems(), hnco.getPeaks())
    grped2 = groupBy(snd, pk_ss, fst)
    for (peak, sss) in grped2.items():
        print peak.tags, '    ', sss, '    ', peak.getPosition()


def findHNCACBStrangeness(proj):
    """
    Project -> ([Int], [(Int, [Int])])
    """
    hncacb = proj.spectra['hncacb']
    spins = proj.spinsystems
    
    hncacb_peaks = model.fmap_dict(lambda _: [], hncacb.peaks)
    
    for (ssid, ss) in spins.items():
        cnt = 0
        for (spec, pkid) in ss.pkids:
            if spec == 'hncacb':
                cnt += 1
                hncacb_peaks[pkid].append(ssid)
        if cnt == 0:
            print 'spin system', ssid, 'has 0 hncacb peaks'
    
    for (pid, ssids) in hncacb_peaks.items():
        if hncacb.peaks[pid].tags != []:
            continue
        if len(ssids) != 1:
            print 'hncacb peak', pid, 'in spin systems:', ssids
        else:
            print pid, 'is good to go' 
    # HNCACB peaks in more or less than 1 spin system
    # spin systems with 0 HNCACB peaks
    # spin systems with lots of HNCACB peaks?  what's a good threshold?


def inspectHNCACBPeaks(proj):
    """ look at whether they're tagged or not"""
    pks = proj._spectra['hncacb'].getPeaks()
    (good, bad) = split(lambda pk: pk.tags == [], pks)
    
    for i, p in enumerate(good):
        if p.height > 0:
            print 'C-Alpha? ', p, i
        elif p.height < 0:
            print 'C-Beta? ', p, i
        else:
            raise ValueError('peak with 0 height ???')
    
    print '\n\n'
    
    for j, pk in enumerate(bad):
        print 'bad peak', pk, j


def findHNCACBPeaksInSSAndTaggedBackbone():
    """
    HNCACB peaks in a spin system but not tagged 'backbone' or vice versa
    """
    proj = getData()
    # map (get 1) . filter ((==) "hncacb" . get 0) $ concatMap (get pkids) (getSpinSystems proj)
    inss = map(lambda x: x[1], filter(lambda n: n[0] == 'hncacb', concatMap(lambda ss: ss.pkids, proj.getSpinSystems())))
    back = map(lambda x: x.id, filter(lambda pk: pk.tags == ['backbone'], proj._spectra['hncacb'].getPeaks()))
    return inss, back


def hncacbWeirdness():
    """hncacb peaks not tagged 'backbone' but in a spin system,
       or those tagged 'backbone' but not in a spin system"""
    proj = getData()
    peaks = proj._spectra['hncacb'].getPeaks()
    inss = set(map(lambda x: x[1].id, joined2(proj.getSpinSystems(), peaks)))
    ssid_to_pkids = groupBy(lambda x:x[0].id, joined2(proj.getSpinSystems(), peaks), lambda y: y[1].id)
    pkid_to_ssids = {}
    for (s, ps) in ssid_to_pkids.items():
        for p in ps:
            if not pkid_to_ssids.has_key(p):
                pkid_to_ssids[p] = []
            pkid_to_ssids[p].append(s)
    counts = ([], [], [], [])
    tags = set([])
    for peak in peaks:
        tags = tags.union(set(peak.tags))
        if 'edge' in peak.tags or 'processing artifact' in peak.tags:
            if peak.id in inss:
                counts[2].append(peak.id)
                print 'badly tagged but in spin system', peak.id, peak.tags, ' spins: ', pkid_to_ssids[peak.id]
            else:
                counts[3].append(peak.id)
                pass # good
        else:
            if peak.id in inss:
                counts[0].append(peak.id)
            else:
                counts[1].append(peak.id)
                print (' ' * 75) + 'nicely tagged but not in spin system', peak.id, peak.tags
    print map(len, counts)
    print 'tags: ', tags
    return counts
#    print counts[2]
#    proj2 = getData(simple=False)
#    bads = map(lambda q: ['hncacb', q], counts[2])
#    for (_, ss) in proj2.spinsystems.iteritems():
#        if ['hncacb', 769] in ss.pkids:
#            print 'wow', ss
#        ss.pkids = filter(lambda x: x not in bads, ss.pkids)
#    pt.json_out(ROOT + "project8.txt", proj2)


def hncacbOverlap():
    proj = getData()
#    jed = joined2(proj.getSpinSystems(), getAllPeaks(proj))
#    return groupBy(fst, jed, snd)
#    return jed
    peaks = proj._spectra['hncacb'].getPeaks()
    simple = [(i, x.dims[2][0], x.dims[2][1][0]) for (i, x) in enumerate(peaks, start=1) if 'backbone' in x.tags]
    selfjoined = inner_join(lambda p, q: p[2] == q[2] and abs(p[1] - q[1]) < 0.01 and p[0] < q[0], simple, simple)
    return selfjoined


def getCACBPairs():
    proj = getData()
    goodss = [s for s in proj.getSpinSystems() if 'backbone' in s.tags]
    jed = joined2(goodss, proj._spectra['hncacb'].getPeaks())
    gred = groupBy(lambda x: x[0].id, jed, lambda y: y[1].getPosition()) # RIGHT HERE !!!
    return gred


def findAlaPeaks():
    proj = getData()
    pks = []
    for x in proj._spectra['hncacb'].getPeaks():
        if x.getPosition()[2] < 25 and 'backbone' in x.tags:
            pks.append(x)
    return sorted(pks, key=lambda x: x.getPosition()[2])


def findAlaSS():
    proj = getData()
    return sorted([s for s in proj.getSpinSystems() if 'alanine' in s.aatypes], key=lambda ss: ss.id)


def findSSByAAType(aatype='alanine'):
    proj = getData()
    return [ss for ss in proj.getSpinSystems() if aatype in ss.aatypes]


def findSerThrPeaks():
    proj = getData()
    pks = []
    for x in proj._spectra['hncacb'].getPeaks():
        if x.getPosition()[2] > 50 and 'backbone' in x.tags and 'CB' in x.dims[2][1]:
            pks.append(x)
    return sorted(pks, key=lambda x: x.getPosition()[2])


def findVIPeaks():
    """not 100% sure that they'll be Valine/Isoleucine CAs"""
    proj = getData()
    pks = []
    for x in proj._spectra['hncacb'].getPeaks():
        if x.getPosition()[2] > 60 and 'backbone' in x.tags and 'CA' in x.dims[2][1] and not 'i' in x.tags and not 'i-1' in x.tags:
            pks.append(x)
    return sorted(pks, key=lambda x: x.getPosition()[2])


def findGlyPeaks(upper=50):
    proj = getData()
    pks = []
    for x in proj._spectra['hncacb'].getPeaks():
        if x.getPosition()[2] < upper and 'CA' in x.dims[2][1] and 'backbone' in x.tags:
            pks.append(x)
    return sorted(pks, key=lambda x: x.getPosition()[2])


def findAssignedPeaks(specname='hncacb'):
    proj = getData()
    pks, junk = [], []
    for pk in proj._spectra[specname].getPeaks():
        if 'backbone' in pk.tags:
            if 'i' in pk.tags or 'i-1' in pk.tags:
                pks.append(pk)
            else:
                junk.append(pk)
        else:
            junk.append(pk)
    return (pks, junk)


def findSSFragments():
    proj = getData()
    spins = proj.getSpinSystems()
    frags = {}
    for s in spins:
        if not 'backbone' in s.tags: # not sure how important it is to ditch these ...
            if len(s.ssnexts) > 0: # or len(s.residueids) > 0:  # <-- this is now okay since I assigned a sidechain system
                raise ValueError('something screwy going on here ...')
            continue                 # after all, they don't hurt anything, right? wait, yes they do:  they make it confusing
        if len(s.ssnexts) > 1:
            raise ValueError('oops, not expecting ambiguous nexts')
        frags[s.id] = s.ssnexts[:] # make a new copy ... just b/c I don't like modding data that I don't "own"
    while True:
        changed = False
        for k in frags:
            if len(frags[k]) > 0 and frags[k][-1] in frags:
                # print 'k: ', k, frags
                temp = frags.pop(frags[k][-1])
                frags[k].extend(temp)
                changed = True
                break # umm ... which loop does this break out of?
        if not changed:
            break
    return frags


def findMatchingCACB(ca, cb, tola=0.2, tolb=0.2):
    proj = getData()
    cas, cbs = [], []
    joined = groupBy(lambda x: x[1].id, 
                     joined2([s for s in proj.getSpinSystems() if 'backbone' in s.tags], proj._spectra['hncacb'].getPeaks()), 
                     lambda y: y[0].id)
    for pk in proj._spectra['hncacb'].getPeaks():
        shift, tags = pk.dims[2]
 #       if 'CA' in tags:
        if pk.height > 0: # it's a CA ... usually !!! 
            adiff = abs(shift - ca)
            if adiff <= tola:
                cas.append((adiff, pk, joined[pk.id] if joined.has_key(pk.id) else None))
#        elif 'CB' in tags:
        elif pk.height < 0: # it's a CB ... usually !!!
            bdiff = abs(shift - cb)
            if bdiff <= tolb:
#                print 'found one: ', bdiff, shift, cb
                cbs.append((bdiff, pk, joined[pk.id] if joined.has_key(pk.id) else None))
        else: # nothing to do
            pass
    both = []
    for (_1, a, ss1) in cas:
        for (_2, b, ss2) in cbs:
            if abs(a.getPosition()[0] - b.getPosition()[0]) <= 0.04 and abs(a.getPosition()[1] - b.getPosition()[1]) <= 0.2:
                both.append(((_1, _2), (a, b), (ss1, ss2)))
    out = (sorted(cas, key=fst), 
           sorted(cbs, key=fst), 
           sorted(both, key=lambda x: sum(x[0])))
    return out


def findSequenceAssignments():
    proj = getData()
    res = dict((i, []) for (i, _) in enumerate(proj.molecule.residues, start=1))
    for s in proj.getSpinSystems():
        for i in s.residueids:
            res[i].append(s.id)
    for (k, v) in sorted(res.items(), key=fst):
        aatype = proj.molecule.residues[k - 1]
        if len(v) == 1:
            print k, aatype, '  ', v[0]
        elif len(v) == 0:
            print k, aatype, '   --'
        else:
            print k, aatype, '  ', v
#    return res


def shiftAnalyzer(resids=range(1, 106 + 1)):
    shifts = getResidueShifts()
    sss = getData().getSpinSystems()
    res2ss = {}
    for ss in sss:
        if len(ss.residueids) == 1:
            res2ss[ss.residueids[0]] = ss.id
        elif len(ss.residueids) > 1:
            raise ValueError('oops')
    for x in resids:
        if res2ss.has_key(x):
            print 'peaks for residue', x
            ssid = res2ss[x]
            for l in getSSPeaks(ssid):
                if True: # l[0] == 'cconh':
                    print l
        else:
            print 'no peaks for residue', x
        _d = raw_input('press enter to continue\n')
        print
        print len(shifts[x]), 'assigned shifts for residue', x
        for m in sorted(shifts[x].keys()):
            print m, shifts[x][m]
        print '\n'


def checkAtomtypes():
    proj = getData()
    for ss in proj.getSpinSystems():
        if not 'backbone' in ss.tags:
            print 'skipping -- not backbone -- ', ss.id, ss.tags
            continue
        if len(ss.residueids) == 0:
            print 'skipping -- no residue assignment -- ', ss.id
            continue
        for (specname, pkid) in ss.pkids:
            pk = proj._spectra[specname]._peaks[pkid]
            for (_, ats) in pk.dims:
                if len(ats) == 0:
                    print ('poopy missing atomtype assignments -- %s %s %s' % (ss.tags, ss.id, pk))
                elif len(ats) > 1:
                    print 'lots!: ', ss.id, pk


def removeIndexTags():
    proj = getData(simple=False)
    for ssid in sorted(proj.spinsystems):
        ss = proj.spinsystems[ssid]
        for (specname, pkid) in ss.pkids:
            pk = proj.spectra[specname].peaks[pkid]
            if 'i' in pk.tags: pk.tags.remove('i')
            if 'i-1' in pk.tags: pk.tags.remove('i-1')
    pt.json_out(ROOT + "project2.txt", proj)
    


def getResidueShifts():
    proj = getData()
    shifts = dict([(i, {}) for i in range(1, 106 + 1)])
    def add(n, a, sh):
        if not shifts[n].has_key(a):
            shifts[n][a] = []
        shifts[n][a].append(sh)
    for ss in proj.getSpinSystems():
        if not 'backbone' in ss.tags:
            print 'skipping -- not backbone -- ', ss.id, ss.tags
            continue
        if len(ss.residueids) == 0:
            print 'skipping -- no residue assignment -- ', ss.id
            continue
        resid = ss.residueids[0]
        for (specname, pkid) in ss.pkids:
            pk = proj._spectra[specname]._peaks[pkid]
            for (sh, ats) in pk.dims:
                for a in ats:
                    datum = (specname, pk.id, sh) 
                    if a[-5:] == '(i-1)':
                        add(resid - 1, a[:-5], datum)
                    else:
                        add(resid, a, datum)
    return shifts


def avg(vals):
    return sum(vals) * 1.0 / len(vals)


def selectOneShift(shifts):
    """
    if an atom shows up in multiple peaks, we want to combine those
    multiple estimates of chemical shift into a single value
    some spectra are better than others, e.g. HCCH-Tocsy > HCCONH-Tocsy
    """
    # use all hcchtocsy peaks if there's any,
    # otherwise use all peaks
    hcchtocsy = filter(lambda x: x[0] == 'hcchtocsy', shifts)
    if len(hcchtocsy) > 0:
        return avg([v for (_1, _2, v) in hcchtocsy])
    return avg([u for (_3, _4, u) in shifts])


def removeDuplicateAtoms(residue):
    """
    if we have both ambiguous-style -- QB -- and unambiguous-style
    -- HB[2]/HB[3] -- assignments for the same group of atoms, get
    rid of the QB shift b/c it has less information
    """
    dupes = {
        'QA' : ['HA2' , 'HA3' ],
        'QB' : ['HB2' , 'HB3' ],
        'QG' : ['HG2' , 'HG3' ],
        'QG1': ['HG12', 'HG13'],
        'QQG': ['QG1' , 'QG2' ],
        'QD' : ['HD2' , 'HD3' ],
        'QQD': ['QD1' , 'QD2' ]
    }
    new_residue = {}
    for (atom, data) in residue.items():
        c1 = dupes.has_key(atom)
        if c1:
            c2 = all([k in residue for k in dupes[atom]])
#            print atom, c1, c2
            if c2:
                continue
        new_residue[atom] = data
    return new_residue


def disambiguateAtom(atomname):
    """
    cyana doesn't care which way a stereospecific assignment goes,
    just that there are two possibilities and it'll pick one.
    so get rid of the ambiguity markers: e.g. `HB[2]` -> `HB2`
    """
    return filter(lambda x: x not in '[]', atomname)


def checkShifts():
    shifts = getResidueShifts()
    best = {}
    for (resid, cs) in shifts.items():
        residue = {}
        for (atom, shs) in cs.items():
            vals = [y for (_1, _2, y) in shs]
            shift = sum(vals) * 1.0 / len(vals)
            residue[atom] = shift
            mx, mn, av = max(vals), min(vals), avg(vals)
            if (mx - av) > 0.2 or (av - mn) > 0.2:
                print resid, atom, mx, mn, av, vals
        best[resid] = residue
    return best


def map_keys(f, d):
    """not a real map operation since may be structure-changing => solution: throw exceptions!"""
    out = {}
    for (k, v) in d.items():
        new_key = f(k)
        if new_key in out:
            raise ValueError('oops, duplicate new key %s %s' % (new_key, d))
        out[new_key] = v
    return out


def manipulateShifts():
    shifts = getResidueShifts()
    one = model.fmap_dict(lambda s: model.fmap_dict(selectOneShift, s), shifts)
    two = model.fmap_dict(lambda s: map_keys(disambiguateAtom, s), one)
    three = model.fmap_dict(removeDuplicateAtoms, two)
    return three


def getShifts():
    shifts = manipulateShifts()
    def f1(x):
        return x.items()
    def f2(y):
        resid, pairs = y
        return [(resid, z[0], z[1]) for z in pairs]
    def concat(xs):
        return sum(xs, [])
    lists = map(f2, model.fmap_dict(f1, shifts).items())
    return concat(lists)



def printXEasyShifts():
    shifts = getShifts()
    num = 1
    for (resid, atom, shf) in shifts:
        print '{:4} {:7.3f} {:7} {:4} {:3}'.format(num, round(shf, 3), '0.000', atom, resid)
        num += 1


def setAtomType(specname, peakid, dimno, atomtypes=None):
    proj = getData(simple=False)
    if atomtypes is not None:
        old = proj.spectra[specname].peaks[peakid].dims[dimno].atomtypes
        proj.spectra[specname].peaks[peakid].dims[dimno].atomtypes = atomtypes 
        pt.json_out(ROOT + "project.txt", proj)
        print 'old value:', old
    else:
        print proj.spectra[specname].peaks[peakid].dims[dimno].atomtypes
#    return proj




# need to also assign all the CCONH peaks to be i-1 in the carbon dimension
# then assign atomtypes of C dimension of CCONH
# how about tag every CCONH peak not in a spin system to "unknown"?

def findSSNoCCONH():
    proj = getData()
    ss = [(s.id, s.tags, len(filter(lambda pkid: pkid[0] == 'cconh', s.pkids))) for s in proj.getSpinSystems()]
    return ss


def findCCONHWrongNumberOfSS():
    proj = getData()
    pks = dict((i, []) for i in range(1, 755 + 1))
    for ss in proj.getSpinSystems():
        for (name, pkid) in ss.pkids:
            if name == 'cconh':
                if not pks.has_key(pkid):
                    pks[pkid] = []
                pks[pkid].append(ss.id)
    return pks


def findUnassignedSS():
    yes, no = 0, 0
    for ss in getData().getSpinSystems():
        if len(ss.residueids) == 0:
            print 'ss: ', ss.id, ss.tags
            no += 1
        else:
            yes += 1
    print '\n'
    print 'yes: ', yes, '    no: ', no


def getSSPeaks(ssid, *specnames):
    """
    If no spectra names given, return all peaks of the specified SS.
    Otherwise, return only the peaks with names in specnames.
    """
    my_data = getData()
    ss = my_data._spinsystems[ssid]
    pkids = ss.pkids
    def include(name):
        if len(specnames) == 0:
            return True
        return name in specnames
    def getAts(s_name, pk_id):
#        return [dim for dim in my_data._spectra[s_name]._peaks[pk_id].dims]
        return my_data._spectra[s_name]._peaks[pk_id]
    return [(name, pkid, getAts(name, pkid)) for (name, pkid) in pkids if include(name)]


## ANALYSIS !!!

class Analysis(object):
    
    def __init__(self):
        self.data = getData(simple=False)
    
    def save(self, path=None):
        if path is None:
            savePath = ROOT + 'auto_project.txt'
        else:
            savePath = ROOT + path
        pt.json_out(savePath, self.data)
    
    def addPeakTag(self, specname, tag, *pkids):
        for pkid in pkids:
            peak = self.data.spectra[specname].peaks[pkid]
            if tag in peak.tags:
                raise ValueError('cannot add peak tag, already present')
            peak.tags.append(tag)
    
    def createNewPeak(self, specname, shifts, height):
        spec = self.data.spectra[specname]
        pkid = max(spec.peaks.keys()) + 1
        spec.addPeak(pkid, model.Peak([model.PeakDim(s, []) for s in shifts], [], height))
        return pkid
    
    def addPeakToSS(self, ssid, specname, *pkids):
        proj = self.data
        for pkid in pkids:
            ss = proj.spinsystems[ssid]
            if [specname, pkid] in ss.pkids:
                raise ValueError('unable to add peak to spin system: already present')
            if specname not in proj.spectra:
                raise ValueError('unrecognized spectrum')
            if pkid not in proj.spectra[specname].peaks:
                raise ValueError('invalid peak id')
            ss.pkids.append([specname, pkid])
    
    def addPeakDimAtomType(self, specname, dimNo, *pairs):
        peaks = self.data.spectra[specname].peaks
        for (pkid, atomtypes) in pairs:
            dim = peaks[pkid].dims[dimNo]
            for atom in atomtypes:
                if atom in dim.atomtypes:
                    raise ValueError(('cannot add atomtype: already present', atomtypes, peaks[pkid]))
                dim.atomtypes.append(atom)
    
    def removePeakFromSS(self, ssid, specname, *pkids):
        ss = self.data.spinsystems[ssid]
        for pkid in pkids:
            if not [specname, pkid] in ss.pkids:
                raise ValueError('unable to remove peak from spin system: not found')
            ss.pkids = filter(lambda x: x != [specname, pkid], ss.pkids)

    def assign_HBHACONH_HCCONH_SS(self, ssid, hbhaconh_bad, hbhaconh_new, hcconh_bad, hcconh_new):
        """
        still have to assign atomtypes to peak ids on my own
        """
        self.removePeakFromSS(ssid, 'hbhaconh', *hbhaconh_bad)
        self.addPeakTag('hbhaconh', 'artifact', *hbhaconh_bad)
        self.removePeakFromSS(ssid, 'hcconh', *hcconh_bad)
        self.addPeakTag('hcconh', 'artifact', *hcconh_bad)
        ids, zzz = [], [] # wow, what a horrible variable name
        for x in hbhaconh_new:
            newid = self.createNewPeak('hbhaconh', x[0], x[1])
            ids.append(newid)
        for y in hcconh_new:
            newid = self.createNewPeak('hcconh', y[0], y[1])
            zzz.append(newid)
        self.addPeakToSS(ssid, 'hbhaconh', *ids)
        self.addPeakToSS(ssid, 'hcconh', *zzz)
        return (ids, zzz)


if __name__ == "__main__":
    # findCloseNHSQCPeaks()
    # peakListOut()
    # findJunkPeaksInSpinSystems()
    # addTag()
    # checkTags()
    # analyzeSpinSystems()
    # findNHSQCStrangenessAndAddHNCACBPeaksToSpinSystems()
    pass
