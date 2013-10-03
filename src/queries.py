from . import inout as pt
from . import model
from . import querymodel
from .bmrbshifts import stats
from .algebra import ffilter, inner_join, groupBy, split, concatMap, fmap
from operator import attrgetter


ROOT = '../PeakPicker/'

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


def peakListOut(proj):
    """
    Project -> IO ()
    """
    nhsqc = proj._spectra['nhsqc']
    peaks = nhsqc.getPeaks()
    
    proj.spectra['nhsqc_good'] = filter(lambda pk: 'backbone amide' in pk.tags, peaks)
    proj.spectra['nhsqc_bad'] = filter(lambda pk: 'backbone amide' not in pk.tags, peaks)
    
    pt.xeasy_out(proj, {'nhsqc_good': ROOT + "nhsqc_good.txt", 'nhsqc_bad': ROOT + "nhsqc_bad.txt"})


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
    
    
def addTag(proj):
    """
    Project -> IO Project
    """
    spec = proj.spectra['nhsqc']
    for (_, pk) in spec.peaks.items():
        if len(pk.tags) == 0:
            print '          ', _
            pk.tags.append('backbone amide')
        else:
            print _
    pt.json_out(ROOT + "new_project.txt", proj)


def buildSpinSystems(proj, HTOL=0.025, NTOL=0.2):
    """
    Project -> IO Project
    builds spin systems along the way based on matching
    """
    nhsqc, hnco = proj._spectra['nhsqc'], proj._spectra['hnco']
    
    shf = lambda x: x.shift
    ssid = 1
    
    print hnco.axes, nhsqc.axes
                
    for (nid, npk) in nhsqc.peaks.items():
        ss = model.SpinSystem([['nhsqc', nid]], [], [], [])
        print 'ssid:', ssid
        for (hid, hpk) in hnco.peaks.items():
            nhsqc_n, nhsqc_h = map(shf, npk.dims)
            hnco_h, hnco_n, hnco_c = map(shf, hpk.dims)
            if abs(nhsqc_n - hnco_n) <= NTOL and abs(nhsqc_h - hnco_h) <= HTOL:
                ss.pkids.append(['hnco', hid])
                print 'yes'
        proj.spinsystems[ssid] = ss
        ssid += 1
        
    pt.json_out(ROOT + "project_new.txt", proj)


def dumpHncacb(proj):
    """
    Project -> IO Spectrum
    """
    spec = proj.spectra['hncacb']
    pks = dict(filter(lambda (pid, pk): len(pk.tags) == 0, spec.peaks.items()))
    pt.xeasy_out(model.Project('name', 
                             {'hncacb': model.Spectrum(spec.axes, pks)}, 
                             proj.molecule, proj.spinsystems), {'hncacb': ROOT + 'hncacb_oops.txt'})


def addTags(proj):
    """
    Project -> IO Project
    """
    hncacb = pt.xeasy_in("temp", {'hncacb': ROOT + 'hncacb.xez'}).spectra['hncacb']

    print proj.spectra['hncacb'] == hncacb

    i = 1
    for (pid, pk) in proj.spectra['hncacb'].peaks.items():
        if pid not in hncacb.peaks.keys():
            if len(pk.tags) == 0:
                pk.tags.append('processing artifact')
                print pid, pk, i
                i += 1

    pt.json_out(ROOT + "project2.txt", proj)


def anImport():
    """
    () -> IO Project
    """
    # ba = pt.xeasy_in('mynameismatt', {'nhsqc': 'nhsqc.xez', 'hnco': 'hnco.xez'})
    ba = pt.xeasy_in('mynameismatt', {'nhsqc' : ROOT + 'nhsqc.xez', 
                                      'hnco'  : ROOT + 'hnco.xez' , 
                                      'hncacb': ROOT + 'hncacb.xez'})
    pt.json_out(ROOT + 'new_project.txt', ba)
    
    
def findNHSQCStrangenessAndAddHNCACBPeaksToSpinSystems(proj):
    """
    Project -> ???
    """
    nhsqc, hncacb = proj.spectra['nhsqc'], proj.spectra['hncacb']
    spins = proj.spinsystems
    
    nh_peaks = model.fmap_dict(lambda _: [], nhsqc.peaks)
    for (ssid, ss) in spins.items():
        pks = []
        for (spec, pkid) in ss.pkids:
            if spec == 'nhsqc':
                nh_peaks[pkid].append(ssid)
                pks.append(pkid)
        if len(pks) != 1:
            print 'spin system', ssid, 'has weird peaks:', pks
    
    def action(xs):
        if len(xs) != 1:
            raise ValueError('assuming each NHSQC peak assigned to exactly 1 spin system')
        return xs[0]
    
    # nhsqc_spins :: Map NHSQC_PK_ID SPIN_SYSTEM_ID
    nhsqc_spins = model.fmap_dict(action, nh_peaks)
    
    print nhsqc_spins
    
    HTOL, NTOL = 0.0125, 0.1 # make tolerances smaller b/c spectra are well-aligned
    shf = lambda x: x.shift
    
    assert nhsqc.axes == ["N", "H"]
    assert hncacb.axes == ["H", "N", "C"]
    assnCnt = 0
    
    for (nid, npk) in nhsqc.peaks.items():
        ss = spins[nhsqc_spins[nid]]
        for (hid, hpk) in hncacb.peaks.items():
            nhsqc_n, nhsqc_h = map(shf, npk.dims)
            hncacb_h, hncacb_n, hncacb_c = map(shf, hpk.dims)
            if abs(nhsqc_n - hncacb_n) <= NTOL and abs(nhsqc_h - hncacb_h) <= HTOL:
                ss.pkids.append(['hncacb', hid])
                assnCnt += 1
                print 'assignment', assnCnt, "nhsqc id:", nid, "hncacb id:", hid
                print npk
                print hpk
                print
                
    pt.json_out(ROOT + "project_new2.txt", proj)


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


def classifyHNCACBPeaks():
    """ height > 0 -> C-Alpha, height < 0 -> C-Beta """
    proj = getData(simple=False)
    for (i, (_, pk)) in enumerate(filter(lambda x: x[1].tags == [], proj.spectra['hncacb'].peaks.items())):
        if pk.height > 0:
            t = 'CA'
        else:
            t = 'CB'
        pk.dims[2].atomtypes.append(t)
        pk.tags.append('backbone')
        print i, pk
    pt.json_out(ROOT + "project4.txt", proj)


def dumpDataToXEasy():
    """
    Project -> IO ()
    """
    proj = getData(simple=False)
    spectra = [('nhsqc', 'backbone amide'),
               ('hnco', 'backbone amide'),
               ('hncacb', 'backbone')]
    outs = {}
    for (name, tag) in spectra:
        spec = proj.spectra[name]
        good_peaks = filterSpecPeaks(lambda pk: tag in pk.tags, spec)
        bad_peaks = filterSpecPeaks(lambda pk: tag not in pk.tags, spec)
        proj.spectra[name + "_good"] = good_peaks
        proj.spectra[name + "_bad" ] = bad_peaks
        proj.spectra[name + "_full"] = spec
        outs[name + "_good"] = ROOT + name + "_good.peaks"
        outs[name + "_bad" ] = ROOT + name + "_bad.peaks"
        outs[name + "_full"] = ROOT + name + "_full.peaks"
    
    pt.xeasy_out(proj, outs)


def findHNCACBPeaksInSSAndTaggedBackbone():
    """
    HNCACB peaks in a spin system but not tagged 'backbone' or vice versa
    """
    proj = getData()
    # map (get 1) . filter ((==) "hncacb" . get 0) $ concatMap (get pkids) (getSpinSystems proj)
    inss = map(lambda x: x[1], filter(lambda n: n[0] == 'hncacb', concatMap(lambda ss: ss.pkids, proj.getSpinSystems())))
    back = map(lambda x: x.id, filter(lambda pk: pk.tags == ['backbone'], proj._spectra['hncacb'].getPeaks()))
    return inss, back


def addSequence():
    proj = getData(simple=False)
    proj.molecule.residues = list('GGGRDYKDDDDKGTMELELRFFATFREVVGQKSIYWRVDDDATVGDVLRSLEAEYDGLAGRLIEDGEVKPHVNVLKNGREVVHLDGMATALDDGDAVSVFPPVAGG')
    pt.json_out(ROOT + "project5.txt", proj)


def handleArginineSidechains():
    proj = getData(simple=False)
    for ssid in [99, 111, 141]:
        ss = proj.spinsystems[ssid]
        ss.tags.append('folded: NE once, CZ twice')
        ss.tags.append('sidechain')
        ss.aatypes.append('arginine')
        for (name, pkid) in ss.pkids:
            peak = proj.spectra[name].peaks[pkid]
            if name == 'nhsqc':
                peak.dims[0].atomtypes.append('NE')
                peak.dims[1].atomtypes.append('HE')
            elif name == 'hnco':
                peak.dims[0].atomtypes.append('HE')
                peak.dims[1].atomtypes.append('NE')
                peak.dims[2].atomtypes.append('CZ')
            elif name == 'hncacb':
                peak.dims[0].atomtypes.append('HE')
                peak.dims[1].atomtypes.append('NE')
                at = peak.dims[2].atomtypes
                if at == ['CA']:
                    peak.dims[2].atomtypes = ['CD']
                elif at == ['CB']:
                    peak.dims[2].atomtypes = ['CG']
                else:
                    peak.dims[2].atomtypes = ['CD']
                    print 'atomtype problem with: ', ssid, name, pkid, peak
#                    raise ValueError('inner oops ' + str(at))
                peak.tags = ['arginine sidechain']
            else:
                raise ValueError('outer oops')
#        print ss
    pt.json_out(ROOT + "project7.txt", proj)


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


def removeHNCACBJunkFromSS():
    proj = getData(simple=False)
    pks = proj.spectra['hncacb'].peaks
    def predicate(p):
        name, pkid = p
        if name == 'hncacb':
            if 'processing artifact' in pks[pkid].tags:
                if pks[pkid].tags != ['processing artifact']:
                    raise ValueError('oops')
                else:
                    return False
            else: # hncacb but not an artifact, keep it
                return True 
        else: # not hncacb, keep it
            return True
    dumped = set([])
    for (_, ss) in proj.spinsystems.items():
        before = [i for (_, i) in ss.pkids if _ == 'hncacb']
        ss.pkids = filter(predicate, ss.pkids)
        after = [i for (_, i) in ss.pkids if _ == "hncacb"]
        dumped = dumped.union(set(before) - set(after))
    pt.json_out(ROOT + "project9.txt", proj)
    return dumped


def inOut():
    proj = getData(simple=False)
    pt.json_out(ROOT + "project10.txt", proj)


def retagSpinSystems():
    proj = getData(simple=False)
    for (_, s) in proj.spinsystems.items():
        if s.tags == []:
            s.tags.append('backbone')
    pt.json_out(ROOT + "project11.txt", proj)


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


def getShifts():
    """
    this is kind of wrong ... it's only the CA/CB of HNCACB
    peaks that is i/i-1, the H/N shifts are always i
    """
    proj = getData()
    jed = joined2(proj.getSpinSystems(), getAllPeaks(proj))
    ged = groupBy(lambda x: x[0].id, jed, snd)
    # proceed down the spin systems, adding the i peaks,
    #   then moving to the i-1 peaks
    atomshifts = fmap(lambda pks: concatMap(lambda pk: map(lambda d: (d, [t for t in pk.tags if t in set(['i', 'i-1'])]), pk.dims), pks), ged)
    tags = dict(map(lambda x: (x.id, x.tags), proj.getSpinSystems()))
    for x in tags:
        if not x in atomshifts:
            raise ValueError('oops')
        tags[x] = (tags[x], atomshifts[x])
    things = fmap(lambda y: (y[0], sorted(map(lambda z: (z[0][0], z[0][1], z[1]), y[1]), key=lambda x:x[0])), tags)
    return things


def getResidueShifts():
    def getI(specname, atomtypes, dimname, pktags):
        if (specname, dimname) == ('hncacb', 'C'):
            if 'i-1' in pktags and 'i' in pktags:
                return 'both'
            elif 'i-1' in pktags:
                return 'i-1'
            elif 'i' in pktags:
                return 'i'
            else:
                return None
        if specname == 'hcchtocsy':
            return 'i-1'
        if (specname, dimname) == ('cconh', 'C'):
            return 'i-1'
        if (specname, dimname) == ('hcconh', 'h'):
            return 'i-1'
        if (specname, dimname) == ('hbhaconh', 'h'):
            return 'i-1'
        # what about: cb_he, cb_hd ?? may have them assigned to wrong spin systems
        return 'i' # nhsqc, hnco
            
    ats = {
        ('nhsqc', 'N'): ('N'),
        ('nhsqc', 'H'): ('H'),
        ('hncacb', 'N'): ('N'),
        ('hncacb', 'H'): ('H'),
        ('hnco', 'N'): ('N'),
        ('hnco', 'H'): ('H'),
        ('hnco', 'C'): ('CO'),
        ('cconh', 'N'): ('N'),
        ('cconh', 'H'): ('H'),
        ('hcconh', 'N'): ('N'),
        ('hcconh', 'H'): ('H'),
        ('hbhaconh', 'N'): ('N'),
        ('hbhaconh', 'H'): ('H'),
    }
    def getAt(specname, atomtypes, dimname):
        if len(atomtypes) > 0:
            return tuple(atomtypes)
        if (specname, dimname) in ats:
            return ats[(specname, dimname)]
    # 1. get residue-SS assignment
    proj = getData()
    import collections
    res = dict((i, collections.defaultdict(lambda: [])) for (i, _) in enumerate(proj.molecule.residues, start=1))
    asses = {}
    for s in proj.getSpinSystems():
        if len(s.residueids) == 0:
            continue
        if len(s.residueids) != 1:
            raise ValueError('unexpected number of residues assigned to spin system -- %s' % str(s.residueids))
        asses[s.id] = s.residueids[0]
    print asses, '\n\n', res
#    return 0
    # 2. get SS peaks
    for ss in proj.getSpinSystems():
        if ss.id not in asses:
            continue
        for (specname, pkid) in ss.pkids:
            peak = proj._spectra[specname]._peaks[pkid]
            for d, dimname in zip(peak.dims, proj._spectra[specname].axes):
                try:
                    pos = getI(specname, d[1], dimname, peak.tags)
                except:
                    print peak
                    raise
                atoms = getAt(specname, d[1], dimname)
                record = (specname, d[0])
                if pos == 'i':
                    res[asses[ss.id]    ][atoms].append(record)
                elif pos == 'i-1':
                    res[asses[ss.id] - 1][atoms].append(record)
                elif pos == 'both':
                    res[asses[ss.id]    ][atoms].append(record)
                    res[asses[ss.id] - 1][atoms].append(record)
                elif pos is None:
                    pass
                else:
                    raise ValueError('invalid')
    return model.fmap_dict(dict, res)
    # 3. get residue peaks
    # 4. use i/i-1 experiment knowledge, plus backbone/sidechain tags and aatype, to convert peaks to correct residue
    # 5. since we'll have multiple measurements for each atom, figure out some way to collapse the list to a single value
    # 6. for each peak: atomname, spectrum name, shift


def loadCCONH():
    """
    Project -> IO Spectrum
    """
    proj = getData(simple=False)
    xs = pt.xeasy_in('oops', {'cconh': ROOT + "cconh.xez"})
    proj.spectra['cconh'] = xs.spectra['cconh']
    pt.json_out(ROOT + "another_new.txt", proj)


def groupCCONHPeaksIntoSpinSystems(HTOL=0.025, NTOL=0.2):
    proj = getData(simple=False)
    nhsqc, cconh = proj.spectra['nhsqc'], proj.spectra['cconh']
    
    shf = lambda x: x.shift
    
    print cconh.axes, nhsqc.axes
    if cconh.axes != ['H', 'N', 'C']:
        raise ValueError(('oops, axis order', cconh.axes))
    
    for ssid, ss in proj.spinsystems.items():
        for name, nid in ss.pkids:
            if name != 'nhsqc':
                continue
            for (cid, cpk) in cconh.peaks.items():
                nhsqc_n, nhsqc_h = map(shf, nhsqc.peaks[nid].dims)
                cconh_h, cconh_n, cconh_c = map(shf, cpk.dims)
                if abs(nhsqc_n - cconh_n) <= NTOL and abs(nhsqc_h - cconh_h) <= HTOL:
                    ss.pkids.append(['cconh', cid])
                    print ssid, cid, nhsqc_n, nhsqc_h, cconh_n, cconh_h
        
    pt.json_out(ROOT + "new_project2.txt", proj)


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


def retagUnassignedCCONHPeaks():
    pkAssigns = findCCONHWrongNumberOfSS()
    proj = getData(simple=False)
    for (pkid, ss) in pkAssigns.items():
        if len(ss) == 0:
            proj.spectra['cconh'].peaks[pkid].tags.append('unknown')
        # otherwise, good to go
    pt.json_out(ROOT + "new_project3.txt", proj)


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


def loadHCCONHTocsy():
    """
    Project -> IO Spectrum
    """
    proj = getData(simple=False)
    xs = pt.xeasy_in('oops', {'hcconh': ROOT + "hcconh_tocsy.xez"})
    proj.spectra['hcconh'] = xs.spectra['hcconh']
    pt.json_out(ROOT + "again_new.txt", proj)


def groupHCCONHTocsyPeaksIntoSpinSystems(HTOL=0.025, NTOL=0.2):
    proj = getData(simple=False)
    nhsqc, hcconh = proj.spectra['nhsqc'], proj.spectra['hcconh']
    
    shf = lambda x: x.shift
    
    print hcconh.axes, nhsqc.axes
    if hcconh.axes != ['H', 'N', 'h']:
        raise ValueError(('oops, axis order', hcconh.axes))
    
    for ssid, ss in proj.spinsystems.items():
        for name, nid in ss.pkids:
            if name != 'nhsqc':
                continue
            for (cid, cpk) in hcconh.peaks.items():
                nhsqc_n, nhsqc_h = map(shf, nhsqc.peaks[nid].dims)
                hcconh_h, hcconh_n, hcconh_other_h = map(shf, cpk.dims)
                if abs(nhsqc_n - hcconh_n) <= NTOL and abs(nhsqc_h - hcconh_h) <= HTOL:
                    ss.pkids.append(['hcconh', cid])
                    print ssid, cid, nhsqc_n, nhsqc_h, hcconh_n, hcconh_h
        
    pt.json_out(ROOT + "new_project2.txt", proj)


def loadHbhaconh():
    """
    Project -> IO Spectrum
    """
    proj = getData(simple=False)
    xs = pt.xeasy_in('oops', {'hbhaconh': ROOT + "hbhaconh.xez"})
    proj.spectra['hbhaconh'] = xs.spectra['hbhaconh']
    pt.json_out(ROOT + "again_mew.txt", proj)


def groupHbhaconhPeaksIntoSpinSystems(HTOL=0.025, NTOL=0.2):
    proj = getData(simple=False)
    nhsqc, hbhaconh = proj.spectra['nhsqc'], proj.spectra['hbhaconh']
    
    shf = lambda x: x.shift
    
    print hbhaconh.axes, nhsqc.axes
    if hbhaconh.axes != ['H', 'N', 'h']:
        raise ValueError(('oops, axis order', hbhaconh.axes))
    
    for ssid, ss in proj.spinsystems.items():
        for name, nid in ss.pkids:
            if name != 'nhsqc':
                continue
            for (cid, cpk) in hbhaconh.peaks.items():
                nhsqc_n, nhsqc_h = map(shf, nhsqc.peaks[nid].dims)
                hbhaconh_h, hbhaconh_n, hbhaconh_other_h = map(shf, cpk.dims)
                if abs(nhsqc_n - hbhaconh_n) <= NTOL and abs(nhsqc_h - hbhaconh_h) <= HTOL:
                    ss.pkids.append(['hbhaconh', cid])
                    print ssid, cid, nhsqc_n, nhsqc_h, hbhaconh_n, hbhaconh_h
        
    pt.json_out(ROOT + "new_project2.txt", proj)


def getSSPeaks(ssid, *specnames):
    """
    If no spectra names given, return all peaks of the specified SS.
    Otherwise, return only the peaks with names in specnames.
    """
    ss = getData()._spinsystems[ssid]
    pkids = ss.pkids
    if len(specnames) == 0:
        return pkids
    return [(name, pkid) for (name, pkid) in pkids if name in specnames]


def loadAromaticLinks():
    """
    Project -> IO Spectrum
    """
    proj = getData(simple=False)
    xs = pt.xeasy_in('oops', {'cb_he': ROOT + "hbCBcgcdceHE.txt", 
                              'cb_hd': ROOT + "hbCBcgcdHD.txt"})
    proj.spectra['cb_he'] = xs.spectra['cb_he']
    proj.spectra['cb_hd'] = xs.spectra['cb_hd']
    pt.json_out(ROOT + "again_its_new.txt", proj)


def retagAromaticPeaks():
    proj = getData(simple=False)
    for pk in proj.spectra['cb_hd'].peaks.values():
        pk.dims[0].atomtypes.append('CB')
        pk.dims[1].atomtypes.append('HD*')
    for pk in proj.spectra['cb_he'].peaks.values():
        pk.dims[0].atomtypes.append('CB')
        pk.dims[1].atomtypes.append('HE*')
    pt.json_out(ROOT + "again_were_new.txt", proj)


def loadHCCHTocsy():
    trans = {
        ('HA',)             : ['CA'], 
        ('HA', 'HB')        : ['CA'],
        ('HA2', 'HA3')      : ['CA'],
        ('HB',)             : ['CB'],
        ('HB2', 'HB3')      : ['CB'],
        ('HD1',)            : ['CD1'],
        ('HD1', 'HD2')      : ['CD1', 'CD2'],
        ('HD2', 'HD3')      : ['CD'],
        ('HE2', 'HE3')      : ['CE'],
        ('HG',)             : ['CG'],
        ('HG', 'HD1', 'HD2'): ['CG', 'CD1', 'CD2'],
        ('HG1', 'HG2')      : ['CG1', 'CG2'],
        ('HG12', 'HG13')    : ['CG1', 'CG2'],
        ('HG2',)            : ['CG2'],
        ('HG2', 'HG3')      : ['CG']
    }
    proj = getData(simple=False)
    xs = pt.xeasy_in('?', {'hcchtocsy': ROOT + 'hcchtocsy.xez'})
    proj.spectra['hcchtocsy'] = xs.spectra['hcchtocsy']
    with open(ROOT + 'hcch_tocsy_atomtypes.txt') as infile:
        atomtypes = eval(infile.read())
    um = {}
    import collections
    seen = collections.defaultdict(lambda: 0)
    for (ssid, assns) in atomtypes.iteritems():
        for (peak_id, atoms) in assns:
            seen[peak_id] += 1
            # 1. add ('hcchtocsy', $peak_id$) to the appropriate spin systems' peaks
            proj.spinsystems[ssid].pkids.append(['hcchtocsy', peak_id])
            # 2. assign the 1H dimensions of the peak
            peak = proj.spectra['hcchtocsy'].peaks[peak_id]
            peak.dims[0].atomtypes = atoms + ['i-1']
            peak.dims[2].atomtypes = atoms + ['i-1']
            # 3. translate the 1H assignments to 13C assignments, and assign those to the C-dimension
            peak.dims[1].atomtypes = trans[tuple(atoms)] + ['i-1']
    pt.json_out(ROOT + "incomplete.txt", proj)
    return seen
    # the peaks look fine (they're all on the diagonal)
#    for pk in xs.spectra['hcchtocsy'].peaks.values():
#        if abs(pk.dims[0].shift - pk.dims[2].shift) > 0.02:
#            print 'oops: ', pk


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
