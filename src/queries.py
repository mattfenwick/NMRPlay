from . import inout as pt
from . import model
from .algebra import ffilter, inner_join, groupBy, split, concatMap
from .querymodel import fromModel # bring this into scope b/c it seems convenient ... or something
from operator import attrgetter


ROOT = '../PeakPicker/'

def getData(root=ROOT, simple=True):
    data = pt.json_in(root + 'project.txt')
    if simple:
        return fromModel(data)
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


def dumpGoodDataToXEasy():
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
        proj.spectra[name + "_good"] = good_peaks
        outs[name + "_good"] = ROOT + name + "_good.txt"
    
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


if __name__ == "__main__":
    # findCloseNHSQCPeaks()
    # peakListOut()
    # findJunkPeaksInSpinSystems()
    # addTag()
    # checkTags()
    # analyzeSpinSystems()
    # findNHSQCStrangenessAndAddHNCACBPeaksToSpinSystems()
    pass
