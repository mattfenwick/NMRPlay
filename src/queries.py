from . import inout as pt
from . import model
from .algebra import ffilter, inner_join
from operator import attrgetter


ROOT = '../../PeakPicker/'

def getData(root=ROOT):
    return pt.json_in(root + 'project.txt')


def filterSpecPeaks(pred, spec):
    '''
    Filter peaks from a spectrum based on a predicate,
    returning a new spectrum.
    '''
    peaks = ffilter(pred, spec.peaks)
    return model.Spectrum(spec.axes, peaks)


def getAllPeaks(prj):
    """
    Project -> [(String, Int, Peak)]
    """
    peaks = []
    for (name, spectrum) in prj.spectra.iteritems():
        peaks.extend([(name, pkid, pk) for (pkid, pk) in spectrum.peaks.iteritems()])
    return peaks


def joinSSToPeaks(prj):
    joined, ss = {}, prj.spinsystems
    for (ssid, spinsys) in ss.iteritems():
        pks = [prj.spectra[pk_spectrum].peaks[pk_id] for (pk_spectrum, pk_id) in spinsys.pkids]
        joined[ssid] = pks
    return joined

def joined2(prj):
    """
    Project -> [(Int, SpinSystem, [(Int, Peak)])]
    """
    peaks = getAllPeaks(prj)
    def p(ss, pk):
        return list(pk[:2]) in ss[1].pkids
    joined = inner_join(p, prj.spinsystems.iteritems(), peaks)
    return joined


def findSpinsystemsOfPeak(spins, spec, pid):
    """
    Map Int SpinSystem -> String -> Int -> [Int]
    """
    ss_ids = []
    for (ssid, ss) in spins.items():
        if [spec, pid] in ss.pkids:
            ss_ids.append(ssid)
    return ss_ids


def findCloseNHSQCPeaks(proj, HTOL=0.025, NTOL=0.2):
    """
    Project -> [((Int, Peak, [SpinSystem]), (Int, Peak, [SpinSystem]))]
    """
    nhsqc = proj.spectra['nhsqc']
    spins = proj.spinsystems
    
    if nhsqc.axes != ['N', 'H']:
        raise ValueError('oops, unexpected axis order!')
    
    for (n1id, n1pk) in nhsqc.peaks.items():
        for (n2id, n2pk) in nhsqc.peaks.items():
            n1, h1 = map(lambda x: x.shift, n1pk.dims)
            n2, h2 = map(lambda x: x.shift, n2pk.dims)
            if abs(n2 - n1) <= NTOL and abs(h2 - h1) <= HTOL and n1id < n2id:
                print 'peak1', n1id, n1pk
                print 'peak2', n2id, n2pk
                print 'spin systems1: ', map(lambda x: spins[x], findSpinsystemsOfPeak(spins, 'nhsqc', n1id))
                print 'spin systems2: ', map(lambda x: spins[x], findSpinsystemsOfPeak(spins, 'nhsqc', n2id))
                print '\n'


def analyzeSpinSystems(proj):
    """
    Project -> ([(SpinSystemID, [Float], [(Int, String)])], Map Int Int)
    """
    nhsqc, hnco = proj.spectra['nhsqc'], proj.spectra['hnco']
    spins = proj.spinsystems
    
    # 1. find lone NHSQC peaks
    #   - get list of all NHSQC peak ids -- wait, isn't this unnecessary since I have a spin system for each NHSQC peak?
    #   - figure out which spin systems don't have any hnco peaks
    print 'spin systems with more or less than 1 HNCO peak:'
    for (ssid, ss) in spins.items():
        if len(ss.pkids) != 2 and len(nhsqc.peaks[ss.pkids[0][1]].tags) == 0:
            print ssid, '    ', map(lambda x: x.shift, nhsqc.peaks[ssid].dims), '    ', ss.pkids
    print '\n'
    
    # 2. find ??? HNCO peaks -- member of < or > 1 spin system
    #   - get list of HNCO peak ids
    #   - look through all spin systems, and total up the hnco peak ids in each
    hncoids = dict((pkid, 0) for (pkid, pk) in hnco.peaks.items() if 'backbone amide' in pk.tags)
    for (ssid, ss) in spins.items():
        for (spec, pid) in ss.pkids:
            if spec == 'hnco' and hncoids.has_key(pid):
                hncoids[pid] += 1
    for (key, val) in hncoids.items():
        if val != 1:
            print 'hnco peak', key, 'is interesting: ', val, 'spin systems!'
        else:
            print key, 'is boring'


def checkTags(proj):
    """
    Project -> [(String, PeakID)]
    """
    for (_, pk) in proj.spectra['hnco'].peaks.items():
        if len(pk.tags) != 1:
            print 'problem'
        else:
            print 'good', _


def peakListOut(proj):
    """
    Project -> IO ()
    """
    nhsqc = proj.spectra['nhsqc']
    
    proj.spectra['nhsqc_good'] = filterSpecPeaks(lambda pk: 'backbone amide' in pk.tags, nhsqc)
    proj.spectra['nhsqc_bad'] = filterSpecPeaks(lambda pk: 'backbone amide' not in pk.tags, nhsqc)
    
    pt.xeasy_out(proj, {'nhsqc_good': ROOT + "nhsqc_good.txt", 'nhsqc_bad': ROOT + "nhsqc_bad.txt"})


def findJunkPeaksInSpinSystems(proj):
    """
    Project -> ([(Int, Int, [String])], [(String, String, [String], [Float])])
    """
    hnco = proj.spectra['hnco']
    pks = filterSpecPeaks(lambda pk: pk.tags != ['backbone amide'], hnco).peaks
    for (ssid, ss) in proj.spinsystems.items():
        for (spec, pkid) in ss.pkids:
            if spec == 'hnco' and pkid in pks:
                print 'uh-oh, problem with ', pkid, 'in spin system', ssid, 'with tags', pks[pkid].tags
    hncoids = dict((pkid, []) for (pkid, pk) in hnco.peaks.items())
    for (ssid, ss) in proj.spinsystems.items():
        for (spec, pid) in ss.pkids:
            if spec == 'hnco':
                hncoids[pid].append(ssid)
    for (key, val) in hncoids.items():
        print key, '    ', val, '    ', hnco.peaks[key].tags, '    ', map(lambda x: x.shift, hnco.peaks[key].dims)
    
    
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
    nhsqc, hnco = proj.spectra['nhsqc'], proj.spectra['hnco']
    
    shf = attrgetter('shift')
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
    


if __name__ == "__main__":
    # findCloseNHSQCPeaks()
    # peakListOut()
    # findJunkPeaksInSpinSystems()
    # addTag()
    # checkTags()
    # analyzeSpinSystems()
    # findNHSQCStrangenessAndAddHNCACBPeaksToSpinSystems()
    findHNCACBStrangeness()

