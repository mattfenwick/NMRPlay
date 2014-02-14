from . import inout as pt
from . import model
from . import querymodel
from .bmrbshifts import stats
from .algebra import ffilter, inner_join, groupBy, split, concatMap, fmap
from operator import attrgetter
from . import queries as q # or, maybe -- from .queries import * ??


ROOT = '../PPAbnormal/'


def peakListOut(proj):
    """
    Project -> IO ()
    """
    nhsqc = proj._spectra['nhsqc']
    peaks = nhsqc.getPeaks()
    
    proj.spectra['nhsqc_good'] = filter(lambda pk: 'backbone amide' in pk.tags, peaks)
    proj.spectra['nhsqc_bad'] = filter(lambda pk: 'backbone amide' not in pk.tags, peaks)
    
    pt.xeasy_out(proj, {'nhsqc_good': ROOT + "nhsqc_good.txt", 'nhsqc_bad': ROOT + "nhsqc_bad.txt"})


    
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


def moveTagsIntoAtomtypes():
    proj = getData(simple=False)
    counts = [0, 0]
    for ssid in sorted(proj.spinsystems):
        ss = proj.spinsystems[ssid]
        if 'backbone' in ss.tags:
            if len(ss.residueids) == 0:
                print 'ss ', ssid, ' unassigned'
                continue
            print 'doing ss ', ssid
            counts[0] += 1
            pk_cts = [0, set([])]
            for (specname, pkid) in ss.pkids:
                pk = proj.spectra[specname].peaks[pkid]
                if specname == 'hnco':
                    pk.dims[0].atomtypes.append('H')
                    pk.dims[1].atomtypes.append('N')
                    pk.dims[2].atomtypes.append('C(i-1)')
                elif specname == 'nhsqc':
                    pk.dims[0].atomtypes.append('N')
                    pk.dims[1].atomtypes.append('H')
                elif specname == 'hncacb':
                    pk.dims[0].atomtypes.append('H')
                    pk.dims[1].atomtypes.append('N')
                    ats = pk.dims[2].atomtypes
                    if ats not in [['CA'], ['CB']]:
                        raise ValueError('bad atomtypes -- %s %s %s' % (ats, ssid, pkid))
                    atom, new_ats = ats[0], []
                    if 'i-1' in pk.tags:
                        new_ats.append(atom + '(i-1)')
                    if 'i' in pk.tags: # not `elif` b/c it can be both i and i-1
                        new_ats.append(atom)
                    pk.dims[2].atomtypes = new_ats
                elif specname == 'hcconh':
                    pk.dims[0].atomtypes.append('H')
                    pk.dims[1].atomtypes.append('N')
                    pk.dims[2].atomtypes = [('%s(i-1)' % a) for a in pk.dims[2].atomtypes]
                elif specname == 'hbhaconh':
                    pk.dims[0].atomtypes.append('H')
                    pk.dims[1].atomtypes.append('N')
                    pk.dims[2].atomtypes = [('%s(i-1)' % a) for a in pk.dims[2].atomtypes]
                elif specname == 'cconh':
                    pk.dims[0].atomtypes.append('H')
                    pk.dims[1].atomtypes.append('N')
                    pk.dims[2].atomtypes = [('%s(i-1)' % a) for a in pk.dims[2].atomtypes]
                elif specname == 'hcchtocsy':
                    pk.dims[0].atomtypes = [('%s(i-1)' % a) for a in pk.dims[0].atomtypes]
                    pk.dims[1].atomtypes = [('%s(i-1)' % b) for b in pk.dims[1].atomtypes]
                    pk.dims[2].atomtypes = [('%s(i-1)' % c) for c in pk.dims[2].atomtypes]
                elif specname == 'cb_hd': # 'hbCBcgcdHD':
                     pass
#                    pk.dims[0].atomtypes = [('%s(i-1)' % a) for a in pk.dims[0].atomtypes]
#                    pk.dims[1].atomtypes = [('%s(i-1)' % b) for b in pk.dims[1].atomtypes]
                elif specname == 'cb_he': # 'hbCBcgcdceHE':
                     pass
#                    pk.dims[0].atomtypes = [('%s(i-1)' % a) for a in pk.dims[0].atomtypes]
#                    pk.dims[1].atomtypes = [('%s(i-1)' % b) for b in pk.dims[1].atomtypes]
                else:
                    raise ValueError('invalid spectrum name -- %s' % specname)
                print pk
            # print 'peak counts: ', pk_cts
        else:
            print 'not bothering with: ', ssid, '(', ss.tags, ')'
            counts[1] += 1
    print 'yeses, nos: ', counts
    pt.json_out(ROOT + "project2.txt", proj)


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


def retagUnassignedCCONHPeaks():
    pkAssigns = findCCONHWrongNumberOfSS()
    proj = getData(simple=False)
    for (pkid, ss) in pkAssigns.items():
        if len(ss) == 0:
            proj.spectra['cconh'].peaks[pkid].tags.append('unknown')
        # otherwise, good to go
    pt.json_out(ROOT + "new_project3.txt", proj)


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


def removeRedundantHCCHTocsyIMinusOneAtomtypes():
    proj = getData(simple=False)
    for pk in proj.spectra['hcchtocsy'].peaks.values():
        for dim in pk.dims:
            if 'i-1' in dim.atomtypes:
                dim.atomtypes.remove('i-1')
    pt.json_out(ROOT + "again_without_i-1s.txt", proj)


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
