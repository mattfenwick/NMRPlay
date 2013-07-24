
class MyBase(object):
    '''
    Provides:
     - standard __repr__ serialization
     - standard value-based equality
    '''
    
    def __repr__(self):
        return repr(self.__dict__)
    
    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class Peak(MyBase):
    
    def __init__(self, spectrum_name, pkid, dims, tags, height):
        for d in dims:
            if len(d) != 2:
                raise TypeError(('peak dimension', d))
        if not isinstance(pkid, int):
            raise TypeError(('peak id', pkid))
        self.spectrum_name = spectrum_name
        self.id = pkid
        self.dims = dims
        self.tags = tags
        self.height = height
    
    def getPosition(self):
        return [d[0] for d in self.dims]
    
    def getAtomtypes(self):
        return [d[1] for d in self.dims]


class Spectrum(MyBase):
    '''
    This isn't so much a spectrum as it is a bunch of peaks from the same spectrum
    type.  Tags are responsible for keeping separate 'peaklists'.
    '''
    
    def __init__(self, name, axes, peaks):
        self.name = name
        self.axes = axes
        self._peaks = {}
        for p in peaks:
            self.addPeak(p)
    
    def addPeak(self, peak):
        if not isinstance(peak, Peak):
            raise TypeError(('peak', peak))
        if self._peaks.has_key(peak.id):
            raise ValueError(('peak id', peak))
        if len(peak.dims) != len(self.axes):
            raise ValueError('peak dimensions must match spectral axes')
        self._peaks[peak.id] = peak
    
    def getPeaks(self):
        return self._peaks.values()
    

class Molecule(MyBase):
    
    def __init__(self, residues):
        if not isinstance(residues, dict):
            raise TypeError(('residues', residues))
        for rid in residues.keys():
            if not isinstance(rid, int):
                raise TypeError(('residue id', rid))
        self.residues = residues
    
    
class SpinSystem(MyBase):
    '''
    Not sure about the residueid relationship since I haven't used it yet:
     - what if spinsystem1 is followed by spinsystem2, but they're assigned
       to non-sequential residues?  Should fragments of spinsystems be
       assigned instead?  I can't answer until I actually get this far.
    '''
    
    def __init__(self, ssid, pkids, aatypes, residueids, ssnexts, tags):
        for pkid in pkids:
            if not isinstance(pkid, list):
                raise TypeError(('peak id', pkid))
            if len(pkid) != 2:
                raise ValueError(('peak id', pkid))
            if not isinstance(pkid[0], basestring):
                raise TypeError(('peak id spectrum', pkid))
            if not isinstance(pkid[1], int):
                raise TypeError(('peak id id', pkid))
        self.id = ssid # should be an int
        self.pkids = pkids
        self.aatypes = aatypes
        for rid in residueids:
            if not isinstance(rid, int):
                raise TypeError(('residue id', rid))
        self.residueids = residueids
        for ssid in ssnexts:
            if not isinstance(ssid, int):
                raise TypeError(('spin system id', ssid))
        self.ssnexts = ssnexts
        self.tags = tags
    
    # just a couple example methods that show how I envision instances being used
    #   ssid shouldn't be changed
    def addAAType(self, aatype):
        self.aatypes.append(aatype)
        
    def addTag(self, tag):
        self.tags.append(tag)
        
        
class Project(MyBase):
    
    def __init__(self, name, spectra, molecule, spinsystems):
        '''
        Who knows what 'name' means?
        '''
        if not isinstance(molecule, Molecule):
            raise TypeError(('molecule', molecule))
        self._spectra, self._spinsystems = {}, {}
        for s in spectra:
            self.addSpectrum(s)
        for ss in spinsystems:
            self.addSpinSystem(ss)
        self.name = name
        self.molecule = molecule

    def addSpectrum(self, spectrum):
        if not isinstance(spectrum, Spectrum):
            raise TypeError(('spectrum', spectrum))
        if self._spectra.has_key(spectrum.name):
            raise ValueError(('spectrum name', spectrum))
        self._spectra[spectrum.name] = spectrum
    
    def getSpectra(self):
        return self._spectra.values()
    
    def addSpinSystem(self, ss):
        if not isinstance(ss, SpinSystem):
            raise TypeError(('spin system', ss))
        if self._spinsystems.has_key(ss.id):
            raise ValueError(('spin system id', ss))
        self._spinsystems[ss.id] = ss

    def getSpinSystems(self):
        return self._spinsystems.values()



def buildPeak(spectrum_name, pkid, pk):
    dims = [(d.shift, d.atomtypes) for d in pk.dims]
    return Peak(spectrum_name, pkid, dims, pk.tags, pk.height)


def buildSpectrum(name, spectrum):
    peaks = [buildPeak(name, pkid, pk) for (pkid, pk) in spectrum.peaks.items()]
    return Spectrum(name, spectrum.axes, peaks)


def buildSS(ssid, ss):
    return SpinSystem(ssid, ss.pkids, ss.aatypes, ss.residueids, ss.ssnexts, ss.tags)


def fromModel(project):
    spectra = [buildSpectrum(name, spectrum) for (name, spectrum) in project.spectra.items()]
    spinsystems = [buildSS(ssid, ss) for (ssid, ss) in project.spinsystems.items()]
    return Project(project.name, spectra, Molecule(project.molecule.residues), spinsystems)
