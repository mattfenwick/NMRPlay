from . import model


fmap_dict = model.fmap_dict

def loadPeakDim(dim):
    return model.PeakDim(dim['shift'], dim['atomtypes'])
        
def loadPeak(pk):
    return model.Peak(map(loadPeakDim, pk['dims']), pk['tags'], float(pk['height']))
    
def loadSpectrum(spec):
    return model.Spectrum(spec['axes'], 
                          dict((int(pkid), loadPeak(pk)) for (pkid, pk) in spec['peaks'].items()))

def loadMolecule(mol):
    return model.Molecule(mol['residues'])

def loadSpinSystem(ss):
    return model.SpinSystem(ss['pkids'], 
                            ss['aatypes'], 
                            ss['residueids'], 
                            ss['ssnexts'], 
                            ss['tags'])

def loadProject(prj):
    '''
    Convert a primitive object graph -- composed of dicts, lists, etc., 
    which are presumably from JSON -- to the model.
    '''
    return model.Project(prj['name'], 
                         fmap_dict(loadSpectrum, prj['spectra']), 
                         loadMolecule(prj['molecule']),
                         dict((int(ssid), loadSpinSystem(ss)) for (ssid, ss) in prj['spinsystems'].items()))
