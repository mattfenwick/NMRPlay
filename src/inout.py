import json
from .xeasy import parser
from .xeasy import unparser
from . import xeasyconverter
from . import fromjson


def xeasy_in(projname, paths):
    xpkfls = {}
    for (name, path) in paths.items():
        with open(path, 'r') as infile:
            r = parser.run(parser.xeasy, infile.read())
            if r.status != 'success':
                raise ValueError((name + ' parsing failed', r))
            xpkfls[name] = r.value['result']
    return xeasyconverter.xez2patch(projname, xpkfls)

def xeasy_out(pmodel, paths):
    '''
    Expects a spectrum for each path
    '''
    spectra = xeasyconverter.patch2xez(pmodel)
    for (name, path) in paths.iteritems():
        spectrum = spectra[name]
        with open(path, 'w') as outfile:
            outfile.write(unparser.xeasy(spectrum))

def json_in(jsonpath):
    with open(jsonpath, 'r') as infile:
        return fromjson.loadProject(json.loads(infile.read()))

def json_out(jsonpath, proj):
    with open(jsonpath, 'w') as outfile:
        outfile.write(json.dumps(proj.toJson(), indent=2))
