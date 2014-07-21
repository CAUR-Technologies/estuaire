import os, sys

import enthought.mayavi.core.api as mcore
import enthought.mayavi.sources.api as msources
import enthought.mayavi.modules.api as mmodules
import enthought.mayavi.filters.api as mfilters

import lxml.objectify as objectify

def __build_object__(obj, element):
    subobj = getattr(obj, element.tag)
    try:
        for prop in element.attrib:
            setattr(subobj, prop, getattr(subobj, prop).__class__(element.attrib[prop]))
        for i in element.iterchildren():
            __build_object__(subobj, i)
    except AttributeError, e:
            __build_object__(subobj, element)
    except Exception, e:
        print e

def __build_module_manager__(element, engine, parent = None):
    mmobj = mcore.ModuleManager()
    engine.add_module(mmobj, parent)
    childs = []
    for t in ['Module']:
        for child in element.findall(t):
            childs.append(build_type[t](child, engine, parent = mmobj))
    for i in element.iterchildren():
        if i.tag not in ['Module']:
            print i.tag
            __build_object__(mmobj, i)

    return (mmobj, childs)

def __build_module__(element, engine, parent = None):
    moddesc = element.attrib
    mtype = moddesc.pop("type")
    modobj = getattr(mmodules, mtype)(**moddesc)
    engine.add_module(modobj, parent)

    for i in element.iterchildren():
        __build_object__(modobj, i)

    return modobj

def __build_filter__(element, engine, parent = None):
    filtdesc = element.attrib
    ftype = moddesc.pop("type")
    filtobj = getattr(mfilters, ftype)(**filtdesc)
    engine.add_filter(filtobj, parent)

    for i in element.iterchildren():
        __build_object__(modobj, i)

def __build_scene__(element, engine, parent = None, dirname = None):

    # Extracting the non trivial description of the source
    description = element.attrib
    srctype = description.pop("type")
    srcfile = description.pop("file") if dirname is None else os.path.join(dirname, description.pop("file"))

    # Instanciating the trivial source
    srcobj = getattr(msources, srctype)(**description)

    # Uncomment in the realease version
    srcobj.initialize(srcfile)

    # Adding the source to the current engine
    engine.add_source(srcobj, parent)

    # Processing the rest of the pipeline
    childs = []
    for t in ['Source', 'Module', 'Filter', 'ModuleManager']:
        for child in element.findall(t):
            childs.append(build_type[t](child, engine, parent = srcobj))
    return (srcobj, childs)

build_type = {'Scene' : __build_scene__,
              'ModuleManager' : __build_module_manager__,
              'Module' : __build_module__,
              'Filter' : __build_filter__}


class MayaviPipeline(object):
    def __init__(self, xmlfile, engine, scene, srcdir = None):
        root = objectify.XML(xmlfile)
        self.root = [__build_scene__(e, engine, parent = scene, dirname = srcdir) \
                     for e in root.findall('Source')]



if __name__ == '__main__':
    import enthought.mayavi.api as mayavi
    engine = mayavi.Engine()
    engine.start()
    scene = engine.new_scene(name = "Test")
    pipeline = MayaviPipeline(open('pipeline_examples/simple_surface.xml').read(), engine, scene, srcdir = '/home/jee/work/AGSIS-local/E26/all_data/vtk/d7_2004-08-01/')
