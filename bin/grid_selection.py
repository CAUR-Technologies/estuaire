#import sys
#import numpy as np

#import enthought.mayavi.api as mayavi
#import enthought.mayavi.sources.api as msource
#import enthought.mayavi.modules.api as mmodule


#engine = mayavi.Engine()
#engine.start()
#scene = engine.new_scene(name = "Micro-Seismic Model Viewer")
#scene.background = (0.15, 0.15, 0.15)



#for rayfile in sys.argv[1:]:
    #rays = np_load(rayfile)
    #src = msource.VTKDataSource(data = ray.ray_source(rays['rays']))
    #engine.add_source(src)
    #mod = mmodule.Surface()

    #mod.actor.property.opacity = 0.03
    #engine.add_module(mod)

#import wx
#app = wx.App()
#app.MainLoop()

from enthought.traits.api import HasTraits, Instance, on_trait_change
from enthought.traits.api import HasTraits, Instance, Range
from enthought.traits.ui.api import View, Item
from enthought.mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
            SceneEditor

from enthought.mayavi.core.api import Engine

import enthought.mayavi.mlab as mlab

import numpy as np
from agstd.tools import np_load

def extremum(ary):
    return np.min(ary), np.max(ary)

def boundary_selection(events, stations):
    """
    """
    minX, maxX = extremum(events['position'][:,2])
    minX, maxX = minX - (maxX - minX) * 0.1, maxX + (maxX - minX) * 0.1

    minY, maxY = extremum(events['position'][:,1])
    minY, maxY = minY - (maxY - minY) * 0.1, maxY + (maxY - minY) * 0.1

    minZ, maxZ = extremum(events['position'][:,0])
    minZ, maxZ = minZ - (maxZ - minZ) * 0.1, maxZ + (maxZ - minZ) * 0.1

    class MyModel(HasTraits):
        scene = Instance(MlabSceneModel)
        engine = Instance(Engine, args = ())

        originX = Range(low = float(minX), high = float(maxX))
        originY = Range(low = float(minY), high = float(maxY))
        originZ = Range(low = float(minZ), high = float(maxZ))
        shapeX = Range(low = 16, high = 256)
        shapeY = Range(low = 16, high = 256)
        shapeZ = Range(low = 16, high = 256)

        spacing = Range(low = 1, high = 100)

        view = View(Item('scene', height=400, show_label=False,
                        editor=SceneEditor(scene_class=MayaviScene)),
                    Item('originX'),
                    Item('originY'),
                    Item('originZ'),
                    Item('spacing'),
                    Item('shapeX'),
                    Item('shapeY'),
                    Item('shapeZ'), resizable = True)

        def _scene_default(self):
            self.engine.start()
            scene = MlabSceneModel(engine = self.engine)
            return scene

        @on_trait_change('scene.activated')
        def populate_scene(self):
            self.scene.mlab.points3d(events['position'][:,2],
                                     events['position'][:,1],
                                     events['position'][:,0],
                                    scale_factor = 6, color = (0.5, 0.5, 0.0))
            self.scene.mlab.points3d(stations['position'][:,2],
                                     stations['position'][:,1],
                                     stations['position'][:,0],
                                    scale_factor = 10, color = (0.5, 0.0, 0.0))
            extent = [self.originX, self.originX + self.shapeX * self.spacing, 
                      self.originY, self.originY + self.shapeY * self.spacing,
                      self.originZ, self.originZ + self.shapeZ * self.spacing]
            self.outline = mlab.outline(extent = extent)

        @on_trait_change("originX, originY, originZ, shapeX, shapeY, shapeZ, spacing")
        def modif_outline(self):
            extent = [self.originX, self.originX + self.shapeX * self.spacing, 
                      self.originY, self.originY + self.shapeY * self.spacing,
                      self.originZ, self.originZ + self.shapeZ * self.spacing]

            self.outline.bounds = extent


    viewer = MyModel()
    viewer.configure_traits()


if __name__ == '__main__':
    import agstd.main as main
    main.main(boundary_selection, events = np_load, stations = np.load)
