from enthought.tvtk.tvtk import tvtk
import sys

imp = tvtk.VRMLImporter()
imp.file_name = sys.argv[1]
imp.read
imp.update()

actors = imp.renderer.actors
pfilter = tvtk.AppendPolyData()
for i, a in enumerate(actors):
    pfilter.add_input(a.mapper.input)

writer = tvtk.XMLPolyDataWriter()
writer.input = pfilter.output
writer.file_name = "%d.xml" % i
writer.write()


