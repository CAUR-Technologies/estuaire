import sys

import numpy as np

import agstd.sdb.dbase as dbase

db = dbase.SeismicHDF5DB('/raid/jee/agsis/data/syn_50_new.h5f', 'synthetic', 'a', sanity_check = False)


stdev = float(sys.argv[1])


original = db.root.catalogs.original.events.read()
ott = db.root.catalogs.original.traveltimes.read()

catalog = db.new_catalog("random", override = True)

catalog.traveltimes.append(ott)

evtable = original

for e in ['X', 'Y', 'Z']:
    evtable[e] += np.random.normal(scale = stdev, size = evtable.size)

catalog.events.append(evtable)

catalog.events.flush()
catalog.traveltimes.flush()
    
