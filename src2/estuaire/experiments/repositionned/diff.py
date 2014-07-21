import sys

import numpy as np

import agstd.sdb.dbase as dbase

db = dbase.SeismicHDF5DB('/raid/jee/agsis/data/syn_50_new.h5f', 'synthetic')

evtable = [np.load(s) for s in sys.argv[1:]]

original = db.root.catalogs.original.events.read()

c = original[evtable[0]['id']]
d = evtable[0]['position'][:, 0]

spacing = 2.0

mask = [True] * len(evtable)
for i, e in enumerate(evtable):
    for d in range(3):
        mask[i] = (e['position'][:,d] > 0) & (e['position'][:,d] < 127) & mask[i]


print [(len(m), np.sum(m)) for m in mask]

r = [[(original[e['id']][d] - e['position'][:, i] * spacing) ** 2 for i, d in enumerate(['X', 'Y', 'Z'])] for e in evtable]
r2 = [[(original[e[mask]['id']][d] - e[mask]['position'][:, i] * spacing) ** 2 for i, d in enumerate(['X', 'Y', 'Z'])] for e, m in zip(evtable, mask)]


d = np.average(np.sqrt(np.sum(r, axis = 1)), axis = -1)
d2 = np.average(np.sqrt(np.sum(r2, axis = 1)), axis = -1)

for f, g, g2, h, h2 in zip(sys.argv[1:], d, d2, r, r2):
    print f, g, "Masked : ", g2, len(h[0]), len(h2[0])

import matplotlib.pyplot as plt

plt.hist(r2[0][0], bins = 100, range = (0, 100))
plt.xlabel("Positioning Error (m)")
plt.ylabel("Number of Events")
plt.show()



