import datetime
import dateutil
import dateutil.rrule as rrule
import numpy as np

import matplotlib.pyplot as plt

import sys


import agstd.sdb.sqldbase as dbase
import sqlite3

datestart = datetime.date(2009, 1, 1)
dateend = datetime.date(2011, 12, 1)

mqb = dbase.ModelQueryBuilder()
gen = rrule.rrule(rrule.MONTHLY, dtstart = datestart, until = dateend).__iter__()

conn = sqlite3.connect(sys.argv[1])

current = gen.next()
nevents = []
for d in gen:
    mqb.set_date_filter(current, d)
    nevents.append(len(conn.execute(mqb.event_query()).fetchall()))
    current = d



plt.plot(nevents, '-o')
plt.show()
