import numpy as np

import agstd.sdb.sqldbase as dbase
import agstd.main as main
import sqlite3


def separate(db, nbins):
    conn = sqlite3.connect(db)
    dates = conn.execute("SELECT date FROM original_event_catalog ORDER BY date").fetchall()
    splitted = np.split(dates, range(0, len(dates), len(dates) / nbins))
    for s in splitted[1:]:
        print s[0]

    print "Number of elements per bin : %d" % (len(dates) / nbins)

if __name__ == '__main__':
    main.main(separate, nbins = int)
