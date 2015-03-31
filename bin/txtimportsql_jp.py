#!/usr/bin/env python
import eikonal.data as edata
from IPython.core.debugger import Tracer

evdtype = [('id', 'S100'), ('pos', 'float', 3), ('delta_t', 'float')]
stdtype = evdtype

ttdtype = [('stid', 'S100'), ('evid', 'S100'), ('tt', 'float')]

def create_index(tt):
    ret = {}
    for i in range(tt.size):
        ret[tt[i]['id']] = i
    return ret

import cPickle as pickle
import numpy as np
import os
import datetime


import sqlite3
import agstd.sdb.sqldbase as dbase



whole_dtype = [('event_id', int),
                ('station_name','|S20'),
                ('station_pos', float, 3),
                ('event_pos', float, 3),
                ('traveltime', float),
                ('date', '|S50'),
                ('type', '|S20')]

def show_mlab(events, stations):
    from mayavi import mlab
    x, y, z = events.data['position'][:,0], events.data['position'][:,1], events.data['position'][:,2]
    print x.min(), x.max(), x.size,  y.min(), y.max(), y.size, z.min(), z.max(), z.size
    mlab.points3d(x, y, z, mode = 'cube', color = (0, 0, 1) , scale_factor = 10)

    x, y, z = stations.data['position'][:,0], stations.data['position'][:,1], stations.data['position'][:,2]
    mlab.points3d(x, y, z, mode = 'sphere', color = (0.5, 0.5, 0) , scale_factor = 20)
    mlab.show()


seismogram_query = \
"""
INSERT INTO seismogram (station_id, event_id, date, frame, sample_rate) VALUES(:staid, :evnid, :date, :frame, :sample_rate)
"""

event_query = \
"""
INSERT INTO event(id,X,Y,Z,Moment_P,Moment_S,Energy_P,Energy_S,origin)
    VALUES(:id,:X,:Y,:Z,:P_Mo,:S_Mo,:P_E,:S_E,:filename)
"""

delete_event_query = \
"""
DELETE FROM event WHERE id = :id;
"""

catalog_event_query = \
"""
INSERT INTO original_event_catalog(event_id,date,X,Y,Z)
    VALUES(:id,:date,:X,:Y,:Z);
"""

catalog_traveltime_query = \
"""
INSERT INTO original_traveltime_catalog(pick_id,traveltime)
                                        VALUES(?,?);
"""

pick_query = \
"""
INSERT INTO pick(seismogram_id,type,date,frame) VALUES(?,?,?,?);
"""


def add_station(db, stname):
    cursor = db.execute("SELECT id FROM station WHERE name = :name", stname)
    staid = cursor.fetchone()
    if staid is None:
        db.execute("INSERT INTO station(name, X, Y, Z) VALUES(:name,:X,:Y,:Z)", stname)
        staid = db.execute("SELECT last_insert_rowid() FROM station").fetchone()

    return int(staid[0])

def add_event(db, evkey, pos, discrepencies, date, lineid = None):
    evdict = {"id" : int(evkey),
              "X" : pos[0],
              "Y" : pos[1],
              "Z" : pos[2],
              "P_Mo" : 0,
              "S_Mo" : 0,
              "P_E" : 0,
              "S_E" : 0,
              "filename" : "meh",
              "seismogram_id" : int(evkey),
              "type" : "P",
              "date" : datetime.datetime.strptime(date, "%Y/%m/%d %H:%M:%S.%f"),
              "frame" : 0,
              "datetime" : 0}
    cursor = db.execute("SELECT id FROM event WHERE id = :id AND X = :X AND Y = :Y AND Z = :Z", evdict )
    evid = cursor.fetchone()
    if evid == None:
        try:
            db.execute(event_query, evdict)
            db.execute(catalog_event_query, evdict)
        except:
            if evdict['id'] not in discrepencies:
                discrepencies[evdict['id']] = []
            discrepencies[evdict['id']].append(lineid)

            cursor = db.execute("SELECT X, Y, Z FROM event WHERE id = :id", evdict )
            print "Event position incongruency detected on line ", lineid, " [", evdict['id'], "]", pos, " in DB ", cursor.fetchone()
            #print ""

    return int(evkey)

def add_tt(db, evkey, stkey, tt, catalog, ptype):
    sedict = {"evnid" : int(evkey),
              "staid" : int(stkey),
              "tt" : tt,
              "catalog" : catalog,
              "ptype" : ptype,
              "date" : 0,
              "frame" : 0,
              "sample_rate" : 0}
    cursor = db.execute("SELECT id from seismogram WHERE event_id = :evnid AND station_id = :staid", sedict)
    seismogram_id = cursor.fetchone()
    if seismogram_id == None:
        db.execute(seismogram_query, sedict)
        seismogram_id = db.execute("SELECT last_insert_rowid() FROM seismogram").fetchone()

    seismogram_id = seismogram_id[0]

    db.execute(pick_query, (seismogram_id, ptype, 0, 0))
    pick_id = db.execute("SELECT last_insert_rowid() FROM pick").fetchone()[0]


    db.execute(catalog_traveltime_query, (pick_id, tt))



def setdb(filename):
    if os.path.exists(filename):
        os.remove(filename)
    db = sqlite3.connect(filename)
    db.executescript(dbase.__DEFAULT_SEISMIC_DB__)

    return db



def read_in_one_file(input, sqloutput, output_event = None, output_station = None, output_tt_template = None, show_result = False):
    db = setdb(sqloutput)

    discrepencies = {}

    whole_table = np.loadtxt(input, dtype = whole_dtype, delimiter = ',')
    evdict = {}
    stdict = {}
    ev_table = []
    st_table = []
    tt_table = {}
    stname = {}
    for i, line in enumerate(whole_table):
        evkey = line['event_id']
        stkey = tuple(line['station_pos'])
        if evkey not in evdict:
            evdict[evkey] = len(ev_table)
            ev_table.append((evkey, line['event_pos'], 0.0))
        if stkey not in stdict:
            stdict[stkey] = len(st_table)
            st_table.append((i, stkey, 0.0))
        stname['name'] = line['station_name']
        stname['X'] = line['station_pos'][0]
        stname['Y'] = line['station_pos'][1]
        stname['Z'] = line['station_pos'][2]

        db_stid = add_station(db, stname)
        db_evid = add_event(db, evkey, line['event_pos'], discrepencies, line['date'], lineid = i)
        tt_id = add_tt(db, db_evid, db_stid, line['traveltime'], "txt", line['type'])

        if stkey not in tt_table:
            tt_table[stkey] = []
        tt_table[stkey].append((i, evdict[evkey], line['traveltime']))

    ev_ary = edata.EKEventTable(np.array(ev_table, dtype = edata.ev_dtype))
    st_ary = edata.EKStationTable(np.array(st_table, dtype = edata.st_dtype))

    print "Results"
    if output_event != None:
        pickle.dump(ev_ary, open(output_event, 'w'), protocol = pickle.HIGHEST_PROTOCOL)
        print "\tNumber of Events : %i (%d in SQL)" % (ev_ary.size, db.execute("SELECT COUNT(id) FROM event").fetchone()[0])
    if output_station != None:
        pickle.dump(st_ary, open(output_station, 'w'), protocol = pickle.HIGHEST_PROTOCOL)
        print "\tNumber of Stations : %i (%d in SQL)" % (st_ary.size, db.execute("SELECT COUNT(id) FROM station").fetchone()[0])

    mb = dbase.ModelQueryBuilder()
    if output_tt_template != None:
        mb.set_type_filter("P")
        Plen = len(db.execute(mb.traveltime_query()).fetchall())
        mb.set_type_filter("S")
        Slen = len(db.execute(mb.traveltime_query()).fetchall())
        print "\tTotal Recording : %d P and %d S" % (Plen, Slen)
        print "\tTraveltimes :"
        mb.set_type_filter(None)
        for k, v in tt_table.iteritems():
            ttary = np.array(v, dtype = edata.tt_dtype)
            mb.set_inverse_station_filter([stdict[k] + 1])
            sqllen = len(db.execute(mb.traveltime_query()).fetchall())
            print "\t\tStation ID %i - %i recording (S + P) (%i in SQL)" % (stdict[k] + 1, ttary.size, sqllen)
            filename = output_tt_template % stdict[k]
            tt_ary = edata.EKTTTable(ttary, stdict[k], evnfile = output_event, stafile = output_station)
            pickle.dump(tt_ary, open(filename, 'w'), protocol = pickle.HIGHEST_PROTOCOL)

    print "Discrepencies"
    for k, n in discrepencies.iteritems():
        print "\tEvent [%d] : " % (k), n

    if show_result:
        show_mlab(ev_ary, st_ary)

    db.commit()





if __name__ == "__main__":
    import agstd.main as main
    main.main(read_in_one_file, show_result = bool)


