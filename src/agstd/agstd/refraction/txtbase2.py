import eikonal.data as edata

evdtype = [('id', 'S100'), ('pos', 'float', 3), ('delta_t', 'float')]
stdtype = evdtype

ttdtype = [('stid', 'S100'), ('evid', 'S100'), ('tt', 'float')]

class txtsttable(object):
    def __init__(self, filename):
        self.data = np.loadtxt(filename, dtype = evdtype)
        self.index = create_index(self.data)

    def getTable(self):
        ary = np.array([(i, e['pos'], e['delta_t']) for i, e in
                        enumerate(self.data)],
                       dtype = [('id', 'int'), ('position', 'float', 3),
                                ('delta_t', 'float')])
        return ary

class txttttable(object):
    def __init__(self, filename):
        self.data = np.sort(np.loadtxt(filename, dtype = ttdtype))

    def splitUnique(self, evindex, stindex, evnfile = None, stafile = None):
        uniques = np.unique(self.data['stid'])
        ind = np.searchsorted(self.data['stid'], uniques[1:])
        print uniques, ind
        splitted = np.split(self.data, ind)

        ary = []
        for sid, stt in zip(uniques, splitted):
            print sid, stindex[sid]
            tttable = np.array([(-1, evindex[s['evid']], s['tt'])
                                for s in stt], dtype = [('id', int),
                                                        ('event_id', int),
                                                        ('traveltime', float)])
            ary.append(edata.EKTTTable(tttable, stindex[sid],
                                       evnfile = evnfile, stafile = stafile))
        return ary




def create_index(tt):
    ret = {}
    for i in range(tt.size):
        ret[tt[i]['id']] = i
    return ret

import cPickle as pickle
import numpy as np
import os


def main_txt(input_event, input_station, output_event = None, output_station = None, input_tt = None, output_tt_template = None):
    evdata = txtsttable(input_event)
    event_table = edata.EKEventTable(evdata.getTable())
    if None != output_event:
        pickle.dump(event_table, open(output_event, 'w'), protocol = pickle.HIGHEST_PROTOCOL)

    stdata = txtsttable(input_station)
    station_table = edata.EKStationTable(stdata.getTable())
    if (None != input_station) and (None != output_station):
        pickle.dump(station_table, open(output_station, 'w'), protocol = pickle.HIGHEST_PROTOCOL)

    if (None != input_tt) and (None != output_tt_template):
        ttdata = txttttable(input_tt)
        for tt_table in ttdata.splitUnique(evdata.index, stdata.index, evnfile = os.path.abspath(output_event), stafile = os.path.abspath(output_station)):
            filename = output_tt_template % tt_table.station_id
            pickle.dump(tt_table, open(filename, 'w'), protocol = pickle.HIGHEST_PROTOCOL)



whole_dtype = [ ('event_pos', float, 2),
                ('station_pos', float, 2),
                ('traveltime', float)]


def read_in_one_file(input, output_event = None, output_station = None, output_tt_template = None):
    whole_table = np.loadtxt(input, dtype = whole_dtype)
    evdict = {}
    stdict = {}
    ev_table = []
    st_table = []
    tt_table = {}
    for i, line in enumerate(whole_table):
        evkey = tuple(line['event_pos'])
        stkey = tuple(line['station_pos'])
        if evkey not in evdict:
            evdict[evkey] = len(ev_table)
            ev_table.append((i, evkey + (0.0,), 0.0))
        if stkey not in stdict:
            stdict[stkey] = len(st_table)
            st_table.append((i, stkey + (0.0,), 0.0))
        if stkey not in tt_table:
            tt_table[stkey] = []
        tt_table[stkey].append((i, evdict[evkey], line['traveltime']))

    ev_ary = edata.EKEventTable(np.array(ev_table, dtype = edata.ev_dtype))
    st_ary = edata.EKStationTable(np.array(st_table, dtype = edata.st_dtype))

    if output_event != None:
        pickle.dump(ev_ary, open(output_event, 'w'), protocol = pickle.HIGHEST_PROTOCOL)
    if output_station != None:
        pickle.dump(st_ary, open(output_station, 'w'), protocol = pickle.HIGHEST_PROTOCOL)

    if output_tt_template != None:
        for k, v in tt_table.iteritems():
            filename = output_tt_template % stdict[k]
            tt_ary = edata.EKTTTable(np.array(v, dtype = edata.tt_dtype), stdict[k], evnfile = output_event, stafile = output_station)
            pickle.dump(tt_ary, open(filename, 'w'), protocol = pickle.HIGHEST_PROTOCOL)




if __name__ == "__main__":
    import agstd.main as main
    main.main(read_in_one_file)
    """
    agstd.main.main
    import sys
    import numpy as np
    evdata = txtsttable(sys.argv[1])
    stdata = txtsttable(sys.argv[2])
    ttdata = txttttable(sys.argv[3])

    evtable = edata.EKEventTable(evdata.getTable())
    sttable = edata.EKStationTable(stdata.getTable())

    print evtable.data
    print sttable.data
    print [t.station_id for t in ttdata.splitUnique(evdata.index, stdata.index)]
    """


