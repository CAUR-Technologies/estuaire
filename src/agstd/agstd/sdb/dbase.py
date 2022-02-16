import numpy as np
import scipy as sc
import scipy.sparse

import tables

import  datetime, \
        time, \
        dateutil.parser

from itertools import izip


import logging
log = logging.getLogger(__name__)

def tounixtime(date):
    if isinstance(date, datetime.datetime) or isinstance(date, datetime.date):
        return time.mktime(date.timetuple())
    elif isinstance(date, int):
        return date
    elif isinstance(date, str):
        return time.mktime(dateutil.parser.parse(date).timetuple())
    else:
        raise ValueError("Cannot be converted to unix time")

class pick_type:
    S_Wave, P_Wave = range(2)
    strdict = {"P_Wave" : P_Wave,
               "S_Wave" : S_Wave,
               "P-Wave" : P_Wave,
               "S-Wave" : S_Wave,
               "P"      : P_Wave,
               "S"      : S_Wave}

class SeismogramDescription(tables.IsDescription):
    """
    This object represent the datatype of a seismogram. We limit the size
    of the array to 3x20000. Since this table take almost 90% of the database
    size, this table should be highly compressed to ensure a footprint as
    small as possible.
    """
    seisX = tables.Float64Col(shape = 20000, pos = 0)
    seisY = tables.Float64Col(shape = 20000, pos = 1)
    seisZ = tables.Float64Col(shape = 20000, pos = 2)
    station_id = tables.Int32Col(pos = 3) # Row position
    event_id = tables.Int32Col(pos = 4) # Row position
    sampling_rate = tables.Float64Col(pos = 5)
    
    length = tables.Int32Col(shape = 3, pos = 6)


class TriggerDescription(tables.IsDescription):
    """
    This table represent a basic trigger description. It contains the id of
    the station as well as the event. Also contain the sample rate, the
    Trigger Time and the Position of the trigger in the Seismogram. This
    can also contain the original filename.
    The Field Hypo_distance is for viewing purpose only. Prefer the
    calculation from this value since it strongly depends on the catalog hypo
    center positionning.
    """
    station_id          = tables.Int32Col(pos = 0) # Row position
    event_id            = tables.Int32Col(pos = 1) # Row position

    sample_rate         = tables.Float64Col(pos = 2)
    trigger_time        = tables.Float64Col(pos = 3)
    trigger_frame       = tables.Int32Col(pos = 4)

    Q0                  = tables.Float64Col(pos = 5)
    Q                   = tables.Float64Col(pos = 6)
    natural_frequency   = tables.Float64Col(pos = 7)
    damping_factor      = tables.Float64Col(pos = 8)
    meter_type          = tables.StringCol(32, pos = 9)

    filename            = tables.StringCol(64, pos = 10)

    orig_hypo_distance  = tables.Float64Col(pos = 11)


class PickDescription(tables.IsDescription):
    """
    This table represent a pick in the seismogram. It contains the velocity
    as well as the time of the arrival time as the absolute number of seconds
    pass in the current day. The date can be found in the coresponding event
    table.
    """
    trigger_id          = tables.Int32Col(pos = 0) # Row position

    velocity            = tables.Float64Col(pos = 1)
    time                = tables.Float64Col(pos = 2)

    frame               = tables.Int32Col(pos = 3)

    pick_type           = tables.Int32Col(pos = 4)


class SEventDescription(tables.IsDescription):
    """
    This class represent the basic seismic event identity. It is present only
    to keep the date of the event as well as the original filename from which
    it was imported.

    The fields orig_{X} are present if the original import already contains
    a calculated position and a time. Prefer the position included in the
    the coresponding catalog.
    """
    date                = tables.Int64Col(pos = 0) # IN POSIX TIME
    filename            = tables.StringCol(64, pos = 1)

    Mo_P                = tables.Float64Col(pos = 2)
    E_P                 = tables.Float64Col(pos = 3)
    f_P                 = tables.Float64Col(pos = 4)

    Mo_S                = tables.Float64Col(pos = 5)
    E_S                 = tables.Float64Col(pos = 6)
    f_S                 = tables.Float64Col(pos = 7)

    orig_X              = tables.Float64Col(pos = 8)
    orig_Y              = tables.Float64Col(pos = 9)
    orig_Z              = tables.Float64Col(pos = 10)

    orig_time           = tables.Float64Col(pos = 11)


class StationDescription(tables.IsDescription):
    """
    Simple Seismic Station Desctiption. Only contains Position.
    """
    X                   = tables.Float64Col(pos = 0)
    Y                   = tables.Float64Col(pos = 1)
    Z                   = tables.Float64Col(pos = 2)
    name                = tables.StringCol(32, pos = 3)


class CatalogTraveltimeDescription(tables.IsDescription):
    """
    Simple traveltime table.

    This table should have the same length as the pick tables from the main
    database since it directly refers to the record at the same indice.
    """
    traveltime          = tables.Float64Col(pos = 0)


class CatalogSEventDescription(tables.IsDescription):
    """
    This table contains the event description for a catalog. It contains
    the calculated field from a particular algorithm associated with the
    catalog. For each event, we also include the number of P/S-wave arrival
    from which the position/time were calculated. A field for the residual
    as well as an error flag and a message can also be included.

    This table should be the same length as the event table since it refers
    directly to the coresponding indice in that table.
    """
    X                   = tables.Float64Col(pos = 0)
    Y                   = tables.Float64Col(pos = 1)
    Z                   = tables.Float64Col(pos = 2)

    time                = tables.Float64Col(pos = 3)
    residual            = tables.Float64Col(pos = 4)

    S_observations      = tables.Int32Col(pos = 5)
    P_observations      = tables.Int32Col(pos = 6)

    error               = tables.BoolCol(pos = 7)
    message             = tables.StringCol(64, pos = 8)


class CatalogProblem(tables.IsDescription):
    """
    This table contains a list of trigger_id to use in the a tomography
    problem. This only have a list of trigger that participate in the
    modelization process.
    """
    pick_id          = tables.Int32Col(pos = 0)


class SeismicHDF5DB(object):
    """
    This object manage the Micro-Seismic Event Database With an underlying
    HDF5 File backend.
    """
    def __init__(self, filename, groupname,
                 mode = 'r', description = "Seismic DB",
                 sanity_check = True, override = False):
        # Pretty logging
        if log.isEnabledFor(logging.INFO):
            log.info("Opening %s HDF5 Seismic Event Database [group : %s]" % \
                     (filename, groupname))

        self.h5f = tables.openFile(filename, mode)
        self.rootgroup = groupname
        if (override is True) and (mode == 'r'):
            raise ValueError("Cannot override group oppened in read-only mode")
        if (override is True) and hasattr(self.h5f.root, groupname):
            self.h5f.removeNode('/', groupname, recursive = True)
        if not hasattr(self.h5f.root, groupname):
            if log.isEnabledFor(logging.INFO):
                log.info("Creating the initial group hierarchy")

            # Creating Basic Root Nodes
            self.root = self.h5f.createGroup('/', groupname, description)
            # Creating std Tables
            self.h5f.createTable(self.root , "stations",
                                 StationDescription, "Station description")
            self.h5f.createTable(self.root , "events",
                                 SEventDescription, "Event description")

            # Creating the catalog group
            self.h5f.createGroup(self.root, "catalogs", "Catalog Arrival Time")

            # The seismogram table is highly compressed with zlib algorithm for
            # the sake of hfd5 compatibilty
            sfilter = tables.Filters(complevel=9,
                                     complib='zlib', fletcher32=True)
            self.h5f.createTable(self.root, "seismograms",
                                 SeismogramDescription,
                                 "Seismogram", filters = sfilter)

            # The triggers Groups in the catalog
            triggers = self.h5f.createGroup(self.root,
                                            "triggers", "Trigger Data")
            self.h5f.createTable(triggers, "description",
                                 TriggerDescription, "Trigger description")
            self.h5f.createTable(triggers, "pick",
                                 PickDescription, "P/S-Wave Trigger")
        else:
            self.root = getattr(self.h5f.root, groupname)
            if sanity_check:
                self.sanity_check()


    def get_seismogram(self, trigger_id):
        """
        This method return the seismogram assiciated with the coresponding
        trigger
        """
        full, length = self.root.seismograms[trigger_id]
        return full[:length]

    def event_iterator(self, ptypes = [pick_type.S_Wave, pick_type.P_Wave]):
        """
        This method iterate over the events and returns the event ID, the list
        of triggers
        """
        if not hasattr(ptypes, "__iter__"):
            ptypes = [ptypes]
        evt_trigger_table = self.root.triggers.description.col('event_id')
        pick_table = self.root.triggers.pick.col('trigger_id')
        ptype_table = self.root.triggers.pick.col('pick_type')
        for evt_id in range(self.root.events.nrows):
            trg_ids = range(*np.searchsorted(evt_trigger_table, [evt_id, evt_id + 1]))
            for i, t in enumerate(trg_ids):
                trg_ids[i] = t, [pid for pid in range(*np.searchsorted(pick_table, [t, t+1])) if ptype_table[pid] in ptypes]
            yield evt_id, trg_ids

    def sanity_check(self):
        """
        TODO : Should check for unsorted event_id, sanity of catalog data, one-to-many
        relations, etc...
        """


        for c in self.root.catalogs:
            if c.events.nrows != self.root.events.nrows:
                if log.isEnabledFor(logging.ERROR):
                    log.error("Sanity check failed  [catalog : %s] number of events row mismatch")
                raise RuntimeError("Catalog Inconsistency detected")
            if c.traveltimes.nrows != self.root.triggers.pick.nrows:
                if log.isEnabledFor(logging.ERROR):
                    log.error("Sanity check failed  [catalog : %s] number of traveltime row mismatch")
                raise RuntimeError("Catalog Inconsistency detected")

        if log.isEnabledFor(logging.INFO):
            log.info("Sanity Check Passed for <%s, %s> %s" % \
                    (self.h5f.filename, self.root._v_name, self.root._g_gettitle()))

    def sort_trigger(self, keys = ['event_id', 'station_id']):
        if self.h5f.mode == 'r':
            raise AttributeError("Database is oppened in Read-Only")

        t_table = self.root.triggers.description.read()
        p_table = self.root.triggers.pick.read()


        t_order = np.argsort(t_table, order = keys)
        t_table = t_table[t_order]

        i_order = np.argsort(t_order)

        p_table['trigger_id'] = i_order[p_table['trigger_id']]

        p_order = np.argsort(p_table, order = ['trigger_id'])
        p_table = p_table[p_order]

        self.h5f.enableUndo()
        self.h5f.mark()

        try:
            self.root.triggers.description.modifyRows(rows = t_table)
            self.root.triggers.pick.modifyRows(rows = p_table)

            for c in self.root.catalogs:
                tt = c.traveltimes.read()[p_order]
                c.traveltimes.modifyRows(rows = tt)

        except:
            self.h5f.undo()
        finally:
            self.h5f.disableUndo()



    def new_catalog(self, cname, override = False, description = "Catalog"):
        """
        This method create a simple catalog with the given description and
        create the basic table hierarchy.

        If another catalog with the same name exists, it raises an Exception
        unless the parameter override is True.
        """
        if hasattr(self.root.catalogs, cname):
            if override:
                self.h5f.removeNode(self.root.catalogs,
                                    name = cname, recursive = True)
            else:
                raise AttributeError("Catalog already exists")

        catalog = self.h5f.createGroup(self.root.catalogs,
                                       cname, "Catalog Description")
        self.h5f.createTable(catalog, 'events',
                             CatalogSEventDescription, "Catalog Event Description")

        self.h5f.createTable(catalog, "traveltimes",
                             CatalogTraveltimeDescription, "Catalog Traveltimes")

        return catalog

    def add_station(self, position):
        """
        This Method add a station in the database. It ensure the unique existance
        of a station.

        Returns the station id (the index in the table)
        """
        sTable = self.root.stations.read()
        index = np.where((sTable['X'] == position[0]) & \
                         (sTable['Y'] == position[1]) & \
                         (sTable['Z'] == position[2]))[0]
        if len(index) > 1:
            raise RuntimeError("Station Table Corruption Detected ... ")
        elif len(index) == 0:     # The station do not exists
            index = sTable.size
            self.root.stations.append([position])
            self.root.stations.flush()
        else:
            index = index[0]

        return index

    def __del__(self):
        if hasattr(self, 'h5f'):
            self.h5f.close()


class Index(object):
    def __init__(self, indexname):
        pass

class ModelArray(object):
    """
    This object represent a subset of the original database and make it
    easier to filter and only use part of it.
    """
    def __get_mask__(self):
        return self.__mask__

    def __set_mask__(self, mask):
        self.__mask__ = mask
        self.unique_events = np.unique(self.col_evtid[self.mask])
        self.unique_events.sort()
        self.unique_stations = np.unique(self.col_staid[self.mask])
        self.unique_stations.sort()
        self.size = np.sum(self.__mask__ & self.pmask)

    mask = property(__get_mask__, __set_mask__)


    def __init__(self, seismicdb, catalog_name, pick_ids = None, ptype = pick_type.P_Wave):
        self.seismicdb = seismicdb
        self.catalog_name = catalog_name

        if pick_ids is None:
            pick_ids = np.arange(seismicdb.root.triggers.pick.nrows)

        self.type_mask = { pick_type.P_Wave : seismicdb.root.triggers.pick.col('pick_type') == pick_type.P_Wave,
                           pick_type.S_Wave : seismicdb.root.triggers.pick.col('pick_type') == pick_type.S_Wave}

        self.__ptype__ = pick_type.strdict[ptype] if isinstance(ptype, str) else ptype
        self.pmask = self.type_mask[self.__ptype__]

        self.pick_ids = pick_ids
        self.pick_ids = pick_ids

        self.__updatedb__()
        self.reset_mask()


    def __set_ptype__(self, ptype):
        if isinstance(ptype, str):
            self.__ptype__ = pick_type.strdict[ptype]
            self.pmask = self.type_mask[self.__ptype__]
        else:
            self.pmask = self.type_mask[ptype]
            self.__ptype__ = ptype
        self.size = np.sum(self.__mask__ & self.pmask)

    def __get_ptype__(self):
        return self.__ptype__
    ptype = property(__get_ptype__, __set_ptype__)

    def __get_masked_arys__(self):
        mask = self.pmask & self.mask
        return [m[mask] for m in \
                (self.col_staid, self.col_evtid, self.col_triid, self.traveltimes, np.arange(self.mask.size, dtype = 'int32'))]
    masked_arys = property(__get_masked_arys__)

    def __get_sorted_arys__(self):
        """
        Sort By Station
        """
        masked_arys = self.__get_masked_arys__()
        sorted_index = np.argsort(masked_arys[0])
        return [m[sorted_index] for m in masked_arys]
    sorted_arys = property(__get_sorted_arys__)


    def split_stations(self, ptype = pick_type.P_Wave):
        """
        """
        sta, evt, tri, tt, pick = self.__get_sorted_arys__()

        # Createion of the double difference Q matrix

        split_indices = np.searchsorted(sta, self.unique_stations, side = 'right')
        res = [self.unique_stations] + [np.split(s, split_indices) for s in [evt, tri, tt, pick]]
        return res

    def __updatedb__(self):
        # Reload the pick column, stations
        self.col_triid = self.seismicdb.root.triggers.pick.col('trigger_id')[self.pick_ids]
        self.col_evtid = self.seismicdb.root.triggers.description.col('event_id')[self.col_triid]
        self.col_staid = self.seismicdb.root.triggers.description.col('station_id')[self.col_triid]

        self.stations = self.seismicdb.root.stations.read()

        # Reload the catalog
        print 'voici le nom du catalog', self.catalog_name
        self.catalog = getattr(self.seismicdb.root.catalogs, self.catalog_name)
        self.traveltimes = self.catalog.traveltimes.col('traveltime')[self.pick_ids]

        self.events = self.catalog.events.read()

    def reset_mask(self, include_error = False):
        """
        This method reset pick mask to its original value. (i.e. only remove
        negative traveltimes and error in the events positionning)
        """
        if include_error:
            mask = np.ones(self.traveltimes.size, dtype = 'bool')
            if log.isEnabledFor(logging.WARNING):
                log.warning("Model Could Contain Errors and/or Negative traveltimes, Use with Care")
        else:
            inverse_mask = self.events['error'][self.col_evtid] \
                    | (self.traveltimes < 0) \
                    | np.isnan(self.traveltimes)

            if log.isEnabledFor(logging.INFO):
                log.info("Reseting Mask [%d masked]" % (np.sum(inverse_mask)))
            mask = np.where(inverse_mask, False, True)

        self.mask = mask

    def __get_station_centroid__(self):
        desc = self.stations[self.unique_stations]
        return [desc[col].mean() for col in ('X', 'Y', 'Z')]
    station_centroid = property(__get_station_centroid__)

    def ball(self):
        desc = self.stations[self.unique_stations]
        centroid = [desc[col].mean() for col in ('X', 'Y', 'Z')]
        dmax = [np.abs(desc[col] - c).max() for col, c in izip(('X', 'Y', 'Z'), centroid)]

        return centroid, max(dmax)

    def add_centroid_mask(self, radius_factor = 1):
        """
        This method add a centroid mask to the current object mask. The center
        is the centroid of the station and the radius is exprimed in the term
        of the greatest distance between the centroid and a station.
        """
        centroid, radius = self.ball()
        events = self.events[self.col_evtid]
        mask = self.mask
        for c, col in izip(centroid, ('X', 'Y', 'Z')):
            mask &= (events[col] - c) ** 2 < (radius_factor * radius) ** 2
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Station Centroid Mask <factor %f> [%d masked]" % \
                     (radius_factor, nmasked))
        self.mask = mask
        return nmasked

    def add_station_mask(self, stations):
        """
        This method mask all the triggers for a given list station/or station.
        """
        if not hasattr(stations, "__iter__"):
            stations = [stations]
        mask = self.mask
        for s in stations:
            mask = mask & (self.col_staid != s)
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Station Mask for %d stations [%d masked]" % \
                     (len(stations), nmasked))
        self.mask = mask
        return nmasked

    def add_event_mask(self, events):
        """
        This method add en 
        """
        if not hasattr(events, "__iter__"):
            events = [events]
        mask = self.mask
        for e in events:
            mask &= self.col_evtid != e
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Event Mask for %d events [%d masked]" % \
                     (len(events), nmasked))
        self.mask = mask
        return nmasked

    def add_cuboid_mask_zyx(self, origin, length, padding = 0):
        """
        """
        mask = self.mask
        events = self.events[self.col_evtid]
        stations = self.stations[self.col_staid]
        for o, c, l in izip(origin, ('Z', 'Y', 'X'), length):
            mask &= (events[c] > (o +  l * padding)) & (events[c] < (o + l * (1 - padding)))
            mask &= (stations[c] > (o + l * padding)) & (stations[c] < (o + l * (1 - padding)))
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Cuboid Mask for <origin %f, %f, %f> <size %f, %f, %f> [%d masked]" % \
                     (origin[0], origin[1], origin[2], length[0], length[1], length[2], nmasked))
        self.mask = mask
        return nmasked

    def add_cuboid_mask(self, origin, length, padding = 0):
        origin = [o for o in reversed(origin)]
        length = [l for l in reversed(length)]
        return self.add_cuboid_mask_zyx(origin, length, padding = padding)

    def add_date_mask(self, lower_date, higher_date):
        """
        """
        ulower = tounixtime(lower_date)
        uupper = tounixtime(higher_date)
        evtdate = self.seismicdb.root.events.col('date')[self.col_evtid]
        mask = self.mask & (evtdate >= ulower) & (evtdate <= uupper)
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Date Mask from %s to %s [%d masked]" % \
                     (str(datetime.datetime.fromtimestamp(ulower)),
                      str(datetime.datetime.fromtimestamp(uupper)),
                      nmasked))
        self.mask = mask
        return nmasked

    def add_outlier_mask(self, stdmax):
        """
        """
        masked = self.traveltimes[self.mask]

        mean = np.average(self.mask)
        stddev = np.std(self.mask)

        mask = self.mask & ((self.traveltimes - mean) ** 2 < stdmax * (stddev ** 2))
        nmasked = np.sum(self.mask) - np.sum(mask)
        if log.isEnabledFor(logging.INFO):
            log.info("Applied Outlier Mask for traveltime over %f stddev [%d masked]" % \
                     (stdmax, nmasked))
        self.mask = mask
        return nmasked

    def __get_unique_events_positions__(self):
        evtids = self.unique_events
        desc = self.events[evtids]
        return np.array([desc[c] for c in ['Z', 'Y', 'X']], dtype = 'float64')
    unique_events_positions = property(__get_unique_events_positions__)

    def __get_unique_stations_positions__(self):
        staids = self.unique_stations
        desc = self.stations[staids]
        return np.array([desc[c] for c in ['Z', 'Y', 'X']], dtype = 'float64')
    unique_stations_positions = property(__get_unique_stations_positions__)

    def iter_traveltimes(self):
        for sindex, (staid, evtlst, trglst, ttlst, picklst) in enumerate(izip(*self.split_stations())):
            for evtid, trgid, tt, pickid in izip(evtlst, trglst, ttlst, picklst):
                eindex = np.searchsorted(self.unique_events, evtid)
                yield (staid, sindex, evtid, eindex, trgid, tt, pickid)






