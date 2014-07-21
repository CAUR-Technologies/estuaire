import sqlite3
import string

__DEFAULT_SEISMIC_DB__ = \
"""
PRAGMA foreign_keys = ON;
CREATE TABLE station
    (name       VARCHAR(64) DEFAULT "unamed station",
     id         INTEGER     PRIMARY KEY,
     X          DOUBLE       NOT NULL,
     Y          DOUBLE       NOT NULL,
     Z          DOUBLE       NOT NULL,
     CONSTRAINT uc_position UNIQUE (X, Y, Z));


CREATE TABLE event
    (name       VARCHAR(64) DEFAULT "unnamed event",
    id          INTEGER     PRIMARY KEY,
    X           DOUBLE       NOT NULL,
    Y           DOUBLE       NOT NULL,
    Z           DOUBLE       NOT NULL,
    Moment_P    DOUBLE,
    Energy_P    DOUBLE,
    Moment_S    DOUBLE,
    Energy_S    DOUBLE,
    origin      VARCHAR(64));


CREATE TABLE seismogram
    (id                 INTEGER     PRIMARY KEY,
     event_id           INTEGER     NOT NULL,
     station_id         INTEGER     NOT NULL,
     date               DATE        NOT NULL,
     frame              INTEGER     NOT NULL,
     sample_rate        INTEGER     NOT NULL,
     Q0                 DOUBLE,
     Q                  DOUBLE,
     natural_frequency  DOUBLE,
     damping_factor     DOUBLE,
     meter_type         VARCHAR(32),
     origin             VARCHAR(64),
     data               BLOB,
    FOREIGN KEY(event_id) REFERENCES event(id) ON DELETE CASCADE ON UPDATE CASCADE,
    FOREIGN KEY(station_id) REFERENCES station(id) ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT uc_stev UNIQUE (event_id, station_id));

CREATE INDEX idx_seismogram_station ON seismogram (station_id ASC);
CREATE INDEX idx_seismogram_event ON seismogram (station_id ASC);
CREATE INDEX idx_seismogram_date ON seismogram (date ASC);


CREATE TABLE pick
    (id                 INTEGER     PRIMARY KEY,
     seismogram_id      INTEGER     NOT NULL,
     type               VARCHAR(16),
     date               DATE        NOT NULL,
     frame              INTEGER     NOT NULL,
    FOREIGN KEY(seismogram_id) REFERENCES seismogram(id) ON DELETE CASCADE ON UPDATE CASCADE);




CREATE TABLE catalogs
    (name               VARCHAR(64) UNIQUE);


CREATE TABLE original_event_catalog
    (event_id        INTEGER    PRIMARY KEY,
     date            DATE,
     X               DOUBLE      NOT NULL,
     Y               DOUBLE      NOT NULL,
     Z               DOUBLE,
    FOREIGN KEY(event_id) REFERENCES event(id) ON DELETE CASCADE ON UPDATE CASCADE);


CREATE TABLE original_traveltime_catalog
    (pick_id        INTEGER     PRIMARY KEY,
     traveltime     DOUBLE       NOT NULL,
     residual       DOUBLE,
    FOREIGN KEY(pick_id) REFERENCES pick(id) ON DELETE CASCADE ON UPDATE CASCADE);
"""


__TRAVELTIME_QUERY__ = string.Template(
"""
SELECT pick.id, station.id, event.id, ${catalog}_traveltime_catalog.traveltime,
       pick.type ${ending} ORDER BY ${order}
""")

__QUERY_ENDING__ = string.Template(
"""
        FROM event, station, seismogram, pick, ${catalog}_traveltime_catalog,
        ${catalog}_event_catalog
        WHERE seismogram.station_id                    = station.id AND
           seismogram.event_id                      = event.id AND
           pick.seismogram_id                       = seismogram.id AND
           ${catalog}_traveltime_catalog.pick_id    = pick.id AND
           event.id                                 = seismogram.event_id AND
           ${catalog}_event_catalog.event_id        = event.id
            ${filters}
""")



__EVENT_QUERY__ = string.Template(
"""
SELECT id, Z, Y, X FROM event WHERE id IN
    (SELECT event.id ${ending}) ORDER BY id;
""")

__STATION_QUERY__ = string.Template(
"""
SELECT id, Z, Y, X FROM station WHERE id IN
    (SELECT station.id ${ending}) ORDER BY id;
""")

__ORDER_BY_QUERY__ = string.Template(
"""
ORDER BY pick.type, station.id
""")

ev_dtype = [('id',          'int'),
            ('X',           'float'),
            ('Y',           'float'),
            ('Z',           'float')]

st_dtype = ev_dtype

tt_dtype = [('id',                'int'),
            ('station_id',          'int'),
            ('event_id',        'int'),
            ('traveltime',        'float'),
            ('type',              '|S16')]

def format_date(date):
    if isinstance(date, float):
        date = datetime.datetime.fromtimestamp(date).isoformat()
    if not isinstance(date, str):
        date = date.isoformat()
    return date

class ModelQueryBuilder(object):
    def __init__(self):
        self.filters = {}
        self.ptype_filter = None

    def set_type_filter(self, ptype):
        if ptype is None:
            self.ptype_filter = None
        else:
            self.ptype_filter = string.Template('AND pick.type = "%s"' % ptype)

    def set_event_filter(self, ids):
        if len(ids) == 0:
            return
        self.filters['events'] = string.Template(("(NOT (event.id IN (%s)))" % ",".join([str(i) for i in ids])))

    def set_station_filter(self, ids):
        if len(ids) == 0:
            return
        self.filters['stations'] = string.Template(("(NOT (station.id  IN (%s)))" % ",".join([str(i) for i in ids])))

    def set_inverse_station_filter(self, ids):
        if len(ids) == 0:
            return
        self.filters['stations'] = string.Template(("(station.id  IN (%s))" % ",".join([str(i) for i in ids])))

    def set_inverse_event_filter(self, ids):
        if len(ids) == 0:
            return
        self.filters['events'] = string.Template(("(event.id IN (%s))" % ",".join([str(i) for i in ids])))


    def set_cuboid_filter(self, origin, length, padding = 0):
        if (origin is None) or (length is None):
            del self.filters['cuboid']
        else:
            stemplate ='\n(${catalog}_event_catalog.%s BETWEEN %f AND %f)'
            limit = [o + l for o, l in zip(origin, length)]
            tlist = [stemplate % (c, origin[i], limit[i]) \
                                        for i, c in enumerate(['X', 'Y', 'Z'])]

            sstemplate = "\n(station.%s BETWEEN %f AND %f)"
            tslist = [sstemplate % (c, origin[i], limit[i]) \
                                        for i, c in enumerate(['X', 'Y', 'Z'])]



            self.filters['cuboid'] = string.Template("AND ".join(tlist + tslist))


    def set_date_filter(self, fromdate, todate):
        fromdate = format_date(fromdate)
        todate = format_date(todate)

        if (fromdate is None) or (todate is None):
            del self.filters['date']
        else:
            stemplate = string.Template('(${catalog}_event_catalog.date ' \
                                        'BETWEEN "%s" AND "%s")' % (fromdate, todate))
            self.filters['date'] = stemplate

    def __get_query_ending__(self, catalog, include_ptype = False):
        if len(self.filters) == 0:
            filters = ""
        else:
            lstfilt = [f.substitute(catalog = catalog) for f in self.filters.values()]
            filters = "AND " + " AND ".join(lstfilt)
        if include_ptype and (self.ptype_filter is not None):
            filters += self.ptype_filter.substitute(catalog = catalog)
        ending = __QUERY_ENDING__.substitute(catalog = catalog,
                                             filters = filters)
        return ending

    def traveltime_query(self, catalog = "original", order = "station.id, event.id"):
        ending = self.__get_query_ending__(catalog, include_ptype = True)
        return __TRAVELTIME_QUERY__.substitute(catalog = catalog,
                                               ending = ending,
                                               order = order)

    def event_query(self, catalog = "original"):
        ending = self.__get_query_ending__(catalog)
        return __EVENT_QUERY__.substitute(ending = ending, catalog = catalog)

    def station_query(self, catalog = "original"):
        ending = self.__get_query_ending__(catalog)
        return __STATION_QUERY__.substitute(ending = ending)

    def set_filters(self, filters):
        for f in filters:
            getattr(self, "set_%s_filter" % f['name'] )(*f['args'])


