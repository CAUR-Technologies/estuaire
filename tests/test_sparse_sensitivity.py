import numpy as np
from scipy import sparse

from estuaire.data import EKEventTable, EKImageData, EKStationTable, EKTTTable
from estuaire.frechet import compute_sparse_sensitivity


def _velocity_grid(shape=(11, 11), velocity=3.0):
    data = np.full(shape, velocity, dtype=float)
    return EKImageData(data, origin=(0.0, 0.0), spacing=1.0)


def _tables():
    station_dtype = EKStationTable.dtype
    event_dtype = EKEventTable.dtype
    tt_dtype = EKTTTable.dtype

    stations = EKStationTable(
        np.array([(0, (5.0, 5.0, 0.0), 0.0)], dtype=station_dtype)
    )

    events = EKEventTable(
        np.array(
            [
                (0, (8.0, 5.0, 0.0), 0.0),
                (1, (5.0, 8.0, 0.0), 0.0),
            ],
            dtype=event_dtype,
        )
    )

    tt_entries = np.array(
        [
            (0, 0, 1.0),
            (0, 1, 1.0),
        ],
        dtype=tt_dtype,
    )
    tt_table = EKTTTable(tt_entries, staid=0)
    tt_table.event_table = events
    tt_table.station_table = stations

    distances = np.linalg.norm(
        events.data['position'][:, :2] - stations.data['position'][0][:2],
        axis=1,
    )
    tt_table.data['traveltime'] = distances / 3.0

    return stations, events, tt_table


def test_compute_sparse_sensitivity_shapes():
    velocity_grid = _velocity_grid()
    stations, events, tt_table = _tables()

    result = compute_sparse_sensitivity(
        arrival=velocity_grid.copy(),
        velocity=velocity_grid,
        traveltime=tt_table,
        grid_id="grid",
        rk_step=1.0,
    )

    assert set(result) == {
        'residual',
        'grid',
        'station_time',
        'event_time',
        'event_position',
        'description',
    }

    residual = result['residual']
    grid = result['grid']
    station_matrix = result['station_time']
    event_matrix = result['event_time']
    event_position = result['event_position']
    description = result['description']

    assert residual.shape == (len(tt_table.data),)
    assert sparse.isspmatrix_csr(grid)
    assert grid.shape == (len(residual), velocity_grid.data.size)
    assert grid.nnz > 0

    assert station_matrix.shape == (len(residual), stations.data.size)
    assert event_matrix.shape == (len(residual), events.data.size)
    assert event_position.shape[0] == len(residual)

    assert description['grid'] == velocity_grid.data.shape


def test_compute_sparse_sensitivity_zero_events():
    velocity_grid = _velocity_grid()
    stations, _, tt_table = _tables()
    tt_table.data = tt_table.data[:0]

    result = compute_sparse_sensitivity(
        arrival=velocity_grid.copy(),
        velocity=velocity_grid,
        traveltime=tt_table,
        grid_id="grid",
    )

    assert result['residual'] is None
    assert result['grid'] is None
    assert result['station_time'] is None
    assert result['event_time'] is None
    assert result['event_position'] is None
