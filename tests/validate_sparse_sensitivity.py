"""Interactive validation helper for `compute_sparse_sensitivity`.

Run once the Estuary extension is built::

    python setup.py build_ext --inplace
    python tests/validate_sparse_sensitivity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from estuaire import solver

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from estuaire.data import EKEventTable, EKImageData, EKStationTable, EKTTTable
from estuaire.frechet import compute_frechet, compute_sparse_sensitivity


VELOCITY = 3_000.0
SPACING = 1.0
GRID_SHAPE = (60, 60)


class MemoryTTTable(EKTTTable):
    """In-memory variant wiring event and station tables directly."""

    def __init__(self, data, staid, events, stations):
        super().__init__(data, staid)
        self._events = events
        self._stations = stations

    @property
    def event_table(self):
        return self._events

    @property
    def station_table(self):
        return self._stations


def _build_velocity_grid() -> EKImageData:
    data = np.full(GRID_SHAPE, VELOCITY, dtype=float)
    return EKImageData(data, origin=(0.0, 0.0), spacing=SPACING)


def _build_tables() -> tuple[EKStationTable, EKEventTable, MemoryTTTable]:
    stations = EKStationTable(
        np.array([(0, (30.0, 30.0, 0.0), 0.0)], dtype=EKStationTable.dtype)
    )

    events = EKEventTable(
        np.array(
            [
                (0, (45.0, 30.0, 0.0), 0.0),
                (1, (30.0, 45.0, 0.0), 0.0),
                (2, (20.0, 40.0, 0.0), 0.0),
            ],
            dtype=EKEventTable.dtype,
        )
    )

    traveltimes = []
    station_pos = stations.data['position'][0][:2]
    station_id = stations.data['id'][0]
    for event in events.data:
        distance = np.linalg.norm(event['position'][:2] - station_pos)
        traveltimes.append((station_id, event['id'], distance / VELOCITY))

    tt_table = MemoryTTTable(
        np.array(traveltimes, dtype=EKTTTable.dtype),
        staid=station_id,
        events=events,
        stations=stations,
    )
    return stations, events, tt_table


def _arrival_grid_for_station(velocity_grid: EKImageData, station_pos: np.ndarray) -> EKImageData:
    dims = velocity_grid.data.ndim
    spacing_value = float(SPACING)

    inner = tuple(slice(2, -2) for _ in range(dims))
    padded_shape = tuple(dim + 4 for dim in velocity_grid.data.shape)

    viscosity = np.zeros(padded_shape, dtype=np.float64)
    viscosity[inner] = 1.0 / velocity_grid.data

    tag = np.full(padded_shape, 2, dtype=np.int32)
    tag[inner] = 0

    arrival = np.empty(padded_shape, dtype=np.float64)
    arrival.fill(np.inf)

    seed_physical = np.asarray(station_pos[:dims], dtype=np.float64)
    seed_grid = (seed_physical - np.asarray(velocity_grid.origin[:dims], dtype=np.float64)) / spacing_value
    seeds = np.ascontiguousarray(seed_grid.reshape(1, dims) + 2.0)
    start_times = np.zeros(seeds.shape[0], dtype=np.float64)

    solver.SFMM(seeds, start_times, tag, viscosity, arrival, spacing_value, second_order=True)

    arrival_data = arrival[inner]
    arrival_grid = EKImageData(arrival_data, origin=velocity_grid.origin, spacing=velocity_grid.spacing)
    arrival_grid.seeds = np.ascontiguousarray(seed_grid.reshape(1, dims))
    return arrival_grid


def _summarise_sparse(result: dict) -> None:
    residual = result['residual']
    grid = result['grid']

    print("Residual statistics:")
    print(f"  count: {residual.size}")
    print(f"  mean : {residual.mean():.6f} s")
    print(f"  std  : {residual.std():.6f} s")

    print("\nGrid sensitivity matrix:")
    print(f"  shape: {grid.shape}")
    print(f"  nnz  : {grid.nnz}")
    nnz_per_row = grid.getnnz(axis=1)
    print(
        "  nnz/row min={:d} max={:d} avg={:.1f}".format(
            nnz_per_row.min(), nnz_per_row.max(), nnz_per_row.mean()
        )
    )

    sums = np.asarray(grid.sum(axis=1)).ravel()
    print("\nPath-length sums (metres):")
    for idx, total in enumerate(sums):
        print(f"  row {idx}: {total:.6f}")


def _compare_dense_sparse(result: dict, frechet_dense: np.ndarray) -> None:
    sparse_rows = result['grid'].toarray().reshape(frechet_dense.shape[0], -1)
    dense_rows = frechet_dense.reshape(frechet_dense.shape[0], -1)

    print("\nMax absolute difference between dense and sparse rows:")
    for idx, (dense_row, sparse_row) in enumerate(zip(dense_rows, sparse_rows)):
        diff = np.max(np.abs(dense_row - sparse_row))
        print(f"  row {idx}: {diff:.6e}")


def main() -> None:
    velocity_grid = _build_velocity_grid()
    stations, events, tt_table = _build_tables()
    arrival_grid = _arrival_grid_for_station(velocity_grid, stations.data['position'][0])

    sparse_result = compute_sparse_sensitivity(
        arrival=arrival_grid,
        velocity=velocity_grid,
        traveltime=tt_table,
        grid_id="grid",
        rk_step=1.0,
    )

    _summarise_sparse(sparse_result)

    frechet_matrix, tt_dense = compute_frechet(
        velocity=velocity_grid,
        sources=stations.data['position'][:, :2],
        receivers=events.data['position'][:, :2],
        spacing=SPACING,
        origin=(0.0, 0.0),
        pairwise=True,
        return_travel_times=True,
    )

    frechet_rows = frechet_matrix.toarray()

    print("\nFrechet CSR matrix shape:", frechet_matrix.shape)
    print("Travel-time matrix shape:", tt_dense.shape)

    _compare_dense_sparse(sparse_result, frechet_rows)


if __name__ == "__main__":
    main()
