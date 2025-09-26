#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
"""Frechet derivative helpers backed by the C++ fast-marching kernels."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sparse

from estuaire import raytrace, solver
from .data import EKImageData

from tqdm import tqdm

ArrayLike = Union[Sequence[float], np.ndarray]


def compute_frechet(
    velocity: Union[EKImageData, np.ndarray],
    sources: ArrayLike,
    receivers: ArrayLike,
    *,
    spacing: Optional[Union[float, Sequence[float]]] = None,
    origin: Optional[Sequence[float]] = None,
    rk_step: float = 1.0,
    second_order: bool = True,
    cell_slowness: bool = True,
    return_travel_times: bool = False,
    pairwise: bool = False,
    return_rays: bool = False,
    dtype: Union[str, np.dtype] = np.float64,
    return_dense: bool = False,
) -> Union[
    sparse.csr_matrix,
    Tuple[sparse.csr_matrix, np.ndarray],
    Tuple[np.ndarray, np.ndarray],
    np.ndarray,
]:
    """Compute Frechet derivatives and return a CSR sensitivity matrix.

    Parameters
    ----------
    velocity, sources, receivers
        Grid and coordinate inputs share the same interpretation as in the
        original fast-marching helper.
    return_dense
        When ``True`` materialise the sensitivity matrix as a dense NumPy array
        (mirroring the legacy behaviour). By default a CSR matrix is returned.

    Returns
    -------
    frechet : scipy.sparse.csr_matrix or ndarray
        Sparse matrix with one row per source/receiver combination. The optional
        ``original_shape`` attribute stores the ``(n_sources, n_receivers)``
        layout when ``pairwise`` is ``False``. When ``return_dense=True`` a
        dense array with the same ordering as the original helper is produced
        instead.
    travel_times : ndarray, optional
        Returned when ``return_travel_times`` is ``True``.
    """
    grid = _ensure_grid(velocity, spacing=spacing, origin=origin)

    vel_data = np.require(grid.data, dtype=np.float64, requirements=("C",))
    dims = vel_data.ndim

    src_points = _prepare_points(sources, dims, "sources")
    rcv_points = _prepare_points(receivers, dims, "receivers")

    origin_vec = _as_vector(grid.origin, dims, "origin")
    spacing_vec = _as_vector(grid.spacing, dims, "spacing")
    _assert_isotropic_spacing(spacing_vec)
    spacing_value = float(spacing_vec[0])
    if spacing_value <= 0:
        raise ValueError("Spacing must be strictly positive.")

    src_grid = (src_points - origin_vec) / spacing_vec
    rcv_grid = (rcv_points - origin_vec) / spacing_vec

    _validate_points(src_grid, vel_data.shape, "source")
    _validate_points(rcv_grid, vel_data.shape, "receiver")

    if pairwise and src_grid.shape[0] != rcv_grid.shape[0]:
        raise ValueError(
            "pairwise=True requires the same number of sources and receivers."
        )

    n_sources = src_grid.shape[0]
    n_receivers = rcv_grid.shape[0]
    n_cells = int(vel_data.size)
    n_rows = n_sources if pairwise else n_sources * n_receivers

    travel = (
        np.zeros(n_sources, dtype=np.float64)
        if pairwise and return_travel_times
        else np.zeros((n_sources, n_receivers), dtype=np.float64)
        if return_travel_times
        else None
    )

    if return_rays:
        raise NotImplementedError("return_rays is not supported in sparse mode.")

    velocity_flat = vel_data.ravel(order="C")
    velocity_flat_sq = velocity_flat * velocity_flat

    max_len = np.linalg.norm(np.asarray(vel_data.shape, dtype=np.float64)) * (
        4 ** dims
    )
    rk_step = float(rk_step)
    if rk_step <= 0:
        raise ValueError("rk_step must be strictly positive.")
    buffer_len = max(1, int(np.ceil(max_len * (8.0 / rk_step))))
    indices_buffer = np.zeros(buffer_len, dtype=np.uintp)
    sens_buffer = np.zeros(buffer_len, dtype=np.float64)

    data_entries: list[float] = []
    index_entries: list[int] = []
    indptr = [0]
    row_counter = 0

    for src_idx, seed in enumerate(src_grid):
        arrival = _solve_arrival(vel_data, seed, spacing_value, second_order)
        arrival = np.require(arrival, dtype=np.float64, requirements=("C",))

        receiver_indices: Iterable[int]
        if pairwise:
            receiver_indices = (src_idx,)
        else:
            receiver_indices = range(n_receivers)

        for rcv_idx in receiver_indices:
            target = rcv_grid[rcv_idx]

            tt_val, used_idx, weights = raytrace.sensivity(
                arrival,
                vel_data,
                tuple(seed.tolist()),
                tuple(target.tolist()),
                indices_buffer,
                sens_buffer,
                spacing_value,
                h=rk_step,
            )

            indices_vals = np.asarray(used_idx[: used_idx.size], dtype=np.intp)
            weights_vals = np.asarray(weights[: indices_vals.size], dtype=dtype)

            if indices_vals.size == 0:
                indptr.append(indptr[-1])
                row_counter += 1
                if travel is not None:
                    if pairwise:
                        travel[src_idx] = tt_val
                    else:
                        travel[src_idx, rcv_idx] = tt_val
                continue

            factors = velocity_flat_sq[indices_vals]
            mask = factors > 0.0
            if cell_slowness:
                weights_vals = np.where(mask, -weights_vals * factors, 0.0)
            else:
                weights_vals = np.where(mask, weights_vals, 0.0)

            weights_vals = np.asarray(weights_vals, dtype=dtype)
            weights_vals[~np.isfinite(weights_vals)] = 0.0
            nonzero_mask = weights_vals != 0.0
            weights_nonzero = weights_vals[nonzero_mask]
            indices_nonzero = indices_vals[nonzero_mask]

            index_entries.extend(indices_nonzero.tolist())
            data_entries.extend(weights_nonzero.tolist())
            indptr.append(indptr[-1] + weights_nonzero.size)
            row_counter += 1

            if travel is not None:
                if pairwise:
                    travel[src_idx] = tt_val
                else:
                    travel[src_idx, rcv_idx] = tt_val

    while row_counter < n_rows:
        indptr.append(indptr[-1])
        row_counter += 1

    indices_array = np.array(index_entries, dtype=np.int32)
    data_array = np.array(data_entries, dtype=dtype)
    indptr_array = np.array(indptr, dtype=np.int32)

    matrix = sparse.csr_matrix(
        (data_array, indices_array, indptr_array),
        shape=(n_rows, n_cells),
    )

    if pairwise:
        matrix.original_shape = (n_sources,)
        if return_dense:
            frechet_output = matrix.toarray()
        else:
            frechet_output = matrix
        if travel is not None:
            return frechet_output, travel
        return frechet_output

    # Non-pairwise: rows are in src-major order with rcv varying fastest.
    matrix.original_shape = (n_sources, n_receivers)

    if return_dense:
        frechet_output = matrix.toarray().reshape(n_sources, n_receivers, n_cells)
    else:
        frechet_output = matrix

    if travel is not None:
        return frechet_output, travel
    return frechet_output





def compute_sparse_sensitivity(
    arrival,
    velocity,
    traveltime,
    grid_id,
    *,
    rk_step: float = 1.0,
) -> dict:
    """Replicate the legacy sparse Frechet dictionary used by ``bin/sensitivity``."""

    if (arrival.spacing != velocity.spacing) or (arrival.shape != velocity.shape):
        raise AttributeError("Arrival and velocity grids must share spacing and shape")

    tttable = traveltime
    stdesc = tttable.station_row
    current_events = tttable.event_rows

    agrid = arrival.data
    vgrid = velocity.data
    velocity_flat_sq = vgrid.ravel(order="C") ** 2

    ndim = arrival.data.ndim
    start = tuple(arrival.transform_to(stdesc['position'][:ndim]))
    if np.any(np.sqrt(np.sum((arrival.seeds - start) ** 2, axis=1)) > 0.0001):
        raise ValueError("Station position does not match arrival grid seeds")

    max_len = np.linalg.norm(np.array(vgrid.shape, dtype=np.float64)) * (4 ** vgrid.ndim)
    rk_step = float(rk_step)
    buffer_len = max(1, int(np.ceil(max_len * (8.0 / rk_step))))
    indices_buffer = np.zeros(buffer_len, dtype=np.uintp)
    sensitivity_buffer = np.zeros(buffer_len, dtype=np.float64)

    grid_indices_x: list[int] = []
    grid_indices_y: list[int] = []
    grid_sens: list[float] = []
    residuals: list[float] = []

    epa_x: list[int] = []
    epa_y: list[int] = []
    epa_v: list[float] = []

    description = {
        grid_id: agrid.shape,
        "event_position": (tttable.event_table.size, ndim),
        "event_time": (tttable.event_table.size,),
        "station_time": (tttable.station_table.size,),
    }

    if len(current_events) == 0:
        return {
            "residual": None,
            grid_id: None,
            "station_time": None,
            "event_time": None,
            "event_position": None,
            "description": description,
        }

    for idx, event_row in tqdm(enumerate(current_events)):
        finish = tuple(arrival.transform_to(event_row['position'][:ndim]))
        rtt, ind, sens = raytrace.sensivity(
            agrid,
            vgrid,
            start,
            finish,
            indices_buffer,
            sensitivity_buffer,
            velocity.spacing,
            h=rk_step,
        )

        tt = (
            tttable.data['traveltime'][idx]
            + stdesc['delta_t']
            - event_row['delta_t']
        )
        residuals.append(tt - rtt)

        for axis_idx, coord in enumerate(finish):
            delta = 0.1
            new_finish = list(finish)
            new_finish[axis_idx] = coord + delta
            new_finish1 = tuple(new_finish)
            new_finish[axis_idx] = coord - delta
            new_finish2 = tuple(new_finish)

            tt1 = estuaire.bin.traveltime(
                agrid,
                vgrid,
                start,
                new_finish1,
                velocity.spacing,
                h=rk_step,
            )
            tt2 = estuaire.bin.traveltime(
                agrid,
                vgrid,
                start,
                new_finish2,
                velocity.spacing,
                h=rk_step,
            )
            epa_y.append(idx)
            epa_x.append(tttable.data['event_id'][idx] * ndim + axis_idx)
            epa_v.append((tt1 - tt2) / (velocity.spacing * delta))

        contributions = sens[: ind.size]
        if cell_slowness:
            contributions = -contributions * velocity_flat_sq[ind[: ind.size]]
        grid_indices_x.extend(ind[: ind.size].tolist())
        grid_sens.extend(contributions.tolist())
        grid_indices_y.extend([idx] * ind.size)

    residual = np.asarray(residuals, dtype=float)
    grid_matrix = sparse.coo_matrix(
        (grid_sens, (grid_indices_y, grid_indices_x)),
        shape=(residual.size, vgrid.size),
    ).tocsr()

    epa_matrix = sparse.coo_matrix(
        (epa_v, (epa_y, epa_x)),
        shape=(residual.size, tttable.event_table.size * ndim),
    ).tocsr()

    station_indices = np.full(residual.size, tttable.station_id, dtype=int)
    station_matrix = sparse.coo_matrix(
        (-np.ones_like(station_indices, dtype=float), (np.arange(residual.size), station_indices)),
        shape=(residual.size, tttable.station_table.size),
    ).tocsr()

    event_matrix = sparse.coo_matrix(
        (
            np.ones(residual.size, dtype=float),
            (np.arange(residual.size), tttable.data['event_id']),
        ),
        shape=(residual.size, tttable.event_table.size),
    ).tocsr()

    return {
        "residual": residual,
        grid_id: grid_matrix,
        "station_time": station_matrix,
        "event_time": event_matrix,
        "event_position": epa_matrix,
        "description": description,
    }


def _ensure_grid(
    velocity: Union[EKImageData, np.ndarray],
    *,
    spacing: Optional[Union[float, Sequence[float]]],
    origin: Optional[Sequence[float]],
) -> EKImageData:
    if isinstance(velocity, EKImageData):
        return velocity

    if spacing is None:
        raise ValueError("spacing must be provided when velocity is an ndarray.")

    velocity_array = np.asarray(velocity, dtype=np.float64)
    if velocity_array.ndim not in (2, 3):
        raise ValueError("The velocity grid must be 2-D or 3-D.")

    spacing_vec = _as_vector(spacing, velocity_array.ndim, "spacing")

    if origin is None:
        origin_vec = np.zeros(velocity_array.ndim, dtype=np.float64)
    else:
        origin_vec = _as_vector(origin, velocity_array.ndim, "origin")

    return EKImageData(velocity_array.copy(), origin=origin_vec, spacing=spacing_vec)


def _prepare_points(points: ArrayLike, dims: int, label: str) -> np.ndarray:
    arr = np.asarray(points, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.ndim != 2:
        raise ValueError(f"{label} must be a 1-D or 2-D array-like of coordinates.")
    if arr.shape[1] != dims:
        raise ValueError(f"Expected {dims}-component {label}, got shape {arr.shape}.")
    return arr


def _as_vector(value: Union[float, Sequence[float]], dims: int, name: str) -> np.ndarray:
    if np.isscalar(value):
        vec = np.full(dims, float(value), dtype=np.float64)
    else:
        vec = np.asarray(value, dtype=np.float64)
        if vec.size != dims:
            raise ValueError(f"{name} must have exactly {dims} elements.")
        vec = vec.reshape(dims)
    return vec


def _assert_isotropic_spacing(spacing: np.ndarray) -> None:
    if not np.allclose(spacing, spacing[0]):
        raise ValueError("Anisotropic spacing is not supported by the C++ backend.")


def _validate_points(points: np.ndarray, shape: Tuple[int, ...], role: str) -> None:
    shape_arr = np.asarray(shape, dtype=np.float64)
    lower = 0.0
    upper = shape_arr - 1.0
    if np.any(points < lower) or np.any(points > upper):
        raise ValueError(f"A {role} lies outside the grid bounds.")


def _solve_arrival(
    velocity: np.ndarray,
    seed: np.ndarray,
    spacing: float,
    second_order: bool,
) -> np.ndarray:
    dims = velocity.ndim
    inner = tuple(slice(2, -2) for _ in range(dims))

    padded_shape = tuple(dim + 4 for dim in velocity.shape)
    viscosity = np.zeros(padded_shape, dtype=np.float64)
    viscosity[inner] = 1.0 / velocity

    tag = np.full(padded_shape, 2, dtype=np.int32)
    tag[inner] = 0

    arrival = np.empty(padded_shape, dtype=np.float64)
    arrival.fill(np.inf)

    seed_array = np.ascontiguousarray(seed.reshape(1, dims), dtype=np.float64)
    solver.SFMM(
        seed_array + 2.0,
        np.zeros(seed_array.shape[0], dtype=np.float64),
        tag,
        viscosity,
        arrival,
        spacing,
        second_order=second_order,
    )

    return arrival[inner]
