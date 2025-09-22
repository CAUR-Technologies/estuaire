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

from . import raytrace, solver
from .data import EKImageData

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
    dtype: Union[str, np.dtype] = np.float64,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Compute Frechet derivatives between sources and receivers.

    Parameters
    ----------
    velocity
        Velocity grid as an :class:`~eikonal.data.EKImageData` or a raw ``ndarray``
        in ``C`` order. When an array is supplied, ``spacing`` (and optionally
        ``origin``) must also be provided.
    sources, receivers
        Iterable of points expressed in the same coordinate frame as the grid
        origin/spacing. Shapes are ``(n_sources, ndim)`` and ``(n_receivers, ndim)``.
    spacing
        Grid spacing (in each dimension) when ``velocity`` is given as an array.
        The low-level solver currently requires uniform spacing; anisotropic values
        raise ``ValueError``.
    origin
        Grid origin when ``velocity`` is a raw array. Defaults to zeros.
    rk_step
        Integration step passed to the Runge-Kutta ray tracer.
    second_order
        Toggle the second-order fast marching stencil.
    cell_slowness
        If ``True`` derivatives are reported with respect to cell slowness; when
        ``False`` they are with respect to velocity.
    return_travel_times
        If ``True`` also return the travel times for each source/receiver pair.
    pairwise
        When ``True`` treat ``sources`` and ``receivers`` as aligned pairs and
        return a ``(n_pairs, n_cells)`` matrix.
    dtype
        Output dtype for the Frechet matrix (travel times are always float64).

    Returns
    -------
    frechet : ndarray
        ``(n_sources, n_receivers, n_cells)`` Frechet array (or ``(n_pairs, n_cells)``
        in pairwise mode).
    travel_times : ndarray, optional
        Returned only when ``return_travel_times`` is ``True``.
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
        raise ValueError("pairwise=True requires the same number of sources and receivers.")

    n_sources = src_grid.shape[0]
    n_receivers = rcv_grid.shape[0]
    n_cells = int(vel_data.size)

    frechet_rows = (
        np.zeros((n_sources, n_cells), dtype=np.float64)
        if pairwise
        else np.zeros((n_sources, n_receivers, n_cells), dtype=np.float64)
    )
    travel = (
        np.zeros(n_sources, dtype=np.float64)
        if pairwise and return_travel_times
        else np.zeros((n_sources, n_receivers), dtype=np.float64)
        if return_travel_times
        else None
    )

    velocity_flat = vel_data.ravel(order="C")
    velocity_flat_sq = velocity_flat * velocity_flat

    max_len = np.linalg.norm(np.asarray(vel_data.shape, dtype=np.float64))
    max_len *= 4 ** dims
    rk_step = float(rk_step)
    if rk_step <= 0:
        raise ValueError("rk_step must be strictly positive.")
    buffer_len = max(1, int(np.ceil(max_len * (8.0 / rk_step))))
    indices_buffer = np.zeros(buffer_len, dtype=np.uintp)
    sens_buffer = np.zeros(buffer_len, dtype=np.float64)

    row_buffer = np.zeros(n_cells, dtype=np.float64)

    for idx, seed in enumerate(src_grid):
        arrival = _solve_arrival(vel_data, seed, spacing_value, second_order)
        arrival = np.require(arrival, dtype=np.float64, requirements=("C",))

        receiver_indices: Iterable[int]
        if pairwise:
            receiver_indices = (idx,)
        else:
            receiver_indices = range(n_receivers)

        for rcv_idx in receiver_indices:
            target = rcv_grid[rcv_idx]
            row_buffer.fill(0.0)

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

            if used_idx.size == 0:
                if pairwise:
                    frechet_rows[idx, :] = 0.0
                    if travel is not None:
                        travel[idx] = tt_val
                else:
                    frechet_rows[idx, rcv_idx, :] = 0.0
                    if travel is not None:
                        travel[idx, rcv_idx] = tt_val
                continue

            row_buffer[used_idx] = weights

            if not cell_slowness:
                mask = velocity_flat_sq > 0.0
                with np.errstate(divide="ignore", invalid="ignore"):
                    np.divide(
                        -row_buffer,
                        velocity_flat_sq,
                        out=row_buffer,
                        where=mask,
                    )
                row_buffer[~mask] = 0.0
                row_buffer[~np.isfinite(row_buffer)] = 0.0

            if pairwise:
                frechet_rows[idx, :] = row_buffer
                if travel is not None:
                    travel[idx] = tt_val
            else:
                frechet_rows[idx, rcv_idx, :] = row_buffer
                if travel is not None:
                    travel[idx, rcv_idx] = tt_val

    frechet_rows = frechet_rows.astype(dtype, copy=False)
    if travel is None:
        return frechet_rows
    return frechet_rows, travel


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
