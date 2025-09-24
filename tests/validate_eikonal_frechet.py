"""Validation harness for :mod:`estuaire.frechet`.

The script mirrors the high-level checks found in ``uquake/tests`` but exercises
the new :func:`estuaire.frechet.compute_frechet` helper directly.  It verifies that

* path-length sensitivities (slowness parameterisation) sum to the Euclidean
  source/receiver distance;
* travel times follow ``distance / velocity``;
* switching to the velocity parameterisation simply applies the
  ``-L / v**2`` scaling;
* pairwise mode matches the dense ``(N_src, N_rcv, N_cells)`` output.

Run the script from ``libraries/estuary/src/estuaire`` after building the C++
extensions, e.g.::

    python setup.py build_ext --inplace
    python tests/validate_estuaire_frechet.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from estuaire.frechet import compute_frechet
except ImportError as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "Unable to import 'estuaire.frechet'. Build the Estuary extensions first by "
        "running 'python setup.py build_ext --inplace' from the project root."
    ) from exc


VELOCITY = 3200.0  # m/s
SPACING = 1.0      # m
ORIGIN = (0.0, 0.0)
GRID_SHAPE = (40, 40)
TOL = 1e-3

def _build_test_grid() -> np.ndarray:
    """Return a constant-velocity grid."""

    return np.full(GRID_SHAPE, VELOCITY, dtype=np.float64)


def _build_points(n_points: int, rng: np.random.Generator) -> np.ndarray:
    """Create random points within the interior of the grid (avoid boundaries)."""

    lower = 5.0
    upper = GRID_SHAPE[0] - 5.0
    coords = rng.uniform(low=lower, high=upper, size=(n_points, 2))
    return coords


def _check_slowness_pairwise(slowness_rows: np.ndarray,
                             travel_times: np.ndarray,
                             sources: np.ndarray,
                             receivers: np.ndarray) -> None:
    """Assert the slowness sensitivities integrate to geometric path lengths."""

    for idx, (row, src, rcv, tt) in enumerate(zip(slowness_rows, sources, receivers, travel_times)):
        distance = math.dist(src, rcv) * SPACING
        length_sum = row.sum()
        expected_tt = distance / VELOCITY

        if not math.isclose(length_sum, distance, rel_tol=TOL, abs_tol=TOL):
            raise AssertionError(
                f"Path length mismatch for pair {idx}: expected {distance:.6f} m, got {length_sum:.6f} m"
            )
        if not math.isclose(tt, expected_tt, rel_tol=TOL, abs_tol=TOL):
            raise AssertionError(
                f"Travel-time mismatch for pair {idx}: expected {expected_tt:.6f} s, got {tt:.6f} s"
            )


def _check_velocity_parameterisation(slowness_rows: np.ndarray,
                                     velocity_rows: np.ndarray) -> None:
    """Verify the velocity derivatives follow the -L / v^2 scaling point-wise."""

    expected = -slowness_rows / (VELOCITY ** 2)
    if not np.allclose(velocity_rows, expected, rtol=1e-4, atol=1e-8):
        diff = np.abs(velocity_rows - expected).max()
        raise AssertionError(f"Velocity-parameterised Frechet mismatch (max |Δ| = {diff:.3e}).")


def _check_pairwise_vs_full(slowness_rows: np.ndarray,
                            slowness_full: np.ndarray,
                            travel_times: np.ndarray,
                            travel_full: np.ndarray) -> None:
    """Ensure pairwise output collapses the dense tensor along the diagonal."""

    diag_rows = np.stack([
        slowness_full[i, i, :]
        for i in range(slowness_rows.shape[0])
    ], axis=0)

    if not np.allclose(diag_rows, slowness_rows, rtol=1e-4, atol=1e-8):
        raise AssertionError("Pairwise sensitivities do not match dense output along the diagonal.")

    if travel_full is not None:
        diag_tt = np.array([travel_full[i, i] for i in range(travel_full.shape[0])])
        if not np.allclose(diag_tt, travel_times, rtol=1e-4, atol=1e-8):
            raise AssertionError("Pairwise travel times do not match diagonal entries of dense tensor.")


def main() -> None:
    velocity = _build_test_grid()
    rng = np.random.default_rng(2025)

    base_sources = _build_points(8, rng)
    base_receivers = _build_points(8, rng)

    # Inject duplicates to exercise caching paths inside compute_frechet
    sources = np.vstack([base_sources, base_sources[:2]])
    receivers = np.vstack([base_receivers, base_receivers[:2]])

    velocity = velocity + np.random.randn(velocity.shape[0], velocity.shape[1]) * 500

    # Pairwise slowness (returns path lengths) with travel times
slowness_matrix, travel_times = compute_frechet(
    velocity,
    sources,
    receivers,
    spacing=SPACING,
    origin=ORIGIN,
    pairwise=True,
    cell_slowness=True,
    return_travel_times=True,
    rk_step=1.0,
    second_order=True,
)
slowness_rows = slowness_matrix.toarray()

_check_slowness_pairwise(slowness_rows, travel_times, sources, receivers)

    # Dense tensor (no pairwise collapse)
slowness_full_matrix, travel_full = compute_frechet(
    velocity,
    sources,
    receivers,
    spacing=SPACING,
    origin=ORIGIN,
    pairwise=False,
    cell_slowness=True,
    return_travel_times=True,
    rk_step=1.0,
    second_order=True,
)

orig_shape = getattr(slowness_full_matrix, 'original_shape', None)
slowness_full = slowness_full_matrix.toarray()
if orig_shape is not None:
    slowness_full = slowness_full.reshape(*orig_shape, slowness_full_matrix.shape[1])

_check_pairwise_vs_full(slowness_rows, slowness_full, travel_times, travel_full)

# Velocity parameterisation should simply scale rows element-wise
velocity_matrix = compute_frechet(
    velocity,
    sources,
    receivers,
    spacing=SPACING,
    origin=ORIGIN,
    pairwise=True,
    cell_slowness=False,
    return_travel_times=False,
    rk_step=1.0,
    second_order=True,
)

velocity_rows = velocity_matrix.toarray()

_check_velocity_parameterisation(slowness_rows, velocity_rows)

print(
    "Validation successful: estuaire.frechet produces consistent sensitivities "
    f"for {sources.shape[0]} sources/receivers (grid size {velocity.size})."
)


if __name__ == "__main__":
    main()
