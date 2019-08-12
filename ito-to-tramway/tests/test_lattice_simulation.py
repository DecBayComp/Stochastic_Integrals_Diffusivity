"""
Test simulation of RW in a space confined by beads in lattice nodes.


isort:skip_file
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simLattice.simulate import get_boundary_intersection, get_bead_intersection, get_reflections

atol = 1e-8


def test_boundary_crossing():
    L = 1
    r = np.array([0.2, 0.2])

    # Horizontal
    dr = np.array([-0.3, 0])
    leaving, old_intersection, intersection, new_dr = get_boundary_intersection(r, dr, L)
    assert leaving
    assert np.allclose(intersection, np.array([L, 0.2]), atol=atol)
    assert np.allclose(new_dr, np.array([-0.1, 0]), atol=atol)

    # Vertical
    dr = np.array([0, 1.5])
    leaving, old_intersection, intersection, new_dr = get_boundary_intersection(r, dr, L)
    assert np.allclose(intersection, np.array([r[0], 0]), atol=atol)
    assert np.allclose(new_dr, np.array([0, 0.7]), atol=atol)
    assert leaving

# print('Hello')


def test_bead_intersection():
    def ar(s): return np.array(s)

    P = ar([0, 0])
    R = 0.1

    # Horizontal intersection
    r = ar([0.2, 0])
    dr = ar([-0.15, atol])

    reachable, intersection, distance = get_bead_intersection(r, dr, P, R)
    assert reachable
    assert np.allclose(intersection, ar([0.1, 0]), atol=atol)
    assert np.isclose(distance, 0.1)

    # Vertical intersection
    r = ar([0, -0.2])
    dr = ar([atol, 0.17])

    reachable, intersection, distance = get_bead_intersection(r, dr, P, R)
    assert reachable
    assert np.allclose(intersection, ar([0, -0.1]), atol=atol)
    assert np.isclose(distance, 0.1)

    # Vertical intersection 2
    P = ar([0.5, 0.5])
    R = 0.1
    r = np.array([0.5, 0.3])
    dr = np.array([0, 0.4])
    reachable, intersection, distance = get_bead_intersection(r, dr, P, R)
    assert reachable
    assert np.allclose(intersection, ar([0.5, 0.4]), atol=atol)
    assert np.isclose(distance, 0.1)


def test_bead_reflection():
    Nx, Ny = 1, 1
    L = 1
    beads = np.array([[0, 0.3, 0.1]])

    def R_func(x): return 0.1

    # Horizontal
    r = np.array([0.2, 0.3])
    dr = np.array([-0.15, 0])
    reflected, new_r, dr_after_intersection, internal_trajectory, new_dr = get_reflections(
        r, dr, beads, R_func, Nx, Ny, L)
    assert reflected
    assert np.allclose(new_r, np.array([0.15, 0.3]), atol=atol)
    assert np.allclose(internal_trajectory, [[0.2, 0.3, -0.1, 0], [0.1, 0.3, 0.05, 0]], atol=atol)
    assert np.allclose(dr_after_intersection, [0.05, 0], atol=atol)

    # Vertical
    beads = np.array([[0.5, 0.5, 0.1]])
    r = np.array([0.5, 0.3])
    dr = np.array([0, 0.4])
    # reachable, intersection, distance = get_bead_intersection(r, dr, beads[0, :2], beads[0, 2])
    # assert reachable

    reflected, new_r, dr_after_intersection, internal_trajectory, new_dr = get_reflections(
        r, dr, beads, R_func, Nx, Ny, L)
    assert reflected
    assert np.allclose(new_r, np.array([0.5, 0.1]), atol=atol)
    assert np.allclose(new_dr, np.array([0, -0.2]), atol=atol)
    assert np.allclose(internal_trajectory, [[0.5, 0.3, 0, 0.1], [0.5, 0.4, 0, -0.3]], atol=atol)

    # 45 degrees NE
    beads = np.array([[0.5, 0.5, 0.1]])
    r = np.array([0.3, 0.3])
    dr = np.array([1, 1]) * 0.5
    reflected, new_r, dr_after_intersection, internal_trajectory, new_dr = get_reflections(
        r, dr, beads, R_func, Nx, Ny, L)

    assert reflected
    assert np.allclose(internal_trajectory, [
                       [0.3, 0.3, 0.2 - 0.1 / np.sqrt(2), 0.2 - 0.1 / np.sqrt(2)], [0.5 - 0.1 / np.sqrt(2), 0.5 - 0.1 / np.sqrt(2), -0.3 - 0.1 / np.sqrt(2), -0.3 - 0.1 / np.sqrt(2)]], atol=atol)

    # 45 degrees NE, non-zero incidence
    beads = np.array([[0.5, 0.5, 0.1]])
    r = np.array([0.4, 0.3])
    dr = np.array([0.3, 0.3])
    reflected, new_r, dr_after_intersection, internal_trajectory, new_dr = get_reflections(
        r, dr, beads, R_func, Nx, Ny, L)

    # print('int2', internal_trajectory)

    assert reflected
    assert np.allclose(internal_trajectory, [
                       [0.4, 0.3, 0.1, 0.1], [0.5, 0.4, 0.2, -0.2]], atol=atol)


def test_correct_full_dr_with_border_crossing():
    """Test that after all reflections, the actual displacement dr is correctly calculated"""
    Nx, Ny = 1, 1
    L = 1
    beads = np.array([[0, 0.3, 0.1]])

    def R_func(x): return 0.1

    # Horizontal
    r = np.array([0.1, 0.1])
    dr = np.array([-0.2, 0])
    reflected, new_r, dr_after_intersection, internal_trajectory, new_dr = get_reflections(
        r, dr, beads, R_func, Nx, Ny, L)
    # print('int', internal_trajectory)
    assert np.allclose(new_dr, np.array([-0.2, 0]), atol=atol)
    assert reflected
    # assert np.allclose(new_r, np.array([0.15, 0.3]), atol=atol)
    # assert np.allclose(new_dr, new_r - r, atol=atol)
    # assert np.allclose(internal_trajectory, [[0.2, 0.3, -0.1, 0], [0.1, 0.3, 0.05, 0]], atol=atol)
