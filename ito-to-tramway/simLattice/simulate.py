"""
This file contains code simulating diffusion of a particle in a 2D region containing impenetrable circles, on the border of which, the particle experiences a fully elastic reflection.
The radius of circles increases with x thus creating a diffusivity gradient on large scales.
Periodic boundary conditions.
"""

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import norm
from tqdm import trange

# try:
#     bl_has_run
# except Exception:

#     %matplotlib
#     %load_ext autoreload
#     %autoreload 2
#     bl_has_run = True


# matplotlib.use('Agg')
atol = 1e-8


# %%
def simulate():
    """Main simulation routine"""

    # Input parameters
    N = np.long(2e3)
    internal_steps_number = 10
    dt = 0.04   # s
    dim = 2
    D = 0.01  # in um^2/s
    L = 1  # um
    beads_dx = beads_dy = L / 10
    Rmin = beads_dx / 5
    Rmax = Rmin * 2

    # Derived parameters
    dt_int = dt / internal_steps_number
    rs = np.full((N + 1, dim), np.nan, dtype=np.float32)
    drs = np.full((N + 1, dim), np.nan, dtype=np.float32)
    ts = np.arange(N + 1) * dt

    # Static beads array
    def R_func(x):
        return Rmax + (Rmin - Rmax) * np.abs(2 * x / L - 1)

    Nx = np.floor(L / beads_dx).astype(int) + 1
    Ny = np.floor(L / beads_dy).astype(int) + 1
    beads = np.zeros((Nx * Ny, 3))   # Columns: x, y, R
    for ix, iy in np.ndindex(Nx, Ny):
        x, y = ix * beads_dx, iy * beads_dy
        beads[ix * Ny + iy, :] = [x, y, R_func(x)]

    # Initialize
    # np.random.seed()

    # r0 = np.random.rand(dim) * L
    # while (not check_outside_beads(r0, beads)):
    #     r0 = np.random.rand(dim) * L

    # Start randomly outside beads
    count = 0
    while(count < 1e3):
        r0 = np.random.rand(dim) * L
        count += 1
        if check_outside_beads(r0, beads):
            print(f'Found starting point in {count} iterations')
            break

    rs[0, :] = r = r0
    noise = np.random.normal(0.0, 1.0, size=(dim, N, internal_steps_number)) * np.sqrt(dt_int)

    all_internal_trajectory = []
    for big_step in trange(N):
        # dr_big = [0, 0]

        for small_step in range(internal_steps_number):

            alpha = np.array([0, 0])  # Ito drift term
            b = 2.0 * D     # diffusivity term
            dW = noise[:, big_step, small_step]

            dr = alpha * dt_int + b * dW

            _, _, dr, internal_trajectory = get_reflections(r, dr, beads, R_func, Nx, Ny)
            all_internal_trajectory += internal_trajectory
            # dr = adjust_periodic(r, dr)

            r = r + dr
        rs[big_step + 1, :] = r
        # dr[:, big_step] = dr

    # print('Trajectory: ', rs)
    return ts, rs, beads, all_internal_trajectory


def get_reflections(r, dr, beads, R_func, Nx, Ny, L=1):
    """Calculate reflections from static beads.
    The first bead is located at (0,0), all others are periodic.
    The radius increases linearly with L

    Returns:
    internal_trajectory - N x 4, array-like, the subdivisions of the initial displacement into sub-intervals of the form [[x, y, dx, dy],...]

    new_dr = new_r - r, the new jump length from the original position
    """

    # # Temporary
    # L = 1  # um
    # beads_dx = beads_dy = L / 10
    # Rmin = beads_dx / 5
    # Rmax = Rmin * 2
    # # R = beads_rmin

    # Initialize
    r = np.array(r)
    dr = np.array(dr)
    reflected = False
    internal_trajectory = []

    def P(i):
        return beads[i, :2]

    # Cycle through beads locations
    distances = np.full(Nx * Ny, np.inf)
    for i in range(Nx * Ny):
        reachable, _, distance = get_bead_intersection(r, dr, P(i), R_func(P(i)[0]))
        if reachable:
            distances[i] = distance

    # Start with the closest one
    i = np.argmin(distances)
    # print('Dists', distances, i)
    reachable, intersection, distance = get_bead_intersection(r, dr, P(i), R_func(P(i)[0]))

    if reachable:
        # calculate the new dr (direction)
        reflected = True

        radial_vec = intersection - P(i)
        cos_value = (dr * radial_vec).sum() / norm(dr) / norm(radial_vec)

        # Correct for round-off error
        if cos_value > 1:
            cos_value = 1
        elif cos_value < -1:
            cos_value = -1

        incidence_angle = np.pi - np.arccos(cos_value)  # incidence angle
        # print('acos', dr, radial_vec, (dr * radial_vec).sum() / norm(dr) / norm(radial_vec))
        # print('incidence_angle, grad', incidence_angle / np.pi * 180)
        phi = np.pi - 2 * incidence_angle  # rotation angle
        reflection_angle = np.pi / 2 + np.arctan(radial_vec[1] / radial_vec[0])
        # print('reflection angle', reflection_angle)
        reflection_matrix = np.array(
            [[np.cos(2 * reflection_angle), np.sin(2 * reflection_angle)], [np.sin(2 * reflection_angle), -np.cos(2 * reflection_angle)]])

        dr_before_intersct = intersection - r
        dr_after_intersct = r + dr - intersection
        new_dr_after = reflection_matrix @ dr_after_intersct.transpose()
        # print('a', reflection_matrix, dr_after_intersct.transpose(), new_dr_after)
        internal_trajectory.append([*r, *dr_before_intersct])

        # Call function again to check for other reflections
        _, _, new_dr_after, internal_trajectory2 = get_reflections(
            intersection, new_dr_after, beads, R_func, Nx, Ny)
        new_r = intersection + new_dr_after
        internal_trajectory += internal_trajectory2

    else:
        # Check if reaching the border
        # print('Calling boundary intersection with ', r, dr, L)
        leaving, intersection, new_dr_after = get_boundary_intersection(r, dr, L)
        # print('boundary', r, dr, L, ">>>", leaving, intersection, new_dr_after)

        dr_temp = intersection - r
        # print('leave', leaving)
        if leaving:
            internal_trajectory.append([*r, *dr_temp])
            _, _, new_dr_after, internal_trajectory2 = get_reflections(
                intersection, new_dr_after, beads, R_func, Nx, Ny)
            new_r = intersection + new_dr_after
            internal_trajectory += internal_trajectory2
            # new_dr = new_r - r
        else:
            reflected = False
            new_r = r + dr
            internal_trajectory.append([*r, *dr])
        # new_dr = dr
    new_dr = new_r - r

    return reflected, new_r, new_dr, internal_trajectory


def get_bead_intersection(r, dr, P, R):
    """Get closest intersection between the trajectory and the given static bead.

    Input:
    P - bead location
    R - bead radius

    Returns:
    reachabel   -   boolean
    intersection    -   location of the intersection
    distance    -   distance to the intersection
    """

    reachable = False
    intersection, distance = [np.nan] * 2

    # Calculate the intersection point by taking into account the singularity of line description.
    # Use x or y substitutions when it gives a more numerically stable solution

    if np.abs(dr[0]) >= np.abs(dr[1]):
        method = 'y'
        # y = k*x + m
        k = dr[1] / dr[0]
        m = r[1] - k * r[0]

        # Quadratic equation
        a = 1
        b = (-2 * P[0] + 2 * k * (m - P[1])) / (1 + k**2)
        c = (P[0]**2 + (m - P[1])**2 - R**2) / (1 + k**2)
        discriminant = b**2 - 4 * a * c
        intersection_xs = (np.array([-1, 1]) * np.sqrt(discriminant) - b) / 2 / a
        intersections = np.array(list(zip(intersection_xs, intersection_xs * k + m)))
    else:
        method = 'x'
        # x = k*y + m
        k = dr[0] / dr[1]
        m = r[0] - k * r[1]

        # Quadratic equation
        a = 1
        b = (-2 * P[1] + 2 * k * (m - P[0])) / (1 + k**2)
        c = (P[1]**2 + (m - P[0])**2 - R**2) / (1 + k**2)
        discriminant = b**2 - 4 * a * c
        intersection_ys = (np.array([-1, 1]) * np.sqrt(discriminant) - b) / 2 / a
        intersections = np.array(list(zip(intersection_ys * k + m, intersection_ys)))

    # print(dr, method)
    # print('disc', discriminant)

    if discriminant > 0:
        # intersection_xs = (np.array([-1, 1]) * np.sqrt(discriminant) - b) / 2 / a
        # print('intersections', intersections)
        # print("%.20E" % intersection_xs[0], intersection_xs[0] * k + m, 0.5 * k + m)

        # get the closest intersection
        distances = np.sqrt(np.sum((intersections - r)**2, axis=1))

        # print('dist', [r, dr, P, R], distances)
        # Filter round-off errors
        filter = distances > atol
        distances = distances[filter]
        intersections = intersections[filter, :]
        # print('dist-f', distances)
        index = np.argmin(distances)
        intersection = intersections[index]
        # intersection = np.array([intersection_x, k * intersection_x + m])

        # print("Intersection", intersection)

        d_intersection = intersection - r
        distance = distances[index]
        # print('D', norm(d_intersection))
        # print('dint', d_intersection, dr)
        reachable = (
            d_intersection * dr).sum() >= 0 and norm(d_intersection) <= norm(dr) and norm(d_intersection) > 0
        # print('reach', reachable)

    return reachable, intersection, distance


def get_boundary_intersection(r, dr, L):
    """Calculate whether the particle is going to leave the box.
    Take into account only the closest reflection.
    Return a new starting location and displacement that takes into account the reflection.
    """

    # Check if r to r+dr goes through a wall, kx+m
    k = dr[1] / dr[0]
    m = r[1] - k * r[0]
    r_end = r + dr

    # Coordinates of the corner of crossing
    corner = (np.sign(dr) + 1) / 2 * L
    wall_intersections = np.array([[corner[0], corner[0] * k + m],
                                   [(corner[1] - m) / k, corner[1]]])
    distances_x = np.abs(r[0] - wall_intersections[:, 0])
    closest = np.argmin(distances_x)
    reachable = distances_x[closest] < np.abs(dr[0])

    if reachable:
        new_intersection = wall_intersections[closest, :]
        new_intersection[closest] -= np.sign(dr[closest]) * L
        new_r_end = r_end
        new_r_end[closest] -= np.sign(dr[closest]) * L
        new_dr = new_r_end - new_intersection
        leaving = True
    else:
        leaving = False
        new_intersection, new_dr = [np.nan] * 2

    return leaving, new_intersection, new_dr


def plot_beads(beads):
    patches = []
    for i in range(beads.shape[0]):
        bead = beads[i, :]
        patches.append(plt.Circle(bead[:2], bead[2], facecolor=None, edgecolor='r'))

    ax = plt.gca()
    ax.add_collection(matplotlib.collections.PatchCollection(patches))


def check_outside_beads(r, beads):
    """Return True if the point is outside all beads"""
    N = beads.shape[0]
    outside = True

    for i in range(N):
        bead = beads[i, :]
        # print(((bead[:2] - r)**2).sum(), bead[2]**2)
        if ((bead[:2] - r)**2).sum() <= bead[2]**2:
            outside = False
            break

    return outside


# %% Tests
if __name__ == '__main__':
    # get_reflections([0.1, 0.05], [-0.5, 0])
    # np.argmin([np.nan, np.nan])
    np.random.seed()
    ts, rs, beads, all_internal_trajectory = simulate()
    all_internal_trajectory = np.array(all_internal_trajectory)

    # % Plots
    # %matplotlib
    # max_N = np.min([32, rs.shape[0]])
    # max_N = rs.shape[0]
    plt.figure(1, clear=True)
    plt.scatter(rs[0, 0], rs[0, 1], color='g')
    plt.plot(rs[:, 0], rs[:, 1], 'r')
    # plt.plot(all_internal_trajectory[:, 0], all_internal_trajectory[:, 1], 'k')
    plt.axis('square')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plot_beads(beads)
    plt.show()

# list(zip([0, 1], [2, 3]))
print(rs.tolist())
