

import numpy as np
from scipy.spatial import Voronoi


def get_exterior_cells(cells):
    """
    The function returns the indices of the outer layer of cells in the data region.
    The detection is performed by i) identifying voronoi regions with abstract (diverging) vertices - divergent regions and corresponding divergent cells; and ii) regions that have at least two neighboring divergent regions.
    """

    def get_vertices_for_cell(cell):
        return cells.tessellation.cell_vertices[cell]

    cells_len = cells.location_count.size

    # Get corresponding Voronoi regions
    voronoi = Voronoi(cells.tessellation.cell_centers)

    divergent_cells = []
    for cell_ind in range(cells_len):
        vor_region = voronoi.point_region[cell_ind]
        vor_vertices = np.asarray(voronoi.regions[vor_region])

        # Check if any vertex is abstract (-1)
        if np.any(vor_vertices < 0):
            divergent_cells.append(cell_ind)

    # Create a surface cells list, which extends the divergent list by cells in contact with at least 2 divergent cells
    surface_cells = divergent_cells.copy()
    for cell_ind in range(cells_len):
        # Find neighbors (cells that share more than 1 vertex)
        cell_vertices = get_vertices_for_cell(cell_ind)
        cell_neighbors = [i for i in range(cells_len) if i is not cell_ind and len(
            set(get_vertices_for_cell(i)) & set(cell_vertices)) > 0]

        # Identify divergent neighbors
        divergent_neighbors = set(divergent_cells) & set(cell_neighbors)
        if len(divergent_neighbors) > 1:
            surface_cells.append(cell_ind)

    # Remove duplicates
    surface_cells = list(set(surface_cells))

    return surface_cells
