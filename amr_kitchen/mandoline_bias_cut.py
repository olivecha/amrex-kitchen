import sys
import multiprocessing
import meshio
import skimage
import scipy
from scipy.interpolate import NearestNDInterpolator
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import expand_array3d

FIELD_ID = None
LV = None
plane_fun = None

def plane_generator(plane_point: np.ndarray[float],
                    plane_normal: np.ndarray[float]) -> callable:
    """
    This creates a function returning the distance from a plane in 3D
    for arrays of points with shape (n, 3)
    plane_point: a point on the target plane
    pane_normal: normal of the plane
    """
    def plane(coords):
        return np.dot(coords - plane_point, plane_normal)
    return plane


def box_vertices(box: np.ndarray[float]) -> np.ndarray[float]:
    """
    This returns an array with shape (8, 3) of the corner
    points of an AMR box using the faces coordinates
    box: array with shape (3, 2) of the box faces coordinates
    """
    return np.array([[box[0][0], box[1][0], box[2][0]],
                     [box[0][0], box[1][0], box[2][1]],
                     [box[0][0], box[1][1], box[2][0]],
                     [box[0][0], box[1][1], box[2][1]],
                     [box[0][1], box[1][0], box[2][0]],
                     [box[0][1], box[1][0], box[2][1]],
                     [box[0][1], box[1][1], box[2][0]],
                     [box[0][1], box[1][1], box[2][1]]])

def check_intersect(box: np.ndarray[float],
                    plane_fun: callable) -> bool:
    """
    This return true if the plane intersect a box defined
    by its corner vertices.
    box: box vertices with shape (8, 3)
    plane_fun: function returning the distance from a plane
    """
    peq = plane_fun(box_vertices(box))
    if (np.min(peq) < 0) and (np.max(peq) > 0):
        return True
    else:
        return False

# Read the data at the finest level
def mp_slice_box(args):
    pck = args[0]
    bid = args[1]
    lv = 4
    data = pck[FIELD_ID][lv][bid]
    x, y, z = pck.box_points(lv, bid)
    P = plane_fun(np.array([x, y, z]).T).T
    if not (np.max(P) > 0 and np.min(P) < 0):
        return
    points, faces, normals, values = skimage.measure.marching_cubes(P, level=0)
    plane_data = scipy.ndimage.map_coordinates(data, points.T, order=1, mode='nearest')
    indices = pck.cells[lv]['indexes'][bid]
    points += indices[0]
    points[:, 0] /= (pck.grid_sizes[lv][0] - 1)
    points[:, 1] /= (pck.grid_sizes[lv][1] - 1)
    points[:, 2] /= (pck.grid_sizes[lv][2] - 1)
    points[:, 0] *= (pck.geo_high[0] - pck.dx[lv][0]/2 - pck.geo_low[0] + pck.dx[lv][0]/2) 
    points[:, 0] += pck.geo_low[0] + pck.dx[lv][0]/2
    points[:, 1] *= (pck.geo_high[1] - pck.dx[lv][1]/2 - pck.geo_low[1] + pck.dx[lv][1]/2) 
    points[:, 1] += pck.geo_low[1] + pck.dx[lv][1]/2
    points[:, 2] *= (pck.geo_high[2] - pck.dx[lv][2]/2 - pck.geo_low[2] + pck.dx[lv][2]/2) 
    points[:, 2] += pck.geo_low[2] + pck.dx[lv][2]/2
    return data, points, faces

def main():
# INPUTS
    plotfile = sys.argv[1]
    field = sys.argv[2]
    pck = PlotfileCooker(plotfile, ghost=True)
    plane_point = [0.00, 0.0015 + 3 * 7.8125e-05, 0.00]
    plane_normal = [0, -1, 0.577]
    global FIELD_ID
    FIELD_ID = pck.fields[field]

    if np.count_nonzero(np.isclose(plane_normal, 0)) == 0:
        raise ValueError("Plane must be rotated around a single axis")

    # Define the plane function
    global plane_fun
    plane_fun = plane_generator(plane_point, plane_normal)

    # Find all boxes indexes intersecting with plane
    intersect_indices = []
    print("Finding boxes intersecting the plane")
    n_boxes = np.sum([len(pck.boxes[lv]) for lv in range(pck.limit_level + 1)])
    pbar = tqdm(total=n_boxes)
    for lv in range(pck.limit_level + 1):
        lv_indices = []
        for idx, box in enumerate(pck.boxes[lv]):
            pbar.update(1)
            if check_intersect(box, plane_fun):
                lv_indices.append(idx)
        intersect_indices.append(lv_indices)
    pbar.close()

    # Filter out lower level boxes where the region 
    # intersecting with the plane is covered by higher level data
    covering_masks = []
    filtered_indices = []
    print("Filtering out boxes covered by upper level data")
    n_iter = np.sum([len(arr) for arr in intersect_indices[:-1]])
    pbar = tqdm(total=n_iter)
    for lv in range(pck.limit_level): # Last level is not masked
        # use box indices in list
        lv_indices = []
        lv_masks = []
        next_lv_factors = pck.grid_sizes[lv + 1] // pck.box_arrays[lv + 1].shape
        for box_id in intersect_indices[lv]:
            pbar.update(1)
            indices = pck.cells[lv]['indexes'][box_id]
            # Convert to box array indices at lv + 1
            barr_starts = np.array((indices[0] * 2) // next_lv_factors, dtype=int)
            barr_ends = np.array((indices[1] * 2) // next_lv_factors, dtype=int)
            # Get all the intersecting box ids at lv + 1
            next_level_boxes = pck.box_arrays[lv + 1][barr_starts[0]:barr_ends[0]+1,
                                                      barr_starts[1]:barr_ends[1]+1,
                                                      barr_starts[2]:barr_ends[2]+1]
            # Convert the upper level box slice to lower level bool
            bcast_factor = next_lv_factors[0] // 2
            next_lv_map = expand_array3d(next_level_boxes, bcast_factor)
            mask = np.zeros_like(next_lv_map, dtype=bool)
            # mask is true where the value is -1 (means no upper level data)
            mask[next_lv_map == -1] = 1
            # Add only partially or not covered lower level boxes
            if not np.all(~mask):
                x, y, z = pck.box_points(lv, box_id)
                plane_eval = plane_fun(np.vstack([x[mask], y[mask], z[mask]]).T)
                # Evaluate the plane function at the intersecting box
                # and only keep the box if the plane slices it
                if (np.min(plane_eval) < 0) and (np.max(plane_eval) > 0):
                    lv_indices.append(box_id)
                    lv_masks.append(mask)
        filtered_indices.append(lv_indices)
        covering_masks.append(lv_masks)
    pbar.close()
    # Add all intersect indices at limit level
    filtered_indices.append(intersect_indices[pck.limit_level])

    # Read the data at the plane points
    all_points = []
    all_data = []
    all_faces = []
    all_levels = []
    point_count = 0
    # Up until the finest level
    print("Constructing the plane surface")
    for lv in range(pck.limit_level):
        # For each intersecting box at the level
        if len(filtered_indices[lv]) == 0:
            continue
        print(f"Level {lv}")
        for bid, data, msk in tqdm(zip(filtered_indices[lv],
                                       pck[FIELD_ID][lv].iter(filtered_indices[lv]),
                                       covering_masks[lv]),
                             total=len(filtered_indices[lv])):
            # Evaluate the box points in 3D
            x, y, z = pck.box_points(lv, bid)
            # Plane function at the box points
            P = plane_fun(np.array([x, y, z]).T).T
            # Find the plane surface indices
            points, faces, normals, values = skimage.measure.marching_cubes(P, level=0)
            # Interpolate the mask on the surface
            mask_data = scipy.ndimage.map_coordinates(msk, points.T, order=1, cval=1, mode='constant')
            # Interpolate the data on the surface
            plane_data = scipy.ndimage.map_coordinates(data, points.T, order=1, cval=1, mode='constant')
            # Remove masked points from the surface
            removed_points = np.arange(points.shape[0])[mask_data]
            faces_mask = np.any(np.isin(faces, removed_points), axis=1)
            faces = faces[faces_mask, :]
            # Convert box indices to global coordinates
            indices = pck.cells[lv]['indexes'][bid]
            points += indices[0]
            points[:, 0] /= (pck.grid_sizes[lv][0] - 1)
            points[:, 1] /= (pck.grid_sizes[lv][1] - 1)
            points[:, 2] /= (pck.grid_sizes[lv][2] - 1)
            points[:, 0] *= (pck.geo_high[0] - pck.dx[lv][0]/2 - pck.geo_low[0] + pck.dx[lv][0]/2) 
            points[:, 0] += pck.geo_low[0] + pck.dx[lv][0]/2
            points[:, 1] *= (pck.geo_high[1] - pck.dx[lv][1]/2 - pck.geo_low[1] + pck.dx[lv][1]/2) 
            points[:, 1] += pck.geo_low[1] + pck.dx[lv][1]/2
            points[:, 2] *= (pck.geo_high[2] - pck.dx[lv][2]/2 - pck.geo_low[2] + pck.dx[lv][2]/2) 
            points[:, 2] += pck.geo_low[2] + pck.dx[lv][2]/2
            all_points.append(np.transpose(points))
            all_faces.append(np.transpose(faces + point_count))
            point_count += points.shape[0]
            all_levels.append(np.ones(points.shape[0]) * lv)
            all_data.append(plane_data)



    lv = 4
    print(f"Level {lv}")
    pool = multiprocessing.Pool()
    for bid, data in tqdm(zip(filtered_indices[lv],
                              pck[field][lv].iter(filtered_indices[lv])),
                          total=len(filtered_indices[lv])):
        x, y, z = pck.box_points(lv, bid)
        P = plane_fun(np.array([x, y, z]).T).T
        if not (np.max(P) > 0 and np.min(P) < 0):
            continue
        points, faces, normals, values = skimage.measure.marching_cubes(P, level=0)
        plane_data = scipy.ndimage.map_coordinates(data, points.T, order=1, mode='nearest')
        all_data.append(plane_data)
        indices = pck.cells[lv]['indexes'][bid]
        points += indices[0]
        points[:, 0] /= (pck.grid_sizes[lv][0] - 1)
        points[:, 1] /= (pck.grid_sizes[lv][1] - 1)
        points[:, 2] /= (pck.grid_sizes[lv][2] - 1)
        points[:, 0] *= (pck.geo_high[0] - pck.dx[lv][0]/2 - pck.geo_low[0] + pck.dx[lv][0]/2) 
        points[:, 0] += pck.geo_low[0] + pck.dx[lv][0]/2
        points[:, 1] *= (pck.geo_high[1] - pck.dx[lv][1]/2 - pck.geo_low[1] + pck.dx[lv][1]/2) 
        points[:, 1] += pck.geo_low[1] + pck.dx[lv][1]/2
        points[:, 2] *= (pck.geo_high[2] - pck.dx[lv][2]/2 - pck.geo_low[2] + pck.dx[lv][2]/2) 
        points[:, 2] += pck.geo_low[2] + pck.dx[lv][2]/2
        all_points.append(np.transpose(points))
        all_faces.append(np.transpose(faces + point_count))
        point_count += points.shape[0]
        all_levels.append(np.ones(points.shape[0]) * lv)

    all_points = np.hstack(all_points).T
    print("Projecting the slice to 2D coordinates")
    # Rotation axis
    e1 = [1 if i == 0 else 0 for i in plane_normal]
    rot_axis = np.arange(3)[e1]
    # Second vector for the basis
    e2 = np.cross(e1, plane_normal)
    # Project points into the 2D basis
    points_e1 = np.dot(all_points, e1)
    points_e2 = np.dot(all_points, e2)
    # Unform grid to which the data is interpolated
    x_plane = np.linspace(points_e1.min(),
                          points_e1.max(),
                          4096)
    # Keep the mesh resolution for the stretched coordinate
    # TODO: handle cases where dx =/= dy =/= dz by finding
    # the axis the most aligned with the second plane
    # coordinate vector
    # Use absolute value so the second vector of the
    # basis gives positive coordinates
    y_plane = np.arange(np.abs(points_e2).min(),
                        np.abs(points_e2).max() + 1e-16,
                        7.8125e-05)
    print("Interpolating on uniform grid")
    itrp = NearestNDInterpolator(np.transpose([points_e1,
                                               np.abs(points_e2)]),
                                 np.hstack(all_data))

    d_plane = itrp(*np.meshgrid(x_plane, y_plane))
    # Save key for the plane rotation axis
    coords = np.array(['x', 'y', 'z'])
    rot_key = coords[np.flatnonzero(e1)][0]
    # Concatenate the two other coordinates
    pcoord_indices = ~np.array(e1, dtype=bool)
    pcoord_key = ''.join(coords[pcoord_indices])
    output = {rot_key:x_plane,
              pcoord_key:y_plane,
              field:d_plane}
    print("Saving 2D array")
    np.savez(f"{plotfile}_{field}_{rot_key}", **output)

    #cells = {"triangle":np.hstack(all_faces).T}
    #mesh = meshio.Mesh(all_points, cells,
    #                   point_data={"level":np.hstack(all_levels),
    #                               field:np.hstack(all_data)})
    #mesh.write(f"{plotfile}_{field}_bias.vtk")

if __name__ == "__main__":
    main()
