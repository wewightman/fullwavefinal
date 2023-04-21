"""
This module is used to generate 3D intensity maps of muscle fibers 
with axis dependent scaling, angling, and offset

gengeneral must be passed points in the desired output space, as well
as the transformation parameters to transform a model from its given unit configuration
"""
import numpy as np

def gaussian(points):
    """Calculates unit gaussian kernel"""
    return np.exp(-np.sum(points**2, axis=0))

def sphere(points):
    """Returns 1 if within unit sphere, 0 otherwise"""
    vals = np.zeros(points.shape[1])
    inc = np.sum(points**2, axis=0) <= 1
    vals[inc] = 1
    return vals

def cube(points):
    """Returns 1 within unit cube, 0 otherwise"""

    # Calculate all points within 1 on each axis
    within_axis = np.abs(points) <= 1

    # Find points within 1 on all axes
    within = within_axis[0] & within_axis[1] & within_axis[2]

    # return mask
    vals = np.zeros(points.shape[1])
    vals[within] = 1
    return vals

def gengeneral(points, func, sigmas, thetas, x0):
    """find the intensity of the geometry of func at each point
    
    This function works via inversion of a U-space transformation.
    (IDK if thats real, I am just making shit up at this point)
    
    U describes the "unit" 3D space in x_u, y_u, and z_u with a geometry defined between with norm 1 less than 1

    Parameters:
    ----
    `points`: X space points to check bounds of. 3 by N
    `sigmas`: the scaling in x_u, y_u, and z_u
    `thetas`: the tilt about x_u, y_u, and z_u
    `x0`: the offset of the gaussian from X = 0
    `func`: a function that maps N points in R3 to R1 (3xN to N)
    """
    # scale x_u, y_u, and z_u
    T_sigma = np.array(
        [[sigmas[0], 0, 0, 0],
        [0, sigmas[1], 0, 0],
        [0, 0, sigmas[2], 0],
        [0, 0, 0, 1]]
    )

    # Calculate the rotation about the x axis
    theta = thetas[0]
    T_yz = np.array(
        [[1, 0, 0, 0],
        [0, np.cos(theta), -np.sin(theta), 0],
        [0, np.sin(theta), np.cos(theta), 0],
        [0, 0, 0, 1]]
    )

    # Calculate the rotation about the y axis
    theta = thetas[1]
    T_xz = np.array(
        [[np.cos(theta), 0, -np.sin(theta), 0],
        [0, 1, 0, 0],
        [np.sin(theta), 0, np.cos(theta), 0],
        [0, 0, 0, 1]]
    )

    # Calculate the rotation about the z axis
    theta = thetas[2]
    T_xy = np.array(
        [[np.cos(theta), -np.sin(theta), 0, 0],
        [np.sin(theta), np.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]]
    )

    # shift to new origin
    T_x0 = np.array(
        [[1, 0, 0, x0[0]],
        [0, 1, 0, x0[1]],
        [0, 0, 1, x0[2]],
        [0, 0, 0, 1]]
    )

    # matrix to transform u space (output space coordinates) into input space (X) coordinates
    T_ux = T_x0 @ T_xy @ T_xz @ T_yz @ T_sigma

    # Invert the transformation from U -> X to get the transform from X -> U
    T_xu = np.linalg.inv(T_ux)

    # Transform x space to u space
    points = np.concatenate((points, np.ones((1,points.shape[1]))), axis=0)
    points_u = (T_xu @ points)[:3,:]
    
    return func(points_u)
