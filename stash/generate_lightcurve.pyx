cimport cython
import numpy as np
cimport numpy as np
from copy import copy

__all__ = ['generate_lightcurve']

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

DTYPE = np.float32

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def rebin(arr, R_planet_pixels_upper, supersample_factor):
    if supersample_factor == 1:
        return arr

    cdef list new_shape = [2*R_planet_pixels_upper, 2*R_planet_pixels_upper]
    cdef list shape = [new_shape[0], arr.shape[0] // new_shape[0],
                       new_shape[1], arr.shape[1] // new_shape[1]]
    reshaped = arr.reshape(shape).sum(-1).sum(1)
    return reshaped / reshaped.max()


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def generate_lightcurve(image, b_pixels, R_planet_pixels=17.26, background=269,
                        supersample_factor=4):
    """

    Parameters
    ----------
    image : `~np.ndarray`
        SDO HMI image
    R_planet_pixels : float
        Radius of planet in pixels
    background : float
        Backgorund flux in counts
    supersample_factor : int
        Supersample the planet pixelation by a factor ``supersample_factor``

    Returns
    -------
    fluxes : `~np.ndarray`
        Array of fluxes
    """
    cdef int R_planet_pixels_upper = int((R_planet_pixels+1)//1)
    cdef float R_planet_pixels_squared = R_planet_pixels**2
    cdef int n_images = 4096 - 2 * R_planet_pixels_upper
    cdef int step, i, j
    cdef np.ndarray mask = np.zeros((2*R_planet_pixels_upper,
                                     2*R_planet_pixels_upper), dtype=DTYPE)
    cdef np.ndarray mask_supersampled = np.zeros((2 * R_planet_pixels_upper * supersample_factor,
                                                  2 * R_planet_pixels_upper * supersample_factor), dtype=DTYPE)

    cdef np.ndarray fluxes = np.zeros(n_images, dtype=DTYPE)
    cdef float unobscured_flux, obscured_flux, r_pixel

    cdef int planet_center_y = 1680 - b_pixels
    cdef int planet_center_x = 4096 - R_planet_pixels_upper

    cdef int planet_centroid_x, planet_centroid_y

    cdef int y_start, y_end, x_start, x_end
    cdef int delta_y_lower = 0
    cdef int delta_y_upper = 2*R_planet_pixels_upper
    cdef float sum_flux = 0

    unobscured_flux = image.sum()

    cdef int new_shape = 2*R_planet_pixels_upper * supersample_factor
    cdef list shape = [new_shape, mask.shape[0] // new_shape,
                       new_shape, mask.shape[1] // new_shape]

    for i in range(0, 2*R_planet_pixels_upper * supersample_factor):
        for j in range(0, 2*R_planet_pixels_upper * supersample_factor):
            r_pixel = (i - R_planet_pixels_upper * supersample_factor)**2 + (j - R_planet_pixels_upper * supersample_factor)**2
            if r_pixel < (R_planet_pixels * supersample_factor)**2:
                mask_supersampled[i, j] += 1.0

    mask = rebin(mask_supersampled, R_planet_pixels_upper, supersample_factor)

    for step in range(n_images):
        x_start = planet_center_x - step - R_planet_pixels_upper
        x_end = planet_center_x - step + R_planet_pixels_upper

        y_start = planet_center_y - R_planet_pixels_upper
        y_end = planet_center_y + R_planet_pixels_upper

        if y_start < 0:
            delta_y_lower = copy(abs(y_start))
            y_start = 0

        if y_end > 4096:
            delta_y_upper = copy(delta_y_upper - 4096)
            y_end = 4096

        obscured_flux = np.sum(image[y_start:y_end, x_start:x_end] *
                               mask[delta_y_lower:delta_y_upper, :])

        fluxes[step] = (unobscured_flux - obscured_flux)

    return fluxes
