import numpy as np

a = 6378137.0
e2 = 0.00669437999013


def llh2ecef(llh):
    """
    Convert latitude [rad], longitude [rad], and height [m above ellipsoid]
    values into Earth Centered Earth Fixed (ECEF) [m] coordinates

    :param llh 1x3 numpy array of geodetic coordinates
    llh[0] = latitude (radians)
    llh[1] = longitude (radians)
    llh[2] = height (meters above ellipsoid)

    :return ecef 1x3 numpy array of ECEF coordinates [x y z in meters]

    """

    # print("llh =", llh)

    sin_lat = np.sin(llh[0])  # define for multiple reuse

    a_ke = a / np.sqrt(1 - e2 * sin_lat**2)

    rho = (a_ke + llh[2]) * np.cos(llh[0])

    ecef = np.array([rho * np.cos(llh[1]),
                     rho * np.sin(llh[1]),
                     (a_ke * (1 - e2) + llh[2]) * sin_lat])

    return ecef

origin = np.array([ 6.93928948e-01, -1.46945919e+00,  1.00000000e+03])
ecef = llh2ecef(origin)
print(ecef)
