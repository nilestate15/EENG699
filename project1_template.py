import scipy.io as spio
import numpy as np
import llh2ecef as l2e
import matplotlib.pyplot as plt
import math 

# CONSTANTS
A = 6378137.0  # meters
E = 0.0818191908426  # unitless
E2 = 0.00669437999013

# load mat file
mat = spio.loadmat('proj1_flight_trajectory.mat', squeeze_me=True)
trajectory = mat['K']

# Extract variable vectors from struct
t = trajectory['t'].tolist()      # [seconds]
lat = trajectory['lat'].tolist()  # [radians]
lon = trajectory['lon'].tolist()  # [radians]
h = trajectory['h'].tolist()      # [meters] (above ellipsoid)

# Stack list vectors into columns of LLH numPy array
LLH = np.vstack((lat, lon, h)).T  # LLH in 3 columns


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

    a_ke = A / np.sqrt(1 - E2 * sin_lat**2)

    rho = (a_ke + llh[2]) * np.cos(llh[0])

    ecef = np.array([rho * np.cos(llh[1]),
                     rho * np.sin(llh[1]),
                     (a_ke * (1 - E2) + llh[2]) * sin_lat])

    return ecef


# Geodetic to Local Level
# YOUR CODE HERE
# Setting origin of local-level frame
origin = LLH[0]
origin[2] = 0.
enuG2L = []

# Conversion Process from Geodetic to local-level frame (starts from index 1 due to index 0 being origin)
coords = LLH[1:]
for i in range(len(coords)):
    coord = coords[i]
    # differences between coordinates
    diff_lat = coord[0] - origin[0]
    diff_lon = coord[1] - origin[1]
    diff_h = coord[2] - origin[2]

    # Equivalent Radius
    Rm = (A * (1 - E**2)) / ((1 - E**2 * math.sin(coord[0])**2))**(3/2)
    Rp = (A) / ((1 - E**2 * math.sin(coord[0])**2))**(1/2)

    # Conversion to meters
    Pe = (Rp + coord[2]) * math.cos(coord[0]) * diff_lon
    Pn = (Rm + coord[2]) * diff_lat
    Pu = diff_h
    
    # Putting local-level coordinates in array
    PG = np.array([Pe, Pn, Pu])
    enuG2L += [PG]
enuG2L = np.array(enuG2L)


# Geodetic to ECEF to Local Level
# YOUR CODE HERE
# Converting Geodetic origin to ECEF origin
ecef_origin = llh2ecef(origin)
enuE2L = []

# DCM matrix
CGE = np.array([[-(math.sin(origin[1])), math.cos(origin[1]), 0.0],
                    [-(math.sin(origin[0]) * math.cos(origin[1])), -(math.sin(origin[0]) * math.sin(origin[1])),  math.cos(origin[0])],
                    [(math.cos(origin[0]) * math.cos(origin[1])), (math.cos(origin[0]) * math.sin(origin[1])), math.sin(origin[0])]])


# Conversion Process from ECEF to local-level frame (starts from index 1 due to index 0 being origin)
for n in range(len(coords)):
    # Converting Geodetic Coords to ECEF
    coord_4e = coords[n]
    ecef_coord = llh2ecef(coord_4e)

    # Finding deltas of X,Y,Z
    diff_x = ecef_coord[0] - ecef_origin[0]
    diff_y = ecef_coord[1] - ecef_origin[1]
    diff_z = ecef_coord[2] - ecef_origin[2]

    # Making PE array
    PE = np.array([diff_x, diff_y, diff_z]).T

    # Solving to find PG for coord
    ecef_PG = np.dot(CGE, PE)

    # Adding to array 
    enuE2L += [ecef_PG]
enuE2L = np.array(enuE2L)




# Plotting
# YOUR CODE HERE
# Plot 1 (Horizontal Position)
plot1 = plt.figure(1)
plt.title('Horizontal Position')
plt.plot(enuG2L[:,0], enuG2L[:,1], marker = '.', markersize = 5)
plt.plot(enuE2L[:,0], enuE2L[:,1], marker = '.', markersize = 5)
plt.legend(['enuG2L origin', 'enuE2L origin', 'enuG2L', 'enuE2L'])
plt.xlabel('Easting (m)')
plt.ylabel('Northing (m)')
plt.grid(linestyle = '--', linewidth = 0.5)

# Plot 2 (Altitude vs Time)
t = np.array (range(1, len(enuE2L)+1))

plot2 = plt.figure(2)
plt.title('Altitude vs Time')
plt.plot(t, enuG2L[:,2], marker = '.', markersize = 5)
plt.plot(t, enuE2L[:,2], marker = '.', markersize = 5)
plt.legend(['enuG2L', 'enuE2L'])
plt.xlabel('Time (s)')
plt.ylabel('Altitude (m)')
plt.grid(linestyle = '--', linewidth = 0.5)

# Plot 3 (Delta of enuE2L and enuG2L)
fig1, (ax1, ax2, ax3) = plt.subplots(3)

delta_e = enuE2L[:,0] - enuG2L[:,0]
delta_n = enuE2L[:,1] - enuG2L[:,1]
delta_u = enuE2L[:,2] - enuG2L[:,2]

fig1.suptitle('Deltas of enuE2L and enuEGL vs Time')
ax1.plot(t, delta_e, marker = '.', markersize = 5)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Delta Easting (m)')
ax2.plot(t, delta_n, marker = '.', markersize = 5)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Delta Northing (m)')
ax3.plot(t, delta_u, marker = '.', markersize = 5)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Delta Up (m)')


plt.show()  # show all plots


