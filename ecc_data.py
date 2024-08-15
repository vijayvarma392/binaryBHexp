import numpy as np
from matplotlib import pyplot as plt
import sxs
import h5py
from scipy.interpolate import InterpolatedUnivariateSpline

def spline_interp(newX, oldX, oldY, allowExtrapolation=False):
    """ Interpolates using splnes.
        If allowExtrapolation=True, extrapolates to zero.
    """
    if len(oldY) != len(oldX):
        raise Exception('Lengths dont match.')

    if not allowExtrapolation:
        if np.min(newX) - np.min(oldX) < -1e-5 \
                or np.max(newX) > np.max(oldX) > 1e-5:

            print(np.min(newX), np.min(oldX), np.max(newX), np.max(oldX))
            print(np.min(newX) < np.min(oldX))
            print(np.max(newX) > np.max(oldX))
            raise Exception('Trying to extrapolate, but '\
                'allowExtrapolation=False')

    if not np.all(np.diff(oldX) > 0):
        raise Exception('oldX must have increasing values')

    # returns 0 when extrapolating
    newY = InterpolatedUnivariateSpline(oldX, oldY, ext=1)(newX)
    return newY

def get_binary_data_NR(datapath):

    data = sxs.load(f'{datapath}/Strain_N2')
    t = data.t

    # Set t=0 at peak of (2,2) mode amplitude. FIXME: Eventually, use all modes for
    # this.
    peak_time = t[np.argmax(np.abs(data.data.T[4]))]
    t -= peak_time

    # Get waveform in dict format
    ellMax=4
    h_dict = {}
    idx = 0
    for ell in range(2, ellMax+1):
        for m in range(-ell, ell+1):
            h_dict[(ell, m)] = data.data.T[idx]
            idx += 1

    # spin h5file
    fspin = h5py.File(f'{datapath}/Horizons.h5', 'r')

    # Get the BH spins
    t_hor, chiAx, chiAy, chiAz = fspin["AhA.dir/chiInertial.dat"][()].T
    chiA = np.array([chiAx, chiAy, chiAz]).T
    _, chiBx, chiBy, chiBz = fspin["AhB.dir/chiInertial.dat"][()].T
    chiB = np.array([chiBx, chiBy, chiBz]).T
    t_hor -= peak_time  # Apply the same time shift as the waveform

    # Get the BH mass ratio
    _, mA = fspin['AhA.dir/ChristodoulouMass.dat'][()].T
    _, mB = fspin['AhB.dir/ChristodoulouMass.dat'][()].T
    # Get the value from the midpoint of the time series. Because masses don't
    # change significantly during the simulation, this is a good approximation.
    mA = mA[len(mA) // 2]
    mB = mB[len(mB) // 2]
    q = mA / mB

    # Get the coprecessing frame dynamics
    data_copr = data.to_coprecessing_frame()
    quat = data_copr.frame
    # Following Eq.3 of https://arxiv.org/abs/1905.09300
    orbphase = (np.unwrap(np.angle(data_copr.data.T[0])) - np.unwrap(np.angle(data_copr.data.T[4])))/4
    orbphase -= orbphase[0]

    # So far we have data on two different time grids:
    # t: the time grid of the waveform data
    # t_hor: the time grid of the horizons data

    # You shoud associate the following quantities with the ones in
    # https://github.com/vijayvarma392/binaryBHexp/blob/master/binaryBHexp#L498
    # Keep in mind that in the surrogate case, there are time arrays nr_sur.tds
    # and t_binary.
    
    # h_dict: Waveform dictionary on t. This is like h_nrsur.
    # chiA, chiB: Spin vectors on t_hor.
    # q: Mass ratio.
    # quat: Coprecessing frame quaternion on t.
    # orbphase: Orbital phase on t.

    ##### Instructions for Ryan
    # Step 1: Define a new time array t_binary that goes from t[0]+500 to t[-1]
    # in steps of 0.1. The first 500 is dropped to get rid of what is called
    # "junk radiation" during the initial stages where the simulation needs to
    # settle.

    """
    t_binary = np.arange(t[0]+500, t[-1], 0.1)
    """

    # Step 2: Interpolate h_dict, chiA, chiB, quat, and orbphase onto this
    # common time array. Use this function for interpolation:
    # https://github.com/vijayvarma392/binaryBHexp/blob/master/binaryBHexp#L280
    # Here, keep in mind that t_hor ends when the two BHs merge, while t
    # continues in the post-merger. So, to interpolate chiA and chiB onto
    # t_binary, you should set allowExtrapolation=True in spline_interp.
    # Also, there is a possibility that some of the indices for chiA/chiB have
    # nan values, so check for that first
    # (https://numpy.org/doc/stable/reference/generated/numpy.isnan.html) and
    # let me know if that happens.

    # chiA and chiB are in the frame of t_hor
    # h_dict, quat, and orbphase are in the frame of t

    """
    for i in h_dict.keys():
        h_dict[i] = spline_interp(t_binary, t, h_dict[i])
    
    orbphase = spline_interp(t_binary, t, orbphase)

    quat0 = spline_interp(t_binary, t, quat.ndarray.T[0])
    quat1 = spline_interp(t_binary, t, quat.ndarray.T[1])
    quat2 = spline_interp(t_binary, t, quat.ndarray.T[2])
    quat3 = spline_interp(t_binary, t, quat.ndarray.T[3])
    quat = np.array([quat0, quat1, quat2, quat3]).T

    chiA0 = spline_interp(t_binary, t_hor, chiA.T[0], allowExtrapolation = True)
    chiA1 = spline_interp(t_binary, t_hor, chiA.T[1], allowExtrapolation = True)
    chiA2 = spline_interp(t_binary, t_hor, chiA.T[2], allowExtrapolation = True)
    chiA = np.array([chiA0, chiA1, chiA2]).T

    chiB0 = spline_interp(t_binary, t_hor, chiB.T[0], allowExtrapolation = True)
    chiB1 = spline_interp(t_binary, t_hor, chiB.T[1], allowExtrapolation = True)
    chiB2 = spline_interp(t_binary, t_hor, chiB.T[2], allowExtrapolation = True)
    chiB = np.array([chiB0, chiB1, chiB2]).T
    """

    # Step 3: Once you have the data on t_binary, follow the same steps as in
    # get_binary_data() to compute L, BhA_traj, BhB_traj, and separation. Then
    # return everything in the same format. 

    """
    omega = get_omegaOrb_from_sparse_data(t_binary, orbphase)
    LHat = surfinBH._utils.lHat_from_quat(quat).T
    LMag = q/(1.+q)**2 * omega**(-1./3)
    L = LHat*LMag[:, None]

    separation = get_separation_from_omega(omega, mA, mB, chiA, \
        chiB, LHat)
    
    BhA_traj = get_trajectory(separation * mB, quat, orbphase, 'A')
    BhB_traj = get_trajectory(separation * mA, quat, orbphase, 'B')

    return t_binary, chiA, chiB, L, h_dict, BhA_traj, \
        BhB_traj, separation
    """
    
    # Step 4: Once all of this is done, you can switch get_binary_data() for
    # get_binary_data_NR() in the code and it should just work.

    QmAmB = np.array([q, mA, mB])

    
    np.savetxt("t.csv", t, delimiter = ",")
    np.savetxt("t_hor.csv", t_hor, delimiter = ",")
    np.savetxt("orbphase.csv", orbphase, delimiter = ",")
    np.savetxt("quat.csv", quat, delimiter = ",")
    np.savetxt("chiA.csv", chiA, delimiter = ",")
    np.savetxt("chiB.csv", chiB, delimiter = ",")
    np.savetxt("QmAmB.csv", QmAmB, delimiter = ",")

    for first in range(2, 5):
        for second in range(-first, first+1):
            np.savetxt("h_dict(" + str(first) + str(second) + ").csv", h_dict[(first, second)], delimiter = ",")
    

# This is an equal-mass, non-spinning simulation with ecc=0.45
datapath = 'SimulationAnnex/Private/AEI_Eccentric/BBH_SHK_q1_0_0_e01_D16/Lev3'
get_binary_data_NR(datapath)