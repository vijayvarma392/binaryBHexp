desc = """Animations of binary black hole scattering.
Generates an animation of a binary black hole merger and the final remnant.

Example usage:
python black_hole_scattering.py --q 2 --chiA 0.2 0.7 -0.1 --chiB 0.2 0.6 0.1

Note: Time values displayed in the plot are non-uniform and non-linear:
During the inspiral there are 30 frames per orbit.
After the merger each frame corresponds to a time step of 100M.

The precessing frame quaternion and phase is obtained from NRSur7dq2.
The separation is estimated from 3.5PN relation between r and omega, and omega
is obtained from NRSur7dq2.
The remnant properties are obtained from surfinBH.

Links:
NRSur7dq2: https://pypi.org/project/NRSur7dq2/
surfinBH: https://pypi.org/project/surfinBH/
"""

import numpy as np
import matplotlib.pyplot as P
import argparse

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

import surfinBH
import NRSur7dq2
import harmonics

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
P.style.use('seaborn')

use_palettable = False
try:
    import palettable
    use_palettable = True
except:
    print("palettable not found, using regular old colors. Get palettable" \
        + " with 'pip install palettable' to get fancy colors.")

if use_palettable:
    from palettable.wesanderson import Aquatic1_5
    from palettable.wesanderson import Darjeeling2_5
    from palettable.wesanderson import Darjeeling3_5
    from palettable.wesanderson import FantasticFox2_5
    from palettable.wesanderson import GrandBudapest5_5
    from palettable.wesanderson import GrandBudapest1_4
    from palettable.wesanderson import GrandBudapest4_5
    from palettable.wesanderson import GrandBudapest3_6
    from palettable.wesanderson import Zissou_5
    from palettable.wesanderson import Royal1_4
    colors_aq_15 = Aquatic1_5.mpl_colors
    colors_dj_25 = Darjeeling2_5.mpl_colors
    colors_dj_35 = Darjeeling3_5.mpl_colors
    colors_ff_25 = FantasticFox2_5.mpl_colors
    colors_gb_55 = GrandBudapest5_5.mpl_colors
    colors_gb_14 = GrandBudapest1_4.mpl_colors
    colors_gb_45 = GrandBudapest4_5.mpl_colors
    colors_gb_36 = GrandBudapest3_6.mpl_colors
    colors_zs_5 = Zissou_5.mpl_colors
    colors_ry_14 = Royal1_4.mpl_colors

    colors_dict = {
            'BhA_traj': colors_zs_5[2],
            'BhB_traj': colors_gb_14[2],
            'BhA_spin': colors_gb_45[3],
            'BhB_spin': colors_aq_15[1],
            'BhC_spin': colors_zs_5[0],
            'L': colors_gb_55[1],
            'J': colors_aq_15[0],
            'info': colors_aq_15[0],
            'h+': colors_dj_25[3],
            'hx': colors_dj_25[1],
            }
else:

    colors_dict = {
            'BhA_traj': 'indianred',
            'BhB_traj': 'rebeccapurple',
            'BhA_spin': 'goldenrod',
            'BhB_spin': 'steelblue',
            'BhC_spin': 'forestgreen',
            'L': 'orchid',
            'info': 'k',
            'h+': 'tomato',
            'hx': 'steelblue',
            }



# Make very low def video. This is needed for pypi.
LOW_DEF = False

# number of frames per orbit
PTS_PER_ORBIT = 30
if LOW_DEF:
    PTS_PER_ORBIT = 15

# Time at which to freeze video for 5 seconds
FREEZE_TIME = -100


zorder_dict = {
        'contourf': -200,
        'info_text': 200,
        'notice_text': 100,
        'traj': 100,
        'Bh': 110,
        'spin': 200,
        'L': 110,
        'J': 90,
        }


class Arrow3D(FancyArrowPatch):
    def __init__(self, vecs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = vecs

    def set_arrow_at_origin(self, vec):
        xs = [0, vec[0]]
        ys = [0, vec[1]]
        zs = [0, vec[2]]
        self._verts3d = xs, ys, zs

    def set_BH_spin_arrow(self, Bh_loc, mass, chi_vec, scale_factor=20):
        """
        The length of the arrow is proportinal to Kerr parameter (a) of the BH.
        """
        x, y, z =  Bh_loc
        u, v, w =  chi_vec*mass
        xs = [x, x+u*scale_factor]
        ys = [y, y+v*scale_factor]
        zs = [z, z+w*scale_factor]
        self._verts3d = xs, ys, zs

    def set_angular_momentum_arrow(self, L, scale_factor=25):
        #Place the base at center.
        xs = [0, L[0]*scale_factor]
        ys = [0, L[1]*scale_factor]
        zs = [0, L[2]*scale_factor]
        self._verts3d = xs, ys, zs

    def reset(self):
        self._verts3d = None

    def draw(self, renderer):
        if self._verts3d is not None:
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
            FancyArrowPatch.draw(self, renderer)


#----------------------------------------------------------------------------
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

#----------------------------------------------------------------------------
def get_trajectory(separation, quat_nrsur, orbphase_nrsur, bh_label):
    """ Gets trajectory of a component BH in a binary given the separation,
    the coprecessing frame quaternion and orbital phase in the coprecessing
    frame.
    """
    if bh_label == 'A':
        offset = 0
    else:
        offset = np.pi

    x_copr = separation * np.cos(orbphase_nrsur+offset)
    y_copr = separation * np.sin(orbphase_nrsur+offset)
    z_copr = np.zeros(len(x_copr))

    Bh_traj_copr = np.array([x_copr, y_copr, z_copr])
    Bh_traj = surfinBH._utils.transformTimeDependentVector(quat_nrsur, \
        Bh_traj_copr, inverse=0)

    return Bh_traj


#-----------------------------------------------------------------------------
def get_uniform_in_orbits_times(t, phi_orb, PTS_PER_ORBIT):
    """
    returns sparse time array such that there are PTS_PER_ORBIT points
    in each orbit.
    """
    # get numer of orbits
    n_orbits = int(abs((phi_orb[-1] - phi_orb[0])/(2*np.pi)))

    # get sparse times such that there are PTS_PER_ORBIT points in each orbit
    n_pts = int(n_orbits*PTS_PER_ORBIT)
    phi_orb_sparse = np.linspace(phi_orb[0], phi_orb[-1], n_pts)
    t_sparse = np.interp(phi_orb_sparse, phi_orb, t)

    return t_sparse

#----------------------------------------------------------------------------
def get_omegaOrb_from_sparse_data(t_sparse, phiOrb_sparse):
    """ Computes orbital frequency from sparse data using splines.
    """
    # spline interpolant for phase
    phiOrb_spl = UnivariateSpline(t_sparse, phiOrb_sparse, s=0)

    # spline for phase derivative
    omegaOrb_spl = phiOrb_spl.derivative()

    return omegaOrb_spl(t_sparse)

#----------------------------------------------------------------------------
def get_marker_size(mass, chi):
    """ Marker size proportional to the Kerr horizon radius """
    if LOW_DEF:
        scale_factor = 15
    else:
        scale_factor = 30

    chimag = np.sqrt(np.sum(chi**2))
    rplus = mass + mass*np.sqrt(1 - chimag**2)
    return rplus * scale_factor

#----------------------------------------------------------------------------
def get_separation_from_omega(omega, mA, mB, chiA, chiB, LHat, pnorder=3.5):
    """ Roughly 3.5 PN accurate separation. This is not verified or tested,
    so don't use this for real science, only visualization. """

    eta = mA*mB
    deltaM = mA - mB

    Sigma_vec = mB*chiB - mA*chiA
    S_vec = mA**2.*chiA + mB**2.*chiB

    # some dot products
    chiAL = np.sum(LHat*chiA, axis=1)
    chiBL = np.sum(LHat*chiB, axis=1)
    chiAB = np.sum(chiA*chiB, axis=1)
    SigmaL = np.sum(Sigma_vec*LHat, axis=1)
    SL = np.sum(S_vec*LHat, axis=1)


    # Get 3.5 PN accurate gamma=1./r from Eq.(4.3) of
    # https://arxiv.org/pdf/1212.5520v2.pdf, but ignore the
    # log term in x**3 term


    x = omega**(2./3.)
    gamma_by_x = 0

    if pnorder >= 0:
        gamma_by_x += 1
    if pnorder >= 1:
        gamma_by_x += x * (1. - 1./3 *eta)
    if pnorder >= 1.5:
        gamma_by_x += x**(3./2) * (5./3 * SL + deltaM * SigmaL )
    if pnorder >= 2:
        gamma_by_x += x**2 * (1 - 65./12 *eta)
    if pnorder >= 2.5:
        gamma_by_x += x**(5./2) * ( (10./3 + 8./9 * eta)*SL \
            + 2* deltaM * SigmaL)
    if pnorder >= 3:
        gamma_by_x += x**3 * (1. + (-2203./2520 -41./192 * np.pi**2)*eta \
                + 229./36 * eta**2 + 1./81 * eta**3)
    if pnorder >= 3.5:
        gamma_by_x += x**(7./2) * ( (5 - 127./12 *eta - 6 * eta**2)*SL + \
                deltaM * SigmaL * (3 - 61./6 *eta - 8./3 * eta**2) )


    r = 1./gamma_by_x/x

    # To this add the 2PN spin-spin term from Eq.(4.13) of
    # https://arxiv.org/pdf/gr-qc/9506022.pdf
    if pnorder >= 2:
        r += omega**(-2./3) * (-1./2 * eta * chiAB) * omega**(4./3)

    return r


#----------------------------------------------------------------------------
def get_grids_on_planes(num_pts_1d, max_range):
    # generate grid
    x_1d = np.linspace(-max_range, max_range, num_pts_1d)
    y_1d = np.linspace(-max_range, max_range, num_pts_1d)
    z_1d = np.linspace(-max_range, max_range, num_pts_1d)

    [xZ, yZ] = np.meshgrid(x_1d, y_1d)
    [xY, zY] = np.meshgrid(x_1d, z_1d)
    [yX, zX] = np.meshgrid(y_1d, z_1d)

    xX = zZ = -max_range
    yY = max_range

    # Get Euclidean radii and th,ph
    rZ = np.sqrt(xZ**2 + yZ**2 + zZ**2)
    thZ = np.arccos(zZ/rZ)
    phZ = np.arctan2(yZ,xZ)

    rY = np.sqrt(xY**2 + yY**2 + zY**2)
    thY = np.arccos(zY/rY)
    phY = np.arctan2(yY,xY)

    rX = np.sqrt(xX**2 + yX**2 + zX**2)
    thX = np.arccos(zX/rX)
    phX = np.arctan2(yX,xX)

    return [rX, thX, phX], [yX, zX], \
           [rY, thY, phY], [xY, zY], \
           [rZ, thZ, phZ], [xZ, yZ]

#----------------------------------------------------------------------------
def get_waveform_on_grid(t_vals, t_idx, h_dict, sph_grid):
    """ Compute absolute value of strain at each r, th, ph value, using
    the retarded time.
    """
    r, th, ph = sph_grid
    h = np.zeros(r.shape, dtype=complex)
    # find the time index that's closest to t_ret = t-r
    t = t_vals[t_idx]
    t_ret_idx = np.vectorize(lambda r: np.argmin(np.abs(t_vals - t + r)))(r)
    for key in h_dict.keys():
        ell, m = key
        ylm = np.vectorize(harmonics.sYlm)(-2, ell, m, th, ph)
        h += h_dict[key][t_ret_idx]*ylm
    return np.real(h/r)

#----------------------------------------------------------------------------
def get_waveform_timeseries(h_dict, azim, elev):
    """ Compute the timeseries to plot in the lower panel from a given viewpoint
    """
    ph = azim * np.pi/180.
    th = (90. - elev) * np.pi/180.
    h = np.zeros_like(h_dict[h_dict.keys()[0]], dtype=complex)
    for key in h_dict.keys():
        ell, m = key
        ylm = harmonics.sYlm(-2, ell, m, th, ph)
        h += h_dict[key]*ylm
    return h


#----------------------------------------------------------------------------
def get_camera_trajectory(t_binary, period=1000, stop_time=-500, azim_shift=1):
    """ Get predetermined trajectory for camera.
        Varies the elevation between 0 and 90 degrees with given period.
        Varies the azimuthal angle by azim_shift degrees every time step.
        Ensures that elevation = 30 and azimuthal angle = -60 at stop_time, the
        defaults set by matplotlib. update_lines should stop moving the camera
        at this time.
    """
    omega = np.pi/period   # notice the sin**2

    stop_idx = np.argmin(np.abs(t_binary - stop_time))
    t = t_binary - t_binary[stop_idx]

    # phase shift so that we get elev=1./3 * 90 at stop_idx
    phi0 = np.arcsin(1./np.sqrt(3))

    # sin**2 is smoother than abs(sin)
    elev_vec = 90*np.sin(omega*t + phi0)**2

    azim_vec = np.arange(len(t_binary))*azim_shift
    azim_vec -= azim_vec[stop_idx] +60

    return [elev_vec, azim_vec]


#----------------------------------------------------------------------------
def get_binary_data(q, chiA, chiB, omega_ref):

    mA = q/(1.+q)
    mB = 1./(1.+q)

    nr_sur = NRSur7dq2.NRSurrogate7dq2()

    # If omega_ref is not given, set f_ref to None, and t_ref to -100
    f_ref = None if omega_ref is None else omega_ref/np.pi
    t_ref = -100 if omega_ref is None else None

    # get NRSur dynamics
    quat_nrsur, orbphase_nrsur, _, _ \
        = nr_sur.get_dynamics(q, chiA, chiB, omega_ref=omega_ref, t_ref=t_ref, \
        allow_extrapolation=True)

    t_binary = get_uniform_in_orbits_times(nr_sur.tds, orbphase_nrsur, \
        PTS_PER_ORBIT)

    if LOW_DEF:
        t_binary = t_binary[t_binary > -3000]

    # If FREEZE_TIME is not in t_binary, add it
    if np.min(np.abs(t_binary - FREEZE_TIME)) > 0.1:
        t_binary = np.sort(np.append(t_binary, FREEZE_TIME))

    # If t=0 is not in t_binary, add it
    if np.min(np.abs(t_binary - 0)) > 0.1:
        t_binary = np.sort(np.append(t_binary, 0))

    # interpolate dynamics on to t_binary
    quat_nrsur = np.array([spline_interp(t_binary, nr_sur.tds, tmp) \
        for tmp in quat_nrsur])
    orbphase_nrsur = spline_interp(t_binary, nr_sur.tds, orbphase_nrsur)

    omega_nrsur = get_omegaOrb_from_sparse_data(t_binary, orbphase_nrsur)

    h_nrsur, chiA_nrsur, chiB_nrsur = nr_sur(q, chiA, chiB, \
        f_ref=f_ref, t_ref=t_ref, return_spins=True, \
        allow_extrapolation=True, t=t_binary)

    LHat = surfinBH._utils.lHat_from_quat(quat_nrsur).T
    separation = get_separation_from_omega(omega_nrsur, mA, mB, chiA_nrsur, \
        chiB_nrsur, LHat)

    # Newtonian
    LMag = q/(1.+q)**2 * omega_nrsur**(-1./3)
    L = LHat*LMag[:, None]

    # Get component trajectories
    BhA_traj = get_trajectory(separation * mB, quat_nrsur, orbphase_nrsur, 'A')
    BhB_traj = get_trajectory(separation * mA, quat_nrsur, orbphase_nrsur, 'B')

    return t_binary, chiA_nrsur, chiB_nrsur, L, h_nrsur, BhA_traj, \
        BhB_traj, separation



#----------------------------------------------------------------------------
def make_zero_if_small(x):
    if abs(x) < 1e-3:
        return 0
    else:
        return x

#----------------------------------------------------------------------------
def update_lines(num, lines, hist_frames, t, t_binary, dataLines_binary, \
        dataLines_remnant, properties_text, freeze_text, timestep_text, \
        time_text, max_range, BhA_traj, BhB_traj, BhC_traj, L, h_nrsur, \
        sph_gridX, gridX, sph_gridY, gridY, sph_gridZ, gridZ, \
        q, mA, mB, chiA_nrsur, chiB_nrsur, mf, chif, vf, \
        waveform_end_time, freeze_idx, draw_full_trajectory, ax, vmin, vmax, \
        linthresh, camera_traj, height_map, project_on_all_planes, \
        wave_time_series):
    """ The function that goes into animation
    """
    current_time = t[num]
    time_text.set_text('$t=%.1f\,M$'%current_time)

    if num == freeze_idx - 1:
        # Add text about freezing before freezing
        freeze_text.set_text('Freezing video')
    if num == freeze_idx + 1:
        # Clear text about freezing after freezing
        freeze_text.set_text('')

    if current_time < waveform_end_time:
        # Plot the waveform on the back planes
        if project_on_all_planes:
            hplusX = get_waveform_on_grid(t, num-1, h_nrsur, sph_gridX)
            hplusY = get_waveform_on_grid(t, num-1, h_nrsur, sph_gridY)
        hplusZ = get_waveform_on_grid(t, num-1, h_nrsur, sph_gridZ)
        ax.collections = []     # It becomes very slow without this
        norm=colors.SymLogNorm(linthresh=linthresh, linscale=1, \
            vmin=vmin, vmax=vmax)
        if project_on_all_planes:
            ax.contourf(hplusX, gridX[0], gridX[1], zdir='x', \
                offset=-max_range, cmap=cm.coolwarm, \
                zorder=zorder_dict['contourf'], vmin=vmin, vmax=vmax, \
                norm=norm)
            ax.contourf(gridY[0], hplusY, gridY[1], zdir='y', \
                offset=max_range, cmap=cm.coolwarm, \
                zorder=zorder_dict['contourf'], vmin=vmin, vmax=vmax, \
                norm=norm)
        if height_map:
            # TODO Magic constants, fix the color map.
            ax.plot_surface(gridZ[0], gridZ[1], - max_range + 3.*hplusZ/vmax, \
                cmap=cm.coolwarm, zorder=zorder_dict['contourf'])
        else:
            ax.contourf(gridZ[0], gridZ[1], hplusZ, zdir='z', \
                offset=-max_range, cmap=cm.coolwarm, \
                zorder=zorder_dict['contourf'], vmin=vmin, vmax=vmax, norm=norm)
    else:
        timestep_text.set_text('Increased time step to 100M')

    if current_time < 0:        # Show binary until t=0
        if num < 2:
            # Clear remnant stuff
            line = lines[len(dataLines_binary)]
            line.set_data([], [])
            line.set_3d_properties([])
            line = lines[len(dataLines_binary)+1]
            line.reset()
            timestep_text.set_text('')

        properties_text.set_text('$q=%.2f$\n' \
            '$\chi_{A}=[%.2f, %.2f, %.2f]$\n' \
            '$\chi_{B}=[%.2f, %.2f, %.2f]$\n'%(q, \
            make_zero_if_small(chiA_nrsur[num-1][0]), \
            make_zero_if_small(chiA_nrsur[num-1][1]), \
            make_zero_if_small(chiA_nrsur[num-1][2]), \
            make_zero_if_small(chiB_nrsur[num-1][0]), \
            make_zero_if_small(chiB_nrsur[num-1][1]), \
            make_zero_if_small(chiB_nrsur[num-1][2]), \
            ))


        for idx in range(len(dataLines_binary)):

            line = lines[idx]
            data = dataLines_binary[idx]


            if idx < 4:
                if idx < 2:
                    if draw_full_trajectory:
                        start = 0
                    else:
                        start = max(0, num-hist_frames)
                else:
                    start = max(0, num-1)

                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, start:num])
                line.set_3d_properties(data[2, start:num])
            else:
                if idx == 4:
                    line.set_BH_spin_arrow(BhA_traj[:,num-1], mA, \
                        chiA_nrsur[num-1])
                elif idx == 5:
                    line.set_BH_spin_arrow(BhB_traj[:,num-1], mB, \
                        chiB_nrsur[num-1])
                elif idx == 6:
                    line.set_angular_momentum_arrow(L.T[:,num-1])

    else:
        if abs(current_time) < 10:
            # Clear binary stuff
            for idx in range(4):
                line = lines[idx]
                line.set_data([], [])
                line.set_3d_properties([])
            for idx in range(4,7):
                line = lines[idx]
                line.reset()


        properties_text.set_text('$m_f=%.2f\,M$\n' \
            '$\chi_f=[%.2f, %.2f, %.2f]$\n' \
            '$v_f = [%.2f, %.2f, %.2f] \\times 10^{-3} c$'%(mf, \
            chif[0], chif[1], chif[2], vf[0]*1e3, vf[1]*1e3, vf[2]*1e3))

        for idx in range(len(dataLines_remnant)):
            line = lines[len(dataLines_binary)+idx]
            data = dataLines_remnant[idx]

            if idx == 0:
                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, num-1:num])
                line.set_3d_properties(data[2, num-1:num])
            else:
                Bh_loc = BhC_traj[:,num-1]
                chi_vec = chif
                mass = mf
                line.set_BH_spin_arrow(Bh_loc, mass, chi_vec)


    if current_time < -500:
        if camera_traj is not None:
            ax.view_init(elev=camera_traj[0][num], azim=camera_traj[1][num])

    # Plot waveform time series
    if wave_time_series:
        h_viewpoint = get_waveform_timeseries(h_nrsur, ax.azim, ax.elev)
        for idx in range(3):
            line = lines[len(dataLines_binary)+len(dataLines_remnant)+idx]
            if idx == 0:
                line.set_data(t_binary, np.real(h_viewpoint))
            elif idx == 1:
                line.set_data(t_binary, np.imag(h_viewpoint))
            else:
                line.set_xdata(current_time)


    return lines


#----------------------------------------------------------------------------
def BBH_scattering(fig, q, chiA, chiB, omega_ref=None, \
        draw_full_trajectory=False, project_on_all_planes=False, \
        height_map=False, wave_time_series=False, \
        auto_rotate_camera=False, save_file=None, still_time=None,  \
        rescale_fig_for_widgets=False):

    chiA = np.array(chiA)
    chiB = np.array(chiB)
    t_binary, chiA_nrsur, chiB_nrsur, L, h_nrsur, BhA_traj, \
        BhB_traj, separation = get_binary_data(q, chiA, chiB, omega_ref)

    max_range = np.nanmax(separation)

    mA = q/(1.+q)
    mB = 1./(1.+q)

    # Get mesh grid on bottom plane to generate waveform
    sph_gridX, gridX, sph_gridY, gridY, sph_gridZ, gridZ \
        = get_grids_on_planes(11, max_range)

    # evaluate remnant fit
    fit_name = 'surfinBH7dq2'
    fit = surfinBH.LoadFits(fit_name)

    # If omega_ref is None, will assume the spins are given in the
    # coorbital frame at t=-100M
    mf, chif, vf, mf_err, chif_err, vf_err \
        = fit.all(q, chiA, chiB, omega0=omega_ref)

    #print(np.linalg.norm(chif))
    #print(np.linalg.norm(vf))
    #print(np.linalg.norm(vf) * 3 * 10**5)

    # Will stop plotting waveform after this time
    # long enough for waveform pattern to disappear, taking into account
    # propagation delay
    waveform_end_time = 50 + 2*max_range

    # common time array: After waveform_end_time, each step is 100M
    t = np.append(t_binary[t_binary<waveform_end_time], \
        np.arange(waveform_end_time, 10000+waveform_end_time, 100))

    # assume merger is at origin
    BhC_traj = np.array([tmp*t for tmp in vf])

    # Attaching 3D axis to the figure
    ax = axes3d.Axes3D(fig)

    if wave_time_series:
        l, b, w, h = ax.get_position().bounds
        if rescale_fig_for_widgets:
            ax.set_position([l, b + 0.25*h, w*0.7, 0.75*h ])
            # axes to plot waveform time series
            hax = fig.add_axes([0.135, 0.08, 0.5, 0.17])
            hax.clear()
        else:
            ax.set_position([l, b + 0.25*h, w, 0.75*h ])
            # axes to plot waveform time series
            hax = fig.add_axes([0.135, 0.08, 0.83, 0.17])


        # estimate maximum of waveform for scale of timeseries
        hmax_est = np.max(np.abs(get_waveform_timeseries(h_nrsur, 0, 90)))
        hax.set_ylim([ -hmax_est, hmax_est ])

    if auto_rotate_camera:
        camera_traj = get_camera_trajectory(t_binary)
    else:
        camera_traj = None

    markersize_BhA = get_marker_size(mA, chiA)
    markersize_BhB = get_marker_size(mB, chiB)
    markersize_BhC = get_marker_size(mf, chif)

    if LOW_DEF:
        time_fontsize = 5
        properties_fontsize = 5
        properties_text_yloc = 0.75
        freeze_fontsize = 7
        timestep_fontsize = 6
        label_fontsize = 5
        ticks_fontsize = 5
        title_fontsize = 7
        ticks_pad = -5
        label_pad = -11
    else:
        time_fontsize = 12
        properties_fontsize = 10
        properties_text_yloc = 0.8
        freeze_fontsize = 14
        timestep_fontsize = 12
        label_fontsize = 10
        ticks_fontsize = 10
        title_fontsize = 14
        ticks_pad = 0
        label_pad = 0

    time_text = ax.text2D(0.03, 0.05, '', transform=ax.transAxes, \
        fontsize=time_fontsize, zorder=zorder_dict['info_text'])
    properties_text = ax.text2D(0.05, properties_text_yloc, '', \
        transform=ax.transAxes, fontsize=properties_fontsize, \
        zorder=zorder_dict['info_text'])
    freeze_text = ax.text2D(0.6, 0.7, '', transform=ax.transAxes, \
        fontsize=freeze_fontsize, color=colors_dict['info'], \
        zorder=zorder_dict['notice_text'])
    timestep_text = ax.text2D(0.45, 0.7, '', transform=ax.transAxes, \
        fontsize=timestep_fontsize, color=colors_dict['info'], \
        zorder=zorder_dict['notice_text'])


    # NOTE: Can't pass empty arrays into 3d version of plot()
    dataLines_binary = [BhA_traj, BhB_traj, BhA_traj, BhB_traj, 1, 1, 1]

    # get wavefrom at viewpoint
    h_viewpoint = get_waveform_timeseries(h_nrsur, ax.azim, ax.elev)

    if LOW_DEF:
        arrow_mutation_scale = 10
    else:
        arrow_mutation_scale = 20

    marker_alpha = 0.9
    traj_alpha = 0.8
    lines = [\
        # These two are for plotting component tracjectories
        ax.plot(BhA_traj[0,0:1]-1e10, BhA_traj[1,0:1], BhA_traj[2,0:1], \
            color=colors_dict['BhA_traj'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ax.plot(BhB_traj[0,0:1]-1e10, BhB_traj[1,0:1], BhB_traj[2,0:1], \
            color=colors_dict['BhB_traj'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \

        # These two are for plotting component BHs
        ax.plot(BhA_traj[0,0:1]-1e10, BhA_traj[1,0:1], BhA_traj[2,0:1], \
            marker='o', markersize=markersize_BhA, markerfacecolor='k', \
            markeredgewidth=0, alpha=marker_alpha, \
            zorder=zorder_dict['Bh'])[0], \
        ax.plot(BhB_traj[0,0:1]-1e10, BhB_traj[1,0:1], BhB_traj[2,0:1], \
            marker='o', markersize=markersize_BhB, markerfacecolor='k',
            markeredgewidth=0, alpha=marker_alpha, \
            zorder=zorder_dict['Bh'])[0], \

        # These two are for plotting component BH spins
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['BhA_spin'], \
            zorder=zorder_dict['spin'])), \
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['BhB_spin'], \
            zorder=zorder_dict['spin'])), \

        # This is for plotting angular momentum direction
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['L'], \
            zorder=zorder_dict['L'])), \

        # This is for plotting remnant BH
        ax.plot(BhC_traj[0,0:1]-1e10, BhC_traj[1,0:1], BhC_traj[2,0:1], \
            marker='o', markersize=markersize_BhC, markerfacecolor='k', \
            markeredgewidth=0, alpha=marker_alpha, \
            zorder=zorder_dict['Bh'])[0], \
        # This is for plotting remnant spin
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['BhC_spin'], \
            zorder=zorder_dict['spin'])), \
        ]

    if wave_time_series:
        lines += [ \
            # These two is for plotting the waveform time series
            hax.plot(t_binary, np.real(h_viewpoint), label='$h_+$', \
                color=colors_dict['h+'], lw=1.2)[0], \
            hax.plot(t_binary, np.imag(h_viewpoint), label='$h_{\\times}$', \
                color=colors_dict['hx'], lw=1.2)[0], \

            # This is for plotting the slider along the waveform time series
            hax.axvline(x=t_binary[0]), \
            ]

        hax.legend(loc='upper left', ncol=2)
        hax.set_xlabel('$t\,(M)$', fontsize=label_fontsize)
        hax.set_ylabel('$h\,r/M$', fontsize=label_fontsize)
        hax.tick_params(axis='x', which='major', labelsize=ticks_fontsize)
        hax.tick_params(axis='y', which='major', labelsize=ticks_fontsize)

    dataLines_remnant = [BhC_traj, 1]


    # Setting the axes properties

    # This seems to set the actual limits to max_range
    ax.set_xlim3d([-max_range*0.96, max_range*0.96])
    ax.set_ylim3d([-max_range*0.96, max_range*0.96])
    ax.set_zlim3d([-max_range*0.96, max_range*0.96])

    ax.set_xlabel('$x\,(M)$', fontsize=label_fontsize)
    ax.set_ylabel('$y\,(M)$', fontsize=label_fontsize)
    ax.set_zlabel('$z\,(M)$', fontsize=label_fontsize)

    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')

    ax.set_facecolor('white')

    #ax.xaxis._axinfo['tick']['inward_factor'] = 0
    #ax.yaxis._axinfo['tick']['inward_factor'] = 0
    #ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0

    ax.tick_params(axis='x', which='major', pad=ticks_pad, \
        labelsize=ticks_fontsize)
    ax.tick_params(axis='y', which='major', pad=ticks_pad, \
        labelsize=ticks_fontsize)
    ax.tick_params(axis='z', which='major', pad=ticks_pad, \
        labelsize=ticks_fontsize)

    ax.xaxis.labelpad = label_pad
    ax.yaxis.labelpad = label_pad
    ax.zaxis.labelpad = label_pad -3

    ax.set_title('NRSur7dq2 + %s'%fit_name, fontsize=time_fontsize, \
        x=0.74, y=0.99)

    # number of frames to include in orbit trace
    hist_frames = int(3./4*(PTS_PER_ORBIT))

    # Will freeze video at this index
    freeze_idx = np.argmin(np.abs(t - FREEZE_TIME))


    # color range for contourf
    # Get linthresh from first index. With SymLogNorm, whenever the
    # value is less than linthresh, the color scale is linear. Else log.
    linthresh = np.max(np.abs(get_waveform_on_grid(t, 0, h_nrsur, sph_gridZ)))
    # Get vmax from waveform at peak.  Add in propagation delay
    zero_idx = np.argmin(np.abs(t-max_range))
    vmax = np.max(get_waveform_on_grid(t, zero_idx, h_nrsur, \
                                       sph_gridZ))
    # Symmetric about 0
    vmin = -vmax

    #NOTE: There is a glitch if I don't skip the first index
    frames = range(1, len(t))
    # Repeat freeze_idx 75 times, this is a hacky way to freeze the video
    frames = np.sort(np.append(frames, [freeze_idx]*75))

    fargs = (lines, hist_frames, t, t_binary, dataLines_binary, \
            dataLines_remnant, properties_text, freeze_text, timestep_text, \
            time_text, max_range, BhA_traj, BhB_traj, BhC_traj, L, h_nrsur, \
            sph_gridX, gridX, sph_gridY, gridY, sph_gridZ, gridZ, \
            q, mA, mB, chiA_nrsur, chiB_nrsur, mf, chif, vf, \
            waveform_end_time, freeze_idx, draw_full_trajectory, ax, \
            vmin, vmax, linthresh, camera_traj, height_map, \
            project_on_all_planes, wave_time_series)

    # save still and exit
    if still_time is not None:
        time_tag = '%d'%(abs(still_time))
        if still_time < 0:
            time_tag = 'm%s'%time_tag
        update_lines(np.argmin(np.abs(t-still_time)), *fargs)
        still_fnametag = '%s_%s'%(save_file.split('.')[0], time_tag)
        P.savefig('%s.png'%still_fnametag, bbox_inches='tight')
        P.savefig('%s.pdf'%still_fnametag, bbox_inches='tight')
        exit()

    line_ani = animation.FuncAnimation(fig, update_lines, frames, \
        fargs=fargs, \
        interval=50, blit=False, repeat=True, repeat_delay=5e3)

    return line_ani


#############################    main    ##################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--q', type=float, required=True,
        help='Mass ratio.')
    parser.add_argument('--chiA', type=float, required=True, nargs=3,
        help='Dimensionless spin of BhA at omega_ref. List of size 3.')
    parser.add_argument('--chiB', type=float, required=True, nargs=3,
        help='Dimensionless spin of BhB at omega_ref. List of size 3.')
    parser.add_argument('--omega_ref', type=float, default=None,
        help='Orbital frequency at which the spins are specified. In units ' \
            'of rad/M. Currently, >= 0.018. If not specified, assumes the ' \
            'spins are specified at t=-100M from the peak of the waveform.')
    parser.add_argument('--save_file', type=str, default=None,
        help='File to save animation to. If given, will save animation to ' \
                'this file. Else will show animation. Allowed extensions are ' \
                'mp4 and gif. mp4 has the best quality. We use lower quality ' \
                'for gif to reduce file size.')
    parser.add_argument('--project_on_all_planes', default=False, \
        action='store_true', \
        help='If given, projects the waveform on all three back planes. ' \
        'By default only does the x-y plane at the bottom.')
    parser.add_argument('--height_map', default=False, action='store_true', \
        help='Map h to a height to visualize in the xy plane ' \
        'instead of a contour plot.  Turns off project_on_all_planes.')
    parser.add_argument('--wave_time_series', default=False, \
        action='store_true', \
        help='Plots an interactive waveform time series at the bottom. The ' \
        'waveform changes based on viewing angle. This disables ' \
        'pause-on-click.')
    parser.add_argument('--auto_rotate_camera', default=False, \
        action='store_true', \
        help='Auto rotates camera viewing angle. This turns off ' \
        'project_on_all_planes. Particularly nice with wave_time_series set.')
    parser.add_argument('--draw_full_trajectory', default=False, \
        action='store_true', \
        help='If given, draws the entire trajectories of the components. ' \
        'Else only retains the last 3/4th of an orbit.')
    parser.add_argument('--still_time', default=None, type=float, \
        help='If given, saves a plot of the movie at this time and exits.')

    args = parser.parse_args()
    if args.height_map or args.auto_rotate_camera:
        args.project_on_all_planes=False


    if LOW_DEF:
        fig = P.figure(figsize=(2.3,2))
    else:
        if args.wave_time_series:
            fig = P.figure(figsize=(5,5.5))
        else:
            fig = P.figure(figsize=(5,4))


    line_ani = BBH_scattering(fig, args.q, args.chiA, args.chiB,
        omega_ref = args.omega_ref,
        draw_full_trajectory = args.draw_full_trajectory,
        height_map = args.height_map,
        project_on_all_planes = args.project_on_all_planes,
        wave_time_series = args.wave_time_series,
        auto_rotate_camera = args.auto_rotate_camera,
        save_file = args.save_file,
        still_time = args.still_time)

    if args.save_file is not None:
        # Set up formatting for the movie files

        extension = args.save_file.split('.')[-1]
        if extension == 'mp4':
            # Might need: conda install -c conda-forge ffmpeg
            Writer = animation.writers['ffmpeg']
        elif extension == 'gif':
            # Might need: brew install imagemagick
            Writer = animation.writers['imagemagick']
        else:
            raise Exception('Invalid extension')

        metadata = {
            'artist' : 'Vijay Varma',
            'genre' : 'Physics',
            'subject' : 'Animation of binary black hole scattering process.',
            'copyright' : surfinBH.__copyright__,
            }
        writer = Writer(fps=15, metadata=metadata)
        if LOW_DEF or extension == 'gif':
            line_ani.save(args.save_file, writer=writer)
        else:
            line_ani.save(args.save_file, writer=writer, dpi=150)

    else:
        # Pause settings
        pause = False
        def onClick(event):
            global pause
            if pause:
                line_ani.event_source.start()
                pause = False
            else:
                line_ani.event_source.stop()
                pause = True

        # The waveform does not update when you rotate when paused, so
        # disable pausing if plotting waveform time series
        if not args.wave_time_series:
            fig.canvas.mpl_connect('button_press_event', onClick)
        P.show()