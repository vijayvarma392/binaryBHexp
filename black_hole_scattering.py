desc = """Animations of binary black hole scattering.
Generates an animation of a binary black hole merger and the final remnant.

Example usage:
python black_hole_scattering.py --q 2 --omega_ref 1.8e-2 --chiA 0.2 0.7 -0.1 --chiB 0.2 0.6 0.1

Note: Time values displayed in the plot are non-uniform and non-linear:
During the inspiral there are 30 frames per orbit.
After the merger each frame corresponds to a time step of 100M.
"""

import numpy as np
import matplotlib.pyplot as P
import argparse
import time

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

import surfinBH
import NRSur7dq2

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch

P.style.use('seaborn')

colors_dict = {
        'BhA_traj': 'indianred',
        'BhB_traj': 'rebeccapurple',
        'BhA_spin': 'goldenrod',
        'BhB_spin': 'steelblue',
        'BhC_spin': 'forestgreen',
        }

# number of frames per orbit
pts_per_orbit = 30

# Time at which to freeze video for 5 seconds
freeze_time = -100

class Arrow3D(FancyArrowPatch):
    def __init__(self, vecs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = vecs

    def set_BH_spin_arrow(self, Bh_loc, mass, chi_vec, scale_factor=25):
        """ The length of the arrow is proportinal to the Kerr parameter
        a of the BH.
        """
        x, y, z =  Bh_loc
        u, v, w =  chi_vec*mass
        xs = [x, x+u*scale_factor]
        ys = [y, y+v*scale_factor]
        zs = [z, z+w*scale_factor]
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
        if np.min(newX) < np.min(oldX) or np.max(newX) > np.max(oldX):
            print np.min(newX), np.min(oldX), np.max(newX), np.max(oldX)
            print np.min(newX) < np.min(oldX)
            print np.max(newX) > np.max(oldX)
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
def get_uniform_in_orbits_times(t, phi_orb, pts_per_orbit):
    """
    returns sparse time array such that there are pts_per_orbit points
    in each orbit.
    """
    # get numer of orbits
    n_orbits = int(abs((phi_orb[-1] - phi_orb[0])/(2*np.pi)))

    # get sparse times such that there are pts_per_orbit points in each orbit
    n_pts = int(n_orbits*pts_per_orbit)
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
    chimag = np.sqrt(np.sum(chi**2))
    rplus = mass + mass*np.sqrt(1 - chimag**2)
    return rplus * 30

#----------------------------------------------------------------------------
def get_separation_from_omega(omega, mA, mB, chiA, chiB, LHat):
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
    gamma = x * ( 1 + x * (1. - 1./3 *eta) \
            + x**(3./2) * (5./3 * SL + deltaM * SigmaL ) \
            + x**2 * (1 - 65./12 *eta) \
            + x**(5./2) * ( (10./3 + 8./9 * eta)*SL + 2* deltaM * SigmaL) \
            + x**3 * (1. + (-2203./2520 -41./192 * np.pi**2)*eta \
                + 229./36 * eta**2 + 1./81 * eta**3) \
            + x**(7./2) * ( (5 - 127./12 *eta - 6 * eta**2)*SL + \
                deltaM * SigmaL * (3 - 61./6 *eta - 8./3 * eta**2) ) \
            )
    r = 1/gamma

    # To this add the 2PN spin-spin term from Eq.(4.13) of
    # https://arxiv.org/pdf/gr-qc/9506022.pdf
    r += omega**(-2./3) * (-1./2 * eta * chiAB) * omega**(4./3)

    return r


#----------------------------------------------------------------------------
def make_zero_if_small(x):
    if abs(x) < 1e-3:
        return 0
    else:
        return x

#----------------------------------------------------------------------------
def update_lines(num, lines, hist_frames, t, dataLines_binary, \
        dataLines_remnant, time_text, properties_text, freeze_text, \
        BhA_traj, BhB_traj, BhC_traj, \
        q, mA, mB, chiA_nrsur, chiB_nrsur, mf, chif, vf, zero_idx, freeze_idx):
    """ The function that goes into animation
    """
    current_time = t[num]
    time_text.set_text('$t=%.1f\,M$'%current_time)


    if num == freeze_idx:
        # Add text about freezing before freezing
        freeze_text.set_text('Freezing video')
    if num == freeze_idx+1:
        time.sleep(5)
        # Clear text about freezing after freezing
        freeze_text.set_text('')


    if current_time < 0:

        if num < 2:
            # Clear remnant stuff
            line = lines[6]
            line.set_data([], [])
            line.set_3d_properties([])
            line = lines[7]
            line.reset()

        for idx in range(len(dataLines_binary)):
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

            line = lines[idx]
            data = dataLines_binary[idx]


            if idx < 4:
                if idx < 2:
                    start = max(0, num-hist_frames)
                else:
                    start = max(0, num-1)

                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, start:num])
                line.set_3d_properties(data[2, start:num])
            else:
                if idx == 4:
                    Bh_loc = BhA_traj[:,num-1]
                    chi_vec = chiA_nrsur[num-1]
                    mass = mA
                elif idx == 5:
                    Bh_loc = BhB_traj[:,num-1]
                    chi_vec = chiB_nrsur[num-1]
                    mass = mB

                line.set_BH_spin_arrow(Bh_loc, mass, chi_vec)
    else:
        num = num - zero_idx + 1    # Ignore first index to avoid glitch
        if num < 2:
            # Clear binary stuff
            for idx in range(4):
                line = lines[idx]
                line.set_data([], [])
                line.set_3d_properties([])
            for idx in range(4,6):
                line = lines[idx]
                line.reset()


            properties_text.set_text('$m_f=%.2f\,M$\n' \
                '$\chi_f=[%.2f, %.2f, %.2f]$\n' \
                '$v_f = [%.2f, %.2f, %.2f] \\times 10^{-3} c$'%(mf, \
                chif[0], chif[1], chif[2], vf[0]*1e3, vf[1]*1e3, vf[2]*1e3))

        for idx in range(len(dataLines_remnant)):
            line = lines[6+idx]
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

    return lines


#----------------------------------------------------------------------------
def BBH_scattering(q, chiA, chiB, omega_ref, return_fig=False):

    chiA = np.array(chiA)
    chiB = np.array(chiB)

    # evaluate remnant fit
    fit_name = 'surfinBH7dq2'
    fit = surfinBH.LoadFits(fit_name)

    # If omega_ref is None, will assume the spins are given in the
    # coorbital frame at t=-100M
    mf, chif, vf, mf_err, chif_err, vf_err \
        = fit.all(q, chiA, chiB, omega0=omega_ref)

    print np.linalg.norm(chif)
    print np.linalg.norm(vf)
    print np.linalg.norm(vf) * 3 * 10**5

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
        pts_per_orbit)

    # If freeze_time is not in t_binary, add it
    if np.min(np.abs(t_binary - freeze_time)) > 0.1:
        t_binary = np.sort(np.append(t_binary, freeze_time))

    # interpolate dynamics on to t_binary
    quat_nrsur = np.array([spline_interp(t_binary, nr_sur.tds, tmp) \
        for tmp in quat_nrsur])
    orbphase_nrsur = spline_interp(t_binary, nr_sur.tds, orbphase_nrsur)

    omega_nrsur = get_omegaOrb_from_sparse_data(t_binary, orbphase_nrsur)

    h_nrsur, chiA_nrsur, chiB_nrsur = nr_sur(q, chiA, chiB, \
        f_ref=f_ref, t_ref=t_ref, return_spins=True, \
        allow_extrapolation=True, LMax=2, t=t_binary)

    LHat = surfinBH._utils.lHat_from_quat(quat_nrsur).T
    separation = get_separation_from_omega(omega_nrsur, mA, mB, chiA_nrsur, \
        chiB_nrsur, LHat)

    # Get component trajectories
    BhA_traj = get_trajectory(separation * mB, quat_nrsur, orbphase_nrsur, 'A')
    BhB_traj = get_trajectory(separation * mA, quat_nrsur, orbphase_nrsur, 'B')

    # time array for remnant
    t_remnant = np.arange(0, 10000, 100)

    # assume merger is at origin
    BhC_traj = np.array([tmp*t_remnant for tmp in vf])

    # Attaching 3D axis to the figure
    fig = P.figure(figsize=(5,4))
    ax = axes3d.Axes3D(fig)

    markersize_BhA = get_marker_size(mA, chiA)
    markersize_BhB = get_marker_size(mB, chiB)
    markersize_BhC = get_marker_size(mf, chif)

    time_text = ax.text2D(0.03, 0.05, '', transform=ax.transAxes, fontsize=12)
    properties_text = ax.text2D(0.05, 0.8, '', transform=ax.transAxes, \
        fontsize=10)
    freeze_text = ax.text2D(0.6, 0.7, '', transform=ax.transAxes, fontsize=14,
            color='tomato')
    

    # NOTE: Can't pass empty arrays into 3d version of plot()
    dataLines_binary = [BhA_traj, BhB_traj, BhA_traj, BhB_traj, 1, 1]

    marker_alpha = 0.9
    traj_alpha = 0.8
    lines = [\
        # These two are for plotting component tracjectories
        ax.plot(BhA_traj[0,0:1]-1e10, BhA_traj[1,0:1], BhA_traj[2,0:1], \
            color=colors_dict['BhA_traj'], lw=2, alpha=traj_alpha)[0], \
        ax.plot(BhB_traj[0,0:1]-1e10, BhB_traj[1,0:1], BhB_traj[2,0:1], \
            color=colors_dict['BhB_traj'], lw=2, alpha=traj_alpha)[0], \

        # These two are for plotting component BHs
        ax.plot(BhA_traj[0,0:1]-1e10, BhA_traj[1,0:1], BhA_traj[2,0:1], \
            marker='o', markersize=markersize_BhA, markerfacecolor='k', \
            markeredgewidth=0, alpha=marker_alpha)[0], \
        ax.plot(BhB_traj[0,0:1]-1e10, BhB_traj[1,0:1], BhB_traj[2,0:1], \
            marker='o', markersize=markersize_BhB, markerfacecolor='k',
            markeredgewidth=0, alpha=marker_alpha)[0], \

        # These two are for plotting component BH spins
        ax.add_artist(Arrow3D(None, mutation_scale=20, lw=3, arrowstyle="-|>", \
            color=colors_dict['BhA_spin'])), \
        ax.add_artist(Arrow3D(None, mutation_scale=20, lw=3, arrowstyle="-|>", \
            color=colors_dict['BhB_spin'])), \

        # This is for plotting remnant BH
        ax.plot(BhC_traj[0,0:1]-1e10, BhC_traj[1,0:1], BhC_traj[2,0:1], \
            marker='o', markersize=markersize_BhC, markerfacecolor='k', \
            markeredgewidth=0, alpha=marker_alpha)[0], \
        # This is for plotting remnant spin
        ax.add_artist(Arrow3D(None, mutation_scale=20, lw=3, arrowstyle="-|>", \
            color=colors_dict['BhC_spin'])), \

        ]

    dataLines_remnant = [BhC_traj, 1]

    max_range = np.nanmax(separation)

    # Setting the axes properties
    ax.set_xlim3d([-max_range, max_range])
    ax.set_ylim3d([-max_range, max_range])
    ax.set_zlim3d([-max_range, max_range])

    ax.set_xlabel('$x\,(M)$', fontsize=10)
    ax.set_ylabel('$y\,(M)$', fontsize=10)
    ax.set_zlabel('$z\,(M)$', fontsize=10)

    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')

    ax.set_facecolor('white')

    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['inward_factor'] = 0

    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4

    ax.tick_params(axis='x', which='major', pad=0)
    ax.tick_params(axis='y', which='major', pad=0)
    ax.tick_params(axis='z', which='major', pad=0)

    ax.xaxis.labelpad = 0
    ax.yaxis.labelpad = 0
    ax.zaxis.labelpad = -3

    ax.set_title('NRSur7dq2 + %s'%fit_name, fontsize=14, x=0.74, y=0.99)

    # number of frames to include in orbit trace
    hist_frames = int(3./4*(pts_per_orbit))

    # common time array
    t = np.append(t_binary[t_binary<0], t_remnant)

    # Will switch to remant after this index
    zero_idx = np.argmin(np.abs(t))

    # Will freeze for 5 seconds at this index
    freeze_idx = np.argmin(np.abs(t - freeze_time))


    #NOTE: There is a glitch if I don't skip the first index
    line_ani = animation.FuncAnimation(fig, update_lines, range(1, len(t)), \
        fargs=(lines, hist_frames, t, dataLines_binary, dataLines_remnant, \
            time_text, properties_text, freeze_text, \
            BhA_traj, BhB_traj, BhC_traj, \
            q, mA, mB, chiA_nrsur, chiB_nrsur, mf, chif, vf, zero_idx, \
            freeze_idx), \
        interval=50, blit=False, repeat=True, repeat_delay=5e3)

    if return_fig:
        return line_ani, fig
    else:
        return line_ani


#############################    main    ##################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=desc,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--q', type=float, required=True,
        help='Mass ratio.')
    parser.add_argument('--chiA', type=float, required=True, nargs=3,
        help='Spin of BhA at omega_ref. Array of size 3.')
    parser.add_argument('--chiB', type=float, required=True, nargs=3,
        help='Spin of BhB at omega_ref. Array of size 3.')
    parser.add_argument('--omega_ref', type=float, default=None,
        help='Starting orbital frequency at which the spins are specified. ' \
            'Currently, > 0.018. If not specified, assumes the spins are ' \
            'specified at t=-100M from the peak of the waveform.')
    parser.add_argument('--save_file', type=str, default=None,
        help='File (without extension) to save animation to. If given, will' \
            ' save animation to this file. Else will show animation.')
    parser.add_argument('--save_format', type=str, default='mp4',
        help='Format for video file, mp4 or gif.')

    args = parser.parse_args()
    line_ani, fig = BBH_scattering(args.q, args.chiA, args.chiB, \
        args.omega_ref, return_fig=True)

    if args.save_file is not None:
        # Set up formatting for the movie files

        if args.save_format == 'mp4':
            # Might need: conda install -c conda-forge ffmpeg
            Writer = animation.writers['ffmpeg']
        elif args.save_format == 'gif':
            # Might need: brew install imagemagick
            Writer = animation.writers['imagemagick']
        else:
            raise Exception('Invalid save_format')

        metadata = {
            'artist' : 'Vijay Varma',
            'title' :'surfinBH animation',
            'genre' : 'Physics',
            'subject' : 'Animation of binary black hole scattering process.',
            'copyright' : surfinBH.__copyright__,
            }
        writer = Writer(fps=15, metadata=metadata)
        line_ani.save('%s.%s'%(args.save_file, args.save_format), writer=writer)

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
        fig.canvas.mpl_connect('button_press_event', onClick)
        P.show()
