desc = """Animations of binary black hole scattering.
Generates an animation of a binary black hole merger and the final remnant.

Example usage:
python animate_surfinBH.py --q 2 --omega0 1.8e-2 --chiA0 0.2 0.7 -0.1 --chiB0 0.2 0.6 0.1

Note: Time values displayed in the plot are non-uniform and non-linear:
During the inspiral there are 30 frames per orbit.
After the merger each frame corresponds to a time step of 100M.
"""

import numpy as np
import matplotlib.pyplot as P
import argparse

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

import surfinBH
import NRSur7dq2
#from gwtools import rotations

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch

P.style.use('seaborn')

colors_dict ={
        'BhA_traj': 'indianred',
        'BhB_traj': 'rebeccapurple',
        'BhA_spin': 'goldenrod',
        'BhB_spin': 'steelblue',
        'BhC_spin': 'forestgreen',
        }


class Arrow3D(FancyArrowPatch):
    def __init__(self, vecs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = vecs

    def set_BH_spin_arrow(self, Bh_loc, chi_vec, scale_factor=12):
        x, y, z =  Bh_loc
        u, v, w =  chi_vec
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
    #FIXME COM should not be equidistant from both BHs

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
def get_separation_from_omega(omega):
    """ Let's do zeroth order approx: Keppler's law """
    separation = omega**(-2./3)
    return separation


#----------------------------------------------------------------------------
def get_quivers(Bh_loc, chi_vec, scale_factor=10):
    """ Gets quivers for spins on a BH
    """
    X, Y, Z =  Bh_loc
    u, v, w =  chi_vec
    segments = (X, Y, Z, X+u*scale_factor, Y+v*scale_factor, Z+w*scale_factor)
    segments = np.array(segments).reshape(6,-1)
    return [[[x, y, z], [u, v, w]] for x, y, z, u, v, w in zip(*list(segments))]

#----------------------------------------------------------------------------
def update_lines(num, lines, hist_frames, t, dataLines_binary, \
        dataLines_remnant, time_text, properties_text, \
        BhA_traj, BhB_traj, BhC_traj, \
        q, chiA_nrsur, chiB_nrsur, mf, chif, vf, zero_idx):
    """ The function that goes into animation
    """
    current_time = t[num]
    time_text.set_text('$t=%.1f\,M$'%current_time)
    if current_time < 0:

        if num == 0:
            # Clear remnant stuff
            line = lines[6]
            line.set_data([], [])
            line.set_3d_properties([])
            line = lines[7]
            line.reset()

        for idx in range(len(dataLines_binary)):
            properties_text.set_text('$q=%.1f$\n' \
                '$\chi_{A}=[%.2f, %.2f, %.2f]$\n' \
                '$\chi_{B}=[%.2f, %.2f, %.2f]$\n'%(q, \
                chiA_nrsur[num-1][0],chiA_nrsur[num-1][1],chiA_nrsur[num-1][2],\
                chiB_nrsur[num-1][0],chiB_nrsur[num-1][1],chiB_nrsur[num-1][2],\
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
                elif idx == 5:
                    Bh_loc = BhB_traj[:,num-1]
                    chi_vec = chiB_nrsur[num-1]

                line.set_BH_spin_arrow(Bh_loc, chi_vec)
    else:
        num = num - zero_idx
        if num == 0:
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
                line.set_BH_spin_arrow(Bh_loc, chi_vec)

    return lines


#----------------------------------------------------------------------------
def BBH_scattering(q, chiA0, chiB0, omega0, return_fig=False):

    chiA0 = np.array(chiA0)
    chiB0 = np.array(chiB0)

    # evaluate remnant fit
    fit_name = 'surfinBH7dq2'
    fit = surfinBH.LoadFits(fit_name)
    mf, chif, vf, mf_err, chif_err, vf_err \
        = fit.all(q, chiA0, chiB0, omega0=omega0)

    mA = q/(1.+q)
    mB = 1./(1.+q)

    nr_sur = NRSur7dq2.NRSurrogate7dq2()

    # get NRSur dynamics
    quat_nrsur, orbphase_nrsur, _, _ \
        = nr_sur.get_dynamics(q, chiA0, chiB0, omega_ref=omega0, \
        allow_extrapolation=True)

    pts_per_orbit = 30
    t_binary = get_uniform_in_orbits_times(nr_sur.tds, orbphase_nrsur, \
        pts_per_orbit)

    # interpolate dynamics on to t_binary
    quat_nrsur = np.array([spline_interp(t_binary, nr_sur.tds, tmp) \
        for tmp in quat_nrsur])
    orbphase_nrsur = spline_interp(t_binary, nr_sur.tds, orbphase_nrsur)

    omega_nrsur = get_omegaOrb_from_sparse_data(t_binary, orbphase_nrsur)

    h_nrsur, chiA_nrsur, chiB_nrsur = nr_sur(q, chiA0, chiB0, \
        f_ref=omega0/np.pi, return_spins=True, allow_extrapolation=True, LMax=2,
        t=t_binary)

    #LHat = rotations.lHat_from_quat(quat_nrsur)

    separation = get_separation_from_omega(omega_nrsur)
    BhA_traj = get_trajectory(separation * mB, quat_nrsur, orbphase_nrsur, 'A')
    BhB_traj = get_trajectory(separation * mA, quat_nrsur, orbphase_nrsur, 'B')


    # time array for remnant
    t_remnant = np.arange(0, 10000, 100)

    # assume merger is at origin
    BhC_traj = np.array([tmp*t_remnant for tmp in vf])

    # Attaching 3D axis to the figure
    fig = P.figure(figsize=(5,4))
    ax = axes3d.Axes3D(fig)

    # FIXME check that this makes sense
    markersize_BhA = mA*50
    markersize_BhB = mB*50
    markersize_BhC = mf*50

    time_text = ax.text2D(0.03, 0.05, '', transform=ax.transAxes, fontsize=12)
    properties_text = ax.text2D(0.05, 0.8, '', transform=ax.transAxes, \
        fontsize=10)

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
    ax.set_xlabel('X')

    ax.set_ylim3d([-max_range, max_range])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-max_range, max_range])
    ax.set_zlabel('Z')

    #ax.xaxis.pane.fill = False
    #ax.yaxis.pane.fill = False
    #ax.zaxis.pane.fill = False

    #ax.grid(False)
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
    ax.zaxis.labelpad = -1



    ax.set_title('NRSur7dq2 + %s'%fit_name, fontsize=14, x=0.75, y=0.99)

    # Creating the Animation object
    hist_frames = 15


    zero_idx = np.argmin(np.abs(t_binary))

    # common time array
    t = np.append(t_binary[:zero_idx], t_remnant)

    line_ani = animation.FuncAnimation(fig, update_lines, len(t), \
        fargs=(lines, hist_frames, t, dataLines_binary, dataLines_remnant, \
            time_text, properties_text, BhA_traj, BhB_traj, BhC_traj, \
            q, chiA_nrsur, chiB_nrsur, mf, chif, vf, zero_idx), \
        interval=50, blit=False, repeat=True, repeat_delay=5e3)

    if return_fig:
        return line_ani, fig
    else:
        return line_ani


#############################    main    ##################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=desc,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--omega0', type=float, required=True,
        help='Starting orbital frequency. Currently, > 0.018.')
    parser.add_argument('--q', type=float, required=True,
        help='Mass ratio.')
    parser.add_argument('--chiA0', type=float, required=True, nargs=3,
        help='Spin of BhA at omega0. Array of size 3.')
    parser.add_argument('--chiB0', type=float, required=True, nargs=3,
        help='Spin of BhB at omega0. Array of size 3.')
    parser.add_argument('--save_file', type=str, default=None,
        help='File (without extension) to save animation to. If given, will' \
            ' save animation to this file. Else will show animation.')
    parser.add_argument('--save_format', type=str, default='mp4',
        help='Format for video file, mp4 or gif.')

    args = parser.parse_args()
    line_ani, fig = BBH_scattering(args.q, args.chiA0, args.chiB0, \
        args.omega0, return_fig=True)

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
