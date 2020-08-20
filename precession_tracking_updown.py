desc = """Animations of black hole spin evolution for precessing binary black holes.

Example usage:
python precession_tracking_updown.py  -r pn_instability/Case_002

Where Case_002 contains Horizons.h5.
"""

import numpy as np
import matplotlib.pyplot as P
import argparse
import h5py

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

from gwtools import harmonics

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
P.style.use('seaborn')

import binaryBHexp

Arrow3D = binaryBHexp.Arrow3D
zorder_dict = binaryBHexp.zorder_dict
colors_dict = binaryBHexp.colors_dict

#----------------------------------------------------------------------------
def update_lines(num, lines, t, time_text, chiA, chiB):
    """ The function that goes into animation
    """
    current_time = t[num]
    time_text.set_text('$t=%.1f\,M$'%current_time)
    for idx in range(len(lines)):
        line = lines[idx]
        if idx == 0:
            line.set_arrow_at_origin(chiA[num-1])
        elif idx == 1:
            line.set_arrow_at_origin(chiB[num-1])
        else:
            if idx == 2:
                data = chiA
            elif idx == 3:
                data = chiB

            line.set_data(data.T[0:2, 0:num])
            line.set_3d_properties(data.T[2, 0:num])

    return lines

#----------------------------------------------------------------------------
def PrecessionTrack(fig, RunDir, save_file=None, still_time=None, \
        rescale_fig_for_widgets=False):

    h5file = h5py.File('%s/Horizons.h5'%RunDir, 'r')
    mA = h5file['AhA.dir/ChristodoulouMass.dat'].value.T[1][0]
    t_orig, chiAx, chiAy, chiAz = h5file['AhA.dir/chiInertial.dat'].value.T
    chiA = np.array([chiAx, chiAy, chiAz])

    mB = h5file['AhB.dir/ChristodoulouMass.dat'].value.T[1][0]
    t_orig, chiBx, chiBy, chiBz = h5file['AhB.dir/chiInertial.dat'].value.T
    chiB = np.array([chiBx, chiBy, chiBz])

    # Get time steps uniform in this approximate definition of precession
    # phase
    precession_phase = np.unwrap(np.arctan2(chiA[1], chiA[0]))
    t = binaryBHexp.get_uniform_in_orbits_times(t_orig, precession_phase, 100)
    chiA = [binaryBHexp.spline_interp(t, t_orig, tmp) for tmp in chiA]
    chiB = [binaryBHexp.spline_interp(t, t_orig, tmp) for tmp in chiB]
    chiA = np.array(chiA).T
    chiB = np.array(chiB).T

    # Attaching 3D axis to the figure
    ax = axes3d.Axes3D(fig)

    time_fontsize = 12
    properties_fontsize = 10
    properties_text_yloc = 0.75
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


    arrow_mutation_scale = 20
    marker_alpha = 0.9
    traj_alpha = 0.8
    lines = [\
        # These two are for plotting component BH spins
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['BhA_spin'], \
            zorder=zorder_dict['spin'])), \
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['BhB_spin'], \
            zorder=zorder_dict['spin'])), \
        # These four are for tracing the spin paths
        ax.plot(chiA.T[0,0:1]-1e10, chiA.T[1,0:1], chiA.T[2,0:1], \
            color=colors_dict['BhA_spin'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ax.plot(chiB.T[0,0:1]-1e10, chiB.T[1,0:1], chiB.T[2,0:1], \
            color=colors_dict['BhB_spin'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ]

    # Setting the axes properties

    # This seems to set the actual limits to max_range
    max_range = max(np.max(np.linalg.norm(chiA, axis=1)), \
                np.max(np.linalg.norm(chiB, axis=1)))
    ax.set_xlim3d([-max_range, max_range])
    ax.set_ylim3d([-max_range, max_range])
    ax.set_zlim3d([-max_range, max_range])

    ax.set_xlabel('$\chi_x$', fontsize=label_fontsize)
    ax.set_ylabel('$\chi_y$', fontsize=label_fontsize)
    ax.set_zlabel('$\chi_z$', fontsize=label_fontsize)

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

    frames = range(1, len(t))
    fargs = (lines, t, time_text, chiA, chiB)

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
    parser.add_argument('--RunDir', '-r', type=str, required=True,
        help='Dir containing Horizons.h5. Movie will be saved there')
    parser.add_argument('--fname', '-f', type=str, default="spins",
        help='Name to use for saving movie file without the extension.')
    parser.add_argument('--still_time', default=None, type=float, \
        help='If given, saves a plot of the movie at this time and exits.')

    args = parser.parse_args()
    fig = P.figure(figsize=(5,4))
    line_ani = PrecessionTrack(fig, args.RunDir,
        save_file = args.fname,
        still_time = args.still_time)

    # Might need: conda install -c conda-forge ffmpeg
    Writer = animation.writers['ffmpeg']
    metadata = {
        'artist' : 'Vijay Varma',
        'genre' : 'Physics',
        'subject' : 'Black hole spin evolution for precessing BBH.',
        'copyright' : binaryBHexp.__copyright__,
        }
    writer = Writer(fps=15, metadata=metadata)
    line_ani.save('%s.mp4'%args.fname, writer=writer, dpi=150)
