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
from gwtools import harmonics

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
P.style.use('seaborn')

import black_hole_scattering

Arrow3D = black_hole_scattering.Arrow3D
zorder_dict = black_hole_scattering.zorder_dict
colors_dict = black_hole_scattering.colors_dict

#----------------------------------------------------------------------------
def update_lines(num, lines, t, properties_text, time_text, q, mA, mB, \
        sA, sB, L, J, chiA_nrsur, chiB_nrsur, separation, draw_full_trajectory):
    """ The function that goes into animation
    """
    current_time = t[num]
    time_text.set_text('$t=%.1f\,M$'%current_time)
    properties_text.set_text('$q=%.2f$\n' \
        '$\chi_{A}=[%.2f, %.2f, %.2f]$\n' \
        '$\chi_{B}=[%.2f, %.2f, %.2f]$\n' \
        '$r=%.2f\,M$\n'%(q, \
        black_hole_scattering.make_zero_if_small(chiA_nrsur[num-1][0]), \
        black_hole_scattering.make_zero_if_small(chiA_nrsur[num-1][1]), \
        black_hole_scattering.make_zero_if_small(chiA_nrsur[num-1][2]), \
        black_hole_scattering.make_zero_if_small(chiB_nrsur[num-1][0]), \
        black_hole_scattering.make_zero_if_small(chiB_nrsur[num-1][1]), \
        black_hole_scattering.make_zero_if_small(chiB_nrsur[num-1][2]), \
        separation[num-1]
        ))

    for idx in range(len(lines)):
        line = lines[idx]
        if idx == 0:
            line.set_arrow_at_origin(sA[num-1])
        elif idx == 1:
            line.set_arrow_at_origin(sB[num-1])
        elif idx == 2:
            line.set_arrow_at_origin(L[num-1])
        elif idx == 3:
            line.set_arrow_at_origin(J[num-1])
        else:
            if idx == 4:
                data = sA
            elif idx == 5:
                data = sB
            elif idx == 6:
                data = L
            elif idx == 7:
                data = J

            line.set_data(data.T[0:2, 0:num])
            line.set_3d_properties(data.T[2, 0:num])

    return lines


#----------------------------------------------------------------------------
def PrecessionTrack(fig, q, chiA, chiB, omega_ref=None, \
        draw_full_trajectory=False, save_file=None, still_time=None,  \
        rescale_fig_for_widgets=False):

    mA = q/(1.+q)
    mB = 1./(1.+q)
    chiA = np.array(chiA)
    chiB = np.array(chiB)
    t, chiA_nrsur, chiB_nrsur, L, h_nrsur, BhA_traj, \
        BhB_traj, separation = black_hole_scattering.get_binary_data( \
        q, chiA, chiB, omega_ref)

    sA = chiA_nrsur * mA**2
    sB = chiB_nrsur * mB**2
    J = sA + sB + L

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
        # This is for plotting oribtal angular momentum direction
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['L'], \
            zorder=zorder_dict['L'])), \
        # This is for plotting the total angular momentum direction
        ax.add_artist(Arrow3D(None, mutation_scale=arrow_mutation_scale, lw=3, \
            arrowstyle="-|>", color=colors_dict['J'], \
            zorder=zorder_dict['J'])), \

        # These four are for tracing the spin paths
        ax.plot(sA.T[0,0:1]-1e10, sA.T[1,0:1], sA.T[2,0:1], \
            color=colors_dict['BhA_spin'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ax.plot(sB.T[0,0:1]-1e10, sB.T[1,0:1], sB.T[2,0:1], \
            color=colors_dict['BhB_spin'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ax.plot(L.T[0,0:1]-1e10, L.T[1,0:1], L.T[2,0:1], \
            color=colors_dict['L'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \
        ax.plot(J.T[0,0:1]-1e10, J.T[1,0:1], J.T[2,0:1], \
            color=colors_dict['J'], lw=2, alpha=traj_alpha, \
            zorder=zorder_dict['traj'])[0], \

        ]

    # Setting the axes properties

    # This seems to set the actual limits to max_range
    max_range = 0.6
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

    ax.set_title('NRSur7dq2', fontsize=time_fontsize, x=0.74, y=0.99)

    #NOTE: There is a glitch if I don't skip the first index
    frames = range(1, len(t[t<=0]))

    fargs = (lines, t, properties_text, time_text, q, mA, mB, \
        sA, sB, L, J, chiA_nrsur, chiB_nrsur, separation, draw_full_trajectory)

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
    parser.add_argument('--draw_full_trajectory', default=False, \
        action='store_true', \
        help='If given, draws the entire trajectories of the components. ' \
        'Else only retains the last 3/4th of an orbit.')
    parser.add_argument('--still_time', default=None, type=float, \
        help='If given, saves a plot of the movie at this time and exits.')

    args = parser.parse_args()

    fig = P.figure(figsize=(5,4))
    line_ani = PrecessionTrack(fig, args.q, args.chiA, args.chiB,
        omega_ref = args.omega_ref,
        draw_full_trajectory = args.draw_full_trajectory,
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
        if extension == 'gif':
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
        fig.canvas.mpl_connect('button_press_event', onClick)
        P.show()
