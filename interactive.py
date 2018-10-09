""" Interactive widgets for running black_hole_scattering.py
"""
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.widgets import Slider, Button
import matplotlib.pyplot as P
import numpy as np

import black_hole_scattering


def run_animation():

    q = sl_q.val
    chiAmag = sl_chiAmag.val
    chiAth = sl_chiAth.val
    chiAph = sl_chiAph.val
    chiBmag = sl_chiBmag.val
    chiBth = sl_chiBth.val
    chiBph = sl_chiBph.val

    chiA = [chiAmag*np.sin(chiAth)*np.cos(chiAph),
            chiAmag*np.sin(chiAth)*np.sin(chiAph),
            chiAmag*np.cos(chiAth)]
    chiB = [chiBmag*np.sin(chiBth)*np.cos(chiBph),
            chiBmag*np.sin(chiBth)*np.sin(chiBph),
            chiBmag*np.cos(chiBth)]
    ani = black_hole_scattering.BBH_scattering(fig, q, chiA, chiB, \
        wave_time_series=True, auto_rotate_camera=False, \
        rescale_fig_for_widgets=True)
    return ani


fig = P.figure(figsize=(7,5.5))

axcolor = 'lightgoldenrodyellow'
# axes for sliders
ax_q = P.axes([0.75, 0.8, 0.15, 0.03], facecolor=axcolor)
ax_chiAmag = P.axes([0.75, 0.75, 0.15, 0.03], facecolor=axcolor)
ax_chiAth = P.axes([0.75, 0.70, 0.15, 0.03], facecolor=axcolor)
ax_chiAph = P.axes([0.75, 0.65, 0.15, 0.03], facecolor=axcolor)
ax_chiBmag = P.axes([0.75, 0.60, 0.15, 0.03], facecolor=axcolor)
ax_chiBth = P.axes([0.75, 0.55, 0.15, 0.03], facecolor=axcolor)
ax_chiBph = P.axes([0.75, 0.50, 0.15, 0.03], facecolor=axcolor)

# run button
ax_run = P.axes([0.75, 0.4, 0.1, 0.03])
button = Button(ax_run, 'Run', color=axcolor, hovercolor='0.975')

# define sliders
sl_q = Slider(ax_q, '$q$', 1, 2, valinit=1.34, valstep=0.01)
sl_chiAmag = Slider(ax_chiAmag, '$\chi_1$', 0, 0.8, valinit=0.8, valstep=0.01)
sl_chiAth = Slider(ax_chiAth, '$\\theta_1$', 0, 2*np.pi, valinit=1.1, \
    valstep=0.01)
sl_chiAph = Slider(ax_chiAph, '$\phi_1$', -np.pi, np.pi, valinit=-0.4, \
    valstep=0.01)
sl_chiBmag = Slider(ax_chiBmag, '$\chi_2$', 0, 0.8, valinit=0.8, valstep=0.01)
sl_chiBth = Slider(ax_chiBth, '$\\theta_2$', 0, 2*np.pi, valinit=1.1, \
    valstep=0.01)
sl_chiBph = Slider(ax_chiBph, '$\phi_2$', -np.pi, np.pi, valinit=-0.4, \
    valstep=0.01)

def run(event):
    global ani
    ani = run_animation()

button.on_clicked(run)
P.show()
