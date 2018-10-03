#!/usr/bin/env python3
"""
Writes out the commands to execute in order to make a movie of the
sinusoidal dependence of kick velocity on in-plane components of spin
in the left-right configuration.

Example usage:
  ./make_sine_movies_script.py > make_sine_movies.sh
  bash make_sine_movies.sh

TODO add command line arguments to control number of movies etc.
"""

import numpy as np

q = 1.
chi_mag = 0.8

n_mov = 9
delta_alpha = 9.*np.pi/8.

alphas = np.linspace(0., np.pi, n_mov) + delta_alpha

base_filename = "animations/sinus"

cmdline_format = "python black_hole_scattering.py --q {:.2f} --chiA {:.2f} {:.2f} {:.2f} --chiB {:.2f} {:.2f} {:.2f} --save_file {}{}.mp4"

for i, alpha in enumerate(alphas):
    chiA = chi_mag * np.array([np.cos(alpha), np.sin(alpha), 0])
    chiB = -chiA
    cmdline = cmdline_format.format(q, chiA[0], chiA[1], chiA[2], chiB[0], chiB[1], chiB[2], \
                                    base_filename, i)
    print(cmdline)

join_cmdline = "ffmpeg "+ "".join(["-i {}{}.mp4 ".format(base_filename,i) for i in range(n_mov) ]) + \
    " -filter_complex \"" + "".join(["[{}:v:0]".format(i) for i in range(n_mov)]) + \
    "hstack=inputs={}\" animations/sinus_all.mp4".format(n_mov)
print(join_cmdline)
