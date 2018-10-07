desc = """
Generates a movie of the sinusoidal dependence of kick velocity on in-plane
components of spin in the left-right configuration.
"""

import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description=desc,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--num_mov', type=int, default=5,
    help='Number of movies. Should be odd.')

args = parser.parse_args()
if args.num_mov%2 != 1:
    raise Exception('num_mov should be odd.')


q = 1.
chi_mag = 0.8

delta_alpha = 9.*np.pi/8.
alphas = np.linspace(0., np.pi, args.num_mov) + delta_alpha

out_dir = 'animations/sine_kicks'
os.system('mkdir -p %s'%out_dir)

base_filename = "%s/sinus"%out_dir
cmdline_format = "python black_hole_scattering.py --q {:.2f} " \
    "--chiA {:.2f} {:.2f} {:.2f} --chiB {:.2f} {:.2f} {:.2f} " \
    "--save_file {}{}.mp4"

for i, alpha in enumerate(alphas):
    chiA = chi_mag * np.array([np.cos(alpha), np.sin(alpha), 0])
    chiB = -chiA
    cmdline = cmdline_format.format(q, chiA[0], chiA[1], chiA[2], \
        chiB[0], chiB[1], chiB[2], base_filename, i)
    os.system(cmdline)

join_cmdline = "ffmpeg " \
    + "".join(["-i {}{}.mp4 ".format(base_filename,i) \
        for i in range(args.num_mov) ]) \
    + " -filter_complex \"" + "".join(["[{}:v:0]".format(i) \
        for i in range(args.num_mov)]) \
    + "hstack=inputs={}\" animations/sinus_all.mp4".format(args.num_mov)
os.system(join_cmdline)
