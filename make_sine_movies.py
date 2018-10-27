desc = """
Generates a movie of the sinusoidal dependence of kick velocity on in-plane
components of spin in the left-right configuration.
"""

import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description=desc,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--num_mov', type=int, default=5,
    help='Number of movies. Should be odd.')

args = parser.parse_args()
if args.num_mov%2 != 1:
    raise Exception('num_mov should be odd.')


q = 1.
chi_mag = 0.8

delta_alpha = 9.*np.pi/8.
alphas = np.linspace(0., np.pi, args.num_mov) + delta_alpha

temp_dir = 'tmp_sine_kicks'
os.system('mkdir -p %s'%temp_dir)

base_filename = "%s/sinus"%temp_dir
cmdline_format = "python black_hole_scattering.py --q {:.2f} " \
    "--chiA {:.2f} {:.2f} {:.2f} --chiB {:.2f} {:.2f} {:.2f} " \
    "--save_file {}{}.mp4"

# Generate individual movies
for i, alpha in enumerate(alphas):
    chiA = chi_mag * np.array([np.cos(alpha), np.sin(alpha), 0])
    chiB = -chiA
    cmdline = cmdline_format.format(q, chiA[0], chiA[1], chiA[2], \
        chiB[0], chiB[1], chiB[2], base_filename, i)
    os.system(cmdline)

# Combine all movies into a single mp4
join_cmdline = "ffmpeg " \
    + "".join(["-i {}{}.mp4 ".format(base_filename,i) \
        for i in range(args.num_mov) ]) \
    + " -filter_complex \"" + "".join(["[{}:v:0]".format(i) \
        for i in range(args.num_mov)]) \
    + "hstack=inputs={}\" animations/sine_kicks.mp4".format(args.num_mov)
os.system(join_cmdline)

### Make shorter version of video for gif
##os.system("ffmpeg -ss 26 -i animations/sine_kicks.mp4 -c copy " \
##    "%s/short_sine_kicks.mp4"%temp_dir)
##
### Convert to gif
### https://askubuntu.com/questions/648603/how-to-create-an-animated-gif-from-mp4-video-via-command-line
##os.system("ffmpeg -i %s/short_sine_kicks.mp4  -r 5 '%s"%(temp_dir, temp_dir) \
##        + "/frame-%04d.png'")
##os.system("convert -delay 20 -loop 0 %s/*.png animations/sine_kicks.gif"\
##    %temp_dir)
