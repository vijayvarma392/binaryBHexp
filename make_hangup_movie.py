desc = """
Generates a movie demonstrating the orbital hang-up effect.
"""

import numpy as np
import argparse
import os

q = 1.
chi_mag = 0.8

cmdline_format = "./binaryBHexp --q {:.2f} " \
    "--chiA 0 0 {:.2f} --chiB 0 0 {:.2f} --wave_time_series " \
    "--save_file animations/{}.mp4"

# Generate individual movies
cmdline = cmdline_format.format(q, chi_mag, chi_mag, "hangup_aligned")
os.system(cmdline)

cmdline = cmdline_format.format(q, -chi_mag, -chi_mag, "hangup_antialigned")
os.system(cmdline)

# Combine into a single mp4
join_cmdline = "ffmpeg -i animations/hangup_aligned.mp4 " \
    + "-i animations/hangup_antialigned.mp4 " \
    + "-filter_complex \"[0:v:0][1:v:0]hstack=inputs=2\" " \
    + "animations/hangup.mp4"
os.system(join_cmdline)
