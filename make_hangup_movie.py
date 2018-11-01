desc = """
Generates a movie demonstrating the orbital hang-up effect.
"""

import numpy as np
import argparse
import os

q = 1.
chi_mag = 0.8

cmdline_format = "./binaryBHexp --q {:.2f} " \
    "--chiA 0 0 {:.2f} --chiB 0 0 {:.2f} --omega_start 1.8e-2 " \
    "--wave_time_series --no_freeze_near_merger --stop_after_ringdown " \
    "--uniform_time_step_size 1 --save_file animations/{}.mp4"

# Generate individual movies
cmdline = cmdline_format.format(q, chi_mag, chi_mag, "hangup_aligned")
os.system(cmdline)

cmdline = cmdline_format.format(q, 0, 0, "hangup_nonspin")
os.system(cmdline)

cmdline = cmdline_format.format(q, -chi_mag, -chi_mag, "hangup_antialigned")
os.system(cmdline)

# Combine into a single mp4
join_cmdline = "ffmpeg -i animations/hangup_aligned.mp4 " \
    + "-i animations/hangup_nonspin.mp4 " \
    + "-i animations/hangup_antialigned.mp4 " \
    + "-filter_complex \"[0:v:0][1:v:0][2:v:0]hstack=inputs=3\" " \
    + "animations/hangup.mp4"
os.system(join_cmdline)

# save stills
still_times = [3, 31.3, 42, 52.25]
for time in still_times:
    os.system("ffmpeg -i animations/hangup.mp4 -ss 00:00:%.3f -vframes 1 animations/hangup_t%d.jpg"%(time, time))
