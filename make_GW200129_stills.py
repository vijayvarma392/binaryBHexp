#############################################################################
##
##      Filename: make_GW200129_stills.py
##
##      Author: Vijay Varma
##
##      Created: 20-04-2022
##
##      Description:
##
#############################################################################

import numpy as np
import os

base_cmd = "./binaryBHexp --q 2.147 --chiA -0.469 0.806 0.218 --chiB 0.013 -0.413 -0.550 --no_surrogate_label --no_freeze_near_merger"

times = [-3000, -2000, -1000, -500, -100, 0, 500, 1000, 2000]
for idx, t in enumerate(times):
    os.system(base_cmd + f" --save_file GW200129_{idx}.png --still_time {t}")
