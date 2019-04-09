import os
import numpy as np

prec_cmd = "./binaryBHexp --q 2 --chiA 0.3 0.7 -0.1 --chiB 0.3 0.6 0.1 --omega_ref 1.8e-2 --no_surrogate_label --no_freeze_near_merger --save_file animations/precessing_dumbed_down.mp4"

prec_rot_cmd = "./binaryBHexp --q 2 --chiA 0.3 0.7 -0.1 --chiB 0.3 0.6 0.1 --omega_ref 1.8e-2 --auto_rotate_camera --no_surrogate_label --no_freeze_near_merger --save_file animations/precessing_auto_rotate_dumbed_down.mp4"

super_kick_cmd = "./binaryBHexp --q 1.34 --chiA 0.62 -0.27 0.34 --chiB -0.62 0.27 0.34 --no_surrogate_label --no_freeze_near_merger --save_file animations/super_kick_dumbed_down.mp4"

aligned_cmd = "./binaryBHexp --q 2 --chiA 0 0 -0.6 --chiB 0 0 0.8 --omega_ref 1.8e-2 --no_surrogate_label --no_freeze_near_merger --save_file animations/aligned_dumbed_down.mp4"

os.system(aligned_cmd)
os.system(prec_cmd)
os.system(prec_rot_cmd)
os.system(super_kick_cmd)
