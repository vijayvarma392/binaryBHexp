### waveform projection stills
wave_proj_cmd() { eval "python black_hole_scattering.py --q 2 --chiA 0.3 0.7 -0.1 --chiB 0.3 0.6 0.1 --omega_ref 1.8e-2  --wave_time_series --auto_rotate_camera --save_file stills/waveform_projection.mp4 --still_time $1";}
wave_proj_cmd -3100
wave_proj_cmd -3450
wave_proj_cmd -3800

#### Super-kick stills
super_kick_cmd() { eval "python black_hole_scattering.py --q 1.34 --chiA 0.62 -0.27 0.34 --chiB -0.62 0.27 0.34 --project_on_all_planes --save_file stills/super_kick.gif --still_time $1";}
super_kick_cmd -2600
super_kick_cmd -100
super_kick_cmd 15
super_kick_cmd 2280
