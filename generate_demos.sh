python black_hole_scattering.py --q 2 --chiA 0.3 0.7 -0.1 --chiB 0.3 0.6 0.1 --omega_ref 1.8e-2 --wave_time_series --project_on_all_planes --save_file animations/precessing.gif

python black_hole_scattering.py --q 2 --chiA 0.3 0.7 -0.1 --chiB 0.3 0.6 0.1 --omega_ref 1.8e-2 --wave_time_series --auto_rotate_camera --save_file animations/precessing_auto_rotate.gif 

python black_hole_scattering.py --q 1.34 --chiA 0.62 -0.27 0.34 --chiB -0.62 0.27 0.34 --wave_time_series --project_on_all_planes --save_file animations/super_kick.gif

python black_hole_scattering.py --q 2 --chiA 0 0 -0.6 --chiB 0 0 0.8 --omega_ref 1.8e-2 --wave_time_series --project_on_all_planes --save_file animations/aligned.gif
