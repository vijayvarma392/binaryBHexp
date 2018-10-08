<script src="http://api.html5media.info/1.1.8/html5media.min.js"></script>

# Black hole scattering
Visualization of the binary black hole scattering process.

## Introduction

When two black holes (BHs) interact, they spiral towards each other and
eventually merge, radiating gravitational waves in the process. The end state
is a single black hole that is only characterized by its mass, angular momentum
and linear momentum. Therefore, this process can be thought of as a scattering
process.

Binary black hole scattering is a very complicated process that requires
numerical simulations to accurately predict the outcome. These simulations are
very expensive, taking about a month on a super computer to run. In recent
years, several "surrogate" models for these simulations have been developed.
These models are trained against a large number of numerical simulations, using
some fancy interpolation methods. In the end, the surrogate can reproduce the
result of the simulation as accurately as the simulation itself, but in a
fraction of a second, and on your laptop, rather than a month on a super
computer!

In this package, we use these surrogate models to visualize the binary black
hole scattering process. This has been done before using numerical simulations,
but there are only a finite number of these simulations, at discrete points in
the parameter space. With this package, you can generate visualizations at any
point in the parameter space, in a few seconds, without needing to do a full
numerical simulation!

## Visualizations

Enough background, let's see some videos! 

### Precessing binary black hole

<video src="animations/precessing.mp4" width="500" controls preload loop></video>


Here we see a precessing binary black hole merger. The black holes are shown as
circular markers, with arrows indicating their spins. The orbital angular
momentum direction is shown by the pink arrow at the origin. The colors in the
planes indicate the value of the plus polarization of the GW as seen by an
observer at that location; red means positive and blue means negative, notice
the quadrupolar pattern of the radiation. In the plot at the bottom, we show
the plus and cross polarizations as seen from the camera viewing angle.  The BH
spins are misaligned w.r.t the angular momentum direction. Therefore, the
spins, the angular momentum, and the orbital plane all precess during the
inspiral. It is interesting to see how the remnant spin direction is very close
to the direction of the orbital angular momentum near merger.

<video src="animations/precessing_auto_rotate.mp4" width="500" controls preload loop></video>

Here we show the same animation, but with varying camera angle. We project the
wavefrom only on the bottom plane. In the plot at the bottom we see how the
waveform changes based on the viewing angle; this is because the different
modes (spin-weighted spherical harmonic modes) combine with different weights
based on the observer viewing angle. Notice how the GW signal is strongest
along the direction of orbital angular momentum and weakest in perpendicular
directions.

See
[PhysRevD.49.6274](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.49.6274)
to learn more about precessing binary black holes.

### Aligned-spin binary black hole

<video src="animations/aligned.mp4" width="500" controls preload loop></video>

Here the spins are aligned w.r.t the angular momentum direction, and there
is no precession. Due to symmetries, the remnant spin is also aligned,
while the remnant velocity is restricted to the orbital plane.

### Super-kick

<video src="animations/super_kick.mp4" width="500" controls preload loop></video>

In this case, notice how the spins line up in a single plane close to the
merger. This is believed to be necessary to cause super-kicks, or very high
remnant velocities. Notice how fast the final black hole flies away compared
to the above cases.

This phenomenon was discovered through numerical simulations. See
[gr-qc/0702133](https://arxiv.org/abs/gr-qc/0702133),
[gr-qc/0702052](https://arxiv.org/abs/gr-qc/0702052), and
[arxiv:1802.04276](https://arxiv.org/abs/1802.04276) for more details.

### Details about these videos

A black hole is indicated by a circular marker with an arrow attached. The
radius of the marker is proportional to the Kerr horizon radius of the BH. The
arrow direction is along the spin direction of the BH, while the arrow
magnitude is proportional to its Kerr spin parameter, a. The pink arrow at the
origin indicates the angular momentum direction (not magnitude, magnitude of
the arrow is fixed here). The time is shown at the bottom left, t=0 corresponds
to the peak of the total waveform amplitude. All quantities are in units of
total mass M.

Before the merger, the trajectories and spins are generated using the
[NRSur7dq2](https://pypi.org/project/NRSur7dq2/) surrogate waveform model.  The
mass ratio and instantaneous spins are shown at the top left. NRSur7dq2
predicts the coprecessing frame quaternions, orbital phase, and component spins
during the inspiral, along with the waveform. We estimate the separation using
the 3.5 PN accurate relation for separation as a function of orbital frequency.
The movie time-steps during this stage are non-linearly related to the
simulation time: there are 30 frames per each orbit.

After the merger, the remnant properties are generated using the
[surfinBH](https://pypi.org/project/surfinBH/) package. The final BH's mass,
spin and recoil kick velocity is shown at the top left. Note that in an actual
simulation the kick is accumulated during the late inspiral and merger, causing
the center of mass to drift. Here we ignore this and generate the remnant with
the final kick velocity at the origin. After the ringdown is complete, we
switch our time steps to demonstrate the recoil kicks: Each frame of the video
during this stage corresponds to a time step of 100M in simulation time.

## Generating new visualizations

Since these videos don't require numerical simulations, you can generate as
many as you like, and at any point in the parameter space.  You can do so using
the script ```black_hole_scattering.py``` in this repository.

### Usage
```python
python black_hole_scattering.py --q 2 --chiA 0.2 0.7 -0.1 --chiB 0.2 0.6 0.1
```

Do ```python black_hole_scattering.py -h``` for more options.

### Dependencies

All of these can be installed through pip or conda.
* [NRSur7dq2](https://pypi.org/project/NRSur7dq2) (at least 1.0.5)  
* [surfinBH](https://pypi.org/project/surfinBH/)
* [gwtools](https://pypi.org/project/gwtools/)

## Credits
The code is developed and maintained by [Vijay
Varma](http://www.tapir.caltech.edu/~vvarma/).  Please credit me if you use
these visualizations in your work, presentations or outreach.  Please, report
bugs to
[&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;](mailto:&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;).  
