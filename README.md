# BELLHOP
A mirror of the original Fortran BELLHOP underwater acoustics simulator, with
numerical properties and robustness improved and bugs fixed.

### Impressum

Copyright (C) 2021-2022 The Regents of the University of California \
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu \
Copyright (C) 1983-2020 Michael B. Porter

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Summary of changes

### iSegz, iSegr initialization

In `BELLHOP`, `iSegz` and `iSegr` are initialized (to 1) only once globally. In
some cases (a source on a boundary), in conjunction with other subtle issues
(see below), the missing initialization caused the initial state of one ray to
be dependent on the final state of the previous ray, which is obviously
non-physical and makes reproducing certain results impossible. This has been
fixed by reinitializing them before every ray.

### Boundary stepping changes

How edge cases are handled seems like it should only affect a small fraction of
results in any program. However, in `BELLHOP`, the step size for the numerical
integration is typically set to a large value, and is then reduced so that the
step lands on the nearest "boundary". (Here, this refers not only to the ocean
surface and floor, but to any "plane" within the ocean where the SSP changes,
where the bottom has a corner/vertex in its definition, and so on.) Thus, for
typical runs, every step of every ray in `BELLHOP` lands on an edge case.
Depending on the exact set of floating-point operations done and the specific
inputs, these steps "randomly" land just before or just after the boundary. This
means the step "randomly" ends up in different segments, which leads to further
differences down the line. While these differences typically have a very small
effect on the overall ray trajectory, they make comparing the results to [the
C++/CUDA version](https://github.com/A-New-BellHope/bellhopcuda) very difficult,
and may in some cases make reproducing results between runs of the Fortran
program difficult. This has been fixed by overhauling the system which moves the
numerical integration to the boundary to give predictable behavior. This also
involved changes to reflections (landing on the top or bottom boundary); now
consistently only two identical points are generated at reflections, one with
the incoming tangent and one with it outgoing.

### Polarity flipping

In `InfluenceCervenyRayCen`, `rnV` was negated for image 2 and again for image 3.
If `Nimage` was set to 2, this means `rnV` would change sign every step. mbp
had implemented this differently (without the bug) in `InfluenceCervenyCart`,
so that implementation was brought to this function.

### SGB eigenrays and arrivals

`InfluenceSGB` had code to support eigenrays, but the condition to check the
distance from the receiver was commented out, meaning that every ray would be
an eigenray. Furthermore, it assumed that any run other than eigenrays was
coherent TL, ignoring incoherent, semi-coherent, and arrivals. This function
was reworked to use `ApplyContribution` to support all of these run types.

### Shallow angle range

In `InfluenceGeoGaussianCart`, a ray is considered to be at a shallow angle if
its angle is less than or equal to 60 degrees, and a steep ray otherwise. The
handling of these two cases is completely different, i.e. there is a
discontinuity as the ray angle passes 60 degrees. However, some environment
files trace rays at round number angles, including exactly 60 degrees.
Differences in the set of floating-point operations used can cause these rays
to be considered on one side or the other of this boundary, leading to different
results. This cannot be fully fixed without changing the physics of the
simulation (i.e. remove the discontinuity), so instead the discontinuity has
been moved slightly away from the round number to increase the chances of
consistent results.

### Missing initialization on some beam parameters

The beam parameters `Nimage`, `iBeamWindow`, and `Component` were not
initialized in code, and not all `BELLHOP`-provided test `.env` files initialize
`Component`. The value of uninitialized variables in Fortran is not defined.
Reasonable initializations have been added to these parameters.

### Other

Minor fixes which should not affect results unless you're doing something
strange. See commit logs.

# Detailed information on some changes

### iSegz, iSegr initialization

In `BELLHOP`, `iSegz` and `iSegr` are initialized (to 1) only once globally.
This means their initial state for each ray is their final state from the
previous ray. *Normally, this shouldn't matter, because the segment indices are
recomputed if they are wrong; but...*

SSP checks whether the current depth and range are within the current segment,
using an interval closed on both sides (including both endpoints). If the
segment is incorrect, the correct segment is searched for, but this search uses
a half-open interval (including the lower endpoint but excluding the upper one).
It often happens that the depth or range value is exactly on an endpoint of a
segment, because usually both source positions and depth segment boundaries are
artificially set to nice round numbers in the environment files. In these cases,
this algorithm will use the lower segment this endpoint belongs to if the
current segment index (`iSegz` or `iSegr`) is the lower segment, and it will use
the upper segment if the segment index was anything else and it had to search
for the correct segment. So, in these cases, the segment the ray starts in
depends on the previous state of `iSegz` or `iSegr`, i.e. the segment the
previous ray finished in. *Normally, this shouldn't matter, because the position
is exactly on the boundary of these two segments, so it should be equally valid
in either of them; but...*

The numerical integration step algorithm takes a trial step forward of a default
step size. If this step puts the position into a new segment, it reduces the
step size so that it falls on the segment boundary instead. In the cases above,
if the position is on the lower boundary of the upper segment, this step will
remain in the same segment, but if it's on the upper boundary of the lower
segment, this step will cross into the upper segment and therefore be reduced
to land on the boundary again, for a step size of zero. There is a failsafe
which ensures a minimum step size, placing the ray slightly inside the upper
segment. This creates an extra, very small (~1e-5 m), step, or not, depending on
which segment the previous ray happened to end in. *Normally, this shouldn't
matter, because one extra infinitesimally small step should not substantially
change the field results; but...*

The SGB (Bucker's Simple Gaussian Beams) transmission loss influence algorithm
has a bug, acknowledged by mbp in code comments, where every step is assumed to
be of the full, maximum step size, which is often 500 m or larger. Thus, one
extra infinitesimally small step is actually one extra full-size step, which
substantially impacts the field results from this beam. Again, this change
happens "randomly" based on the final state of the previous ray/beam traced.

`iSegz` and `iSegr` have been re-initialized to 1 before every ray/beam is
traced, removing this "random" inconsistent behavior. The SGB bug still stands,
as fixing it would mean changing the physics, but at least it produces
reproducible results now.

### Boundary stepping changes

The numerical integration in `BELLHOP` operates by taking a trial step along
the ray with a large step size, and then reducing the step size to land on the
first "interface or boundary crossing". Then it evaluates the field parameters
at this trial location, and adjusts the angles for the trial step and tries
again. Finally, it interpolates between the two trials based on the field
parameters and commits a step in that resulting direction. Steps are reduced
to land on whichever of the following comes first (along the ray):

- Plane where SSP changes along depth
- Top crossing
- Bottom crossing
- Top segment crossing (i.e. where the top definition changes, either due to a
vertex in its height, or changed parameters; this rarely occurs for the top,
but often for the bottom)
- Bottom segment crossing
- Plane where SSP changes along range

Each of these boundaries produces a step size, and the smallest one is used.
Then this step size is used to update the position (and other parameters). Due
to floating-point numbers having limited precision, in practice it is "random"
whether the resulting position is slightly before or slightly after the actual
boundary. For example, if there is a depth SSP change at 1000 m, the position
may end up being 999.999999994 m or 1000.0000000002 m. This often occurs even
though, in this example (as is often the case for values provided in a config
file), the "correct" value (1000.0) is precisely expressible as a floating-
point number. This is not actually non-deterministic behavior--`BELLHOP` will
produce binary identical outputs for identical input files. However, for some
set of rays shot toward the same boundary, some set of them will land on each
side of the boundary, without any discernible pattern to this behavior. And a
version of `BELLHOP` built with a different compiler or on a different machine,
or of course [our C++/CUDA version](https://github.com/A-New-BellHope/bellhopcuda),
will produce a different pattern.

Despite this behavior only causing an infinitesimal difference in the resulting
ray position, this difference gets amplified by `BELLHOP`, and can lead to a
number of inconsistent results later. First, the ray which landed on the near
side of the boundary will take its next step, hit the same boundary again, but
cross it due to the enforced minimum step size. In addition to this ray now
having an extra step, the difference in its position (compared to the ray which
happened to land on the other side of the boundary) has now been amplified from
the scale of about 1e-16 from the floating-point error, to the minimum step size
of about 1e-5. While this is still not a significant difference in the ray
trajectory, some of the issues this can cause are:
- In ray mode, the number of steps in the ray is different; there is a different
number of superfluous near-duplicated points.
- In TL, eigenrays, or arrivals, a ray being considered close enough to a
receiver to influence it is a hard edge for certain influence types. For runs
with tens of thousands each of rays and receivers, it is not rare that a ray
will be close enough to a receiver to "count" in one case, and too far in the
other case. If this ray happens to be the only ray to influence this receiver,
the relative difference between the finite field and zero is an *infinite* error!
- If the ray approaches a corner between two different boundaries (e.g. SSP
depth and range)--as is not rare in the environment files where everything is
round numbers!--the ray may hit one boundary or the other first depending on
the slight change to its trajectory. And depending on the boundaries, this can
lead to much larger changes to the ray trajectory. For example, suppose a SSP
range interface occurs at 10 km, and a vertex of the bottom is also at 10 km,
with a seamount beginning to rise. If the ray hits the bottom first, it will
reflect and miss the seamount, whereas if it hits the range interface first, it
will then hit the seamount and reflect in a completely different direction.
- There are several stopping conditions for rays, for example the level becoming
too low, but crucially the last step which violated the stopping condition is
still written. If a ray stops after a step to a certain boundary, whether
receivers on that boundary are considered for TL / eigenrays / arrivals is
determined by whether the ray is behind or ahead of that boundary. Of course,
since environment files use round numbers, receivers being exactly on boundaries
is very common! This can substantially change their outputs.

These issues are fixed with two main changes. First, a function `StepToBdry2D`
is added which corresponds to `ReduceStep2D`. This function ensures that after
the step, the resulting position is *exactly* on the boundary (and if the value
is not precisely representable in floating-point, that it is the same value as
already stored in memory for the boundary). Second, the functions which update
which segment a position is in have been changed to use the ray tangent to
resolve edge cases, instead of choosing between less-than-or-equal-to versus
less-than and such. The ray is always placed into the later segment: the
segment such that it can move forward a nontrival distance and remain in the
same segment. The idea is, since every step must be an edge case, they are
*exactly* on the edge, and then that exactly-on-the-edge position are handled in
a consistent manner.

# Other information

See index.htm for information from the original repo.

Code retrieved 12/17/21 from http://oalib.hlsresearch.com/AcousticsToolbox/ ,
and was claimed to have been last updated 11/4/20.

Files pertaining to the other simulators (Krakel, Kraken, KrakenField, Scooter)
have been removed.

The Makefile has been set up to build with gfortran by default (Linux or
Windows with mingw etc.)
