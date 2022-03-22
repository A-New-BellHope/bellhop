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
some cases (a source on a boundary), in conjunction with other subtle issues,
the missing initialization caused the initial state of one ray to be dependent
on the final state of the previous ray, which is obviously non-physical and
makes reproducing certain results impossible. This has been fixed by
reinitializing them before every ray.

### Boundary stepping changes

How edge cases are handled seems like it should only affect a small fraction of
results in any program. However, in `BELLHOP`, the step size for the numerical
integration is typically set to a large value, and is then reduced so that the
step lands on the nearest "boundary". (Here, this refers not only to the ocean
surface and floor, but to any "plane" within the ocean where the SSP changes,
where the bottom has a corner/vertex in its definition, and so on.) Thus,
for typical runs, every step of every ray in `BELLHOP` lands on an edge case.
Depending on the exact set of floating-point operations done and the specific
inputs, these steps "randomly" land just before or just after the boundary.
This means the step ends up in different regions, which leads to further
differences down the line. While these differences typically have a very small
effect on the overall ray trajectory, they make comparing the results to [the
C++/CUDA version](https://github.com/A-New-BellHope/bellhopcuda) very difficult,
and may in some cases make reproducing results between runs of the Fortran
program difficult. This has been fixed by overhauling the system which moves
the numerical integration to the boundary to give predictable behavior. This
also involved changes to reflections (landing on the top or bottom boundary);
now a consistent number of duplicate ray points are produced at the reflection
points.

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

### Polarity flipping

In `InfluenceCervenyRayCen`, `rnV` was negated for image 2 and again for image 3.
If `Nimage` was set to 2, this means `rnV` would change sign every step. mbp
had implemented this differently (without the bug) in `InfluenceCervenyCart`,
so that implementation was brought to this function.

### Missing initialization on some beam parameters

The beam parameters `Nimage`, `iBeamWindow`, and `Component` were not
initialized in code, and not all `BELLHOP`-provided test `.env` files initialize
`Component`. The value of uninitialized variables in Fortran is not defined.
Reasonable initializations have been added to these parameters.

### SGB eigenrays and arrivals

`InfluenceSGB` had code to support eigenrays, but the condition to check the
distance from the receiver was commented out, meaning that every ray would be
an eigenray. Furthermore, it assumed that any run other than eigenrays was
coherent TL, ignoring incoherent, semi-coherent, and arrivals. This function
was reworked to use `ApplyContribution` to support all of these run types.

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
but at least it produces reproducible results now.


# Other information

See index.htm for information from the original repo.

Code retrieved 12/17/21 from http://oalib.hlsresearch.com/AcousticsToolbox/
claimed to have been updated 11/4/20.

Files pertaining to the other simulators (Krakel, Kraken, KrakenField, Scooter)
have been removed.

The Makefile has been set up to build with gfortran by default (Linux or
Windows with mingw etc.)
