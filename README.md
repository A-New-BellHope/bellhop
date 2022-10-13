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

## Edge case issues

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

### Strange behavior when escaping the altimetry/bathymetry beam box in 3D

There are mulitple strange behaviors in the code for when rays escape the region
where altimetry and bathymetry is defined in 3D (`GetTopSeg3D` and
`GetBotSeg3D`).
- If the ray escapes the box, a warning is printed, but segment information
  (e.g. `IsegTopx`, `xTopSeg`, `Top_tri_n`) is not updated. But, these values
  are not fully initialized between rays, meaning that it's possible for the
  first step of one ray to use values from the last step of the last ray.
- If the ray escaped the box to the negative side in either dimension, or
  the results are NaN (see below), the ray is then terminated in `TraceRay`
  without the current step. However, if the ray escaped to the positive side,
  the ray is terminated as part of the normal stopping conditions, including
  the current step.
- There is a check for the depth being NaN, which sets the segment to 0
  (invalid). But it is not clear what code catches that condition, or if the
  errors arising from the invalid segment are at all better than the NaN (which
  could only have arisen from an error in the first place).
- The original logic allowed for steps to the edge of the same segment they are
  in, which has been generally changed as discussed above to step into the next
  segment according to the tangent.
- In 2D, the altimetry / bathymetry are extended to very large values so these
  conditions never get hit, whereas in Nx2D and 3D they are not so almost all
  rays escape the box.

These functions have been rewritten as follows:
- There is an override so that being at the outer edge of a segment is OK for
  the last segment, just for one step.
- If the position is actually outside the segment, print a warning and put it
  in the nearest valid segment.
- Changed the stopping condition to stop if on the boundary and pointing
  outwards.
- Values all initialized between rays.

### `RToIR` step-to-edge-case issue

The code `ir = MAX( MIN( INT( ( r - Pos%Rr( r ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )`
or similar is used in various Influence functions. The `INT( )` operation here
is another edge case which steps often step exactly to. For example, suppose `r`
is 3000 after a step to a particular segment boundary, `Pos%Rr( r ) == 0`, and
`Pos%Delta_r == 500`. The argument to `INT( )` may end up being 5.99999999999998
or 6.000000000000002 depending on the exact value of `r`, but this would result
in different values of `ir` and therefore substantially different results. This
code pattern has been replaced with the function `RToIR`, which will "snap to"
each integer value of the argument, to ensure these cases are handled
consistently.


## Uninitialized values

### iSegz, iSegr initialization

In `BELLHOP`, `iSegz` and `iSegr` are initialized (to 1) only once globally. In
some cases (a source on a boundary), in conjunction with other subtle issues
(see below), the missing initialization caused the initial state of one ray to
be dependent on the final state of the previous ray, which is obviously
non-physical and makes reproducing certain results impossible. This has been
fixed by reinitializing them before every ray.

### Missing initialization on some beam parameters

The beam parameters `Nimage`, `iBeamWindow`, and `Component` were not
initialized in code, and not all `BELLHOP`-provided test `.env` files initialize
`Component`. The value of uninitialized variables in Fortran is not defined.
Reasonable initializations have been added to these parameters.

### Uninitialized read in 3D Step phase 2

`Step3D` phase 2 had `RayNormal_unit` using `ray2%phi` which is uninitialized
memory (or left over from the previous ray). This was changed to `ray1%phi`.
Also, it used `RayNormal_unit` despite this being for cases where the tangent
was already normalized, but computed the normalized tangent immediately before.
So this was changed to `RayNormal`.

### Uninitialized variables in `PickEpsilon`

Multiple codepaths do not set output values in `PickEpsilon`, leading to
undefined behavior. For 2D, `HalfWidth` and `epsilonOpt` are unset for:
- "Double" Cerveny style beams (`BeamType( 2 : 2 ) == 'C'`, in addition to
  `BeamType( 1 : 1 )` being `'C'` or `'R'` for Cerveny beams)
- Geometric hat beams specified as `'^'` instead of `'G'`

For 3D / Nx2D, `epsilonOpt` is unset for:
- double Cerveny beams
- WKB beams

For both 2D and 3D, beams not following any of the beam types or beam width
types are also unset, though this leads to other errors as well.

These issues went undetected because the results of `PickEpsilon` are actually
only used in conjunction with Cerveny beams, and the Cerveny 3D influence
function `Influence3D` was removed from `BELLHOP3D` before the earliest version
we have.


## Missing checks

### SGB eigenrays and arrivals

`InfluenceSGB` had code to support eigenrays, but the condition to check the
distance from the receiver was commented out, meaning that every ray would be
an eigenray. Furthermore, it assumed that any run other than eigenrays was
coherent TL, ignoring incoherent, semi-coherent, and arrivals. This function
was reworked to use `ApplyContribution` to support all of these run types.

### Missing divide by zero check in Nx2D reflect

In `Reflect2D` (Nx2D only version), for `Beam%Type( 4 : 4 ) == 'S'` by "Seongil",
the check for `si == 0` preventing a divide-by-zero is missing compared to the
2D version. It has been readded.

### Missing monotonicity checks

The original search algorithms, as well as the new step-to-edge-case-aware
algorithms, both require that input data arrays are monotonically increasing in
relevant dimensions, and will produce incorrect results if this is violated.
BELLHOP did not check for this in altimetry X and Y, and bathymetry X and Y,
but this has been added.

### 3D SSP Hexahedral size overflow check

In 3D SSP Hexahedral, `SSP%Nz` is assigned to `SSP%NPts` and `SSP%Seg%z` to
`SSP%z`. However, the latter has a maximum of `MaxSSP` depth values, whereas
the former did not have a check for exceeding this size. A check has been added.


## Other issues

### Polarity flipping

In `InfluenceCervenyRayCen`, `rnV` was negated for image 2 and again for image 3.
If `Nimage` was set to 2, this means `rnV` would change sign every step. mbp
had implemented this differently (without the bug) in `InfluenceCervenyCart`,
so that implementation was brought to this function.

### Precomputing phase in 2D hat ray-centered influence

In `InfluenceGeoHatRayCen` (2D), the original code did `IncPhaseIfCaustic` (or
rather, the non-function equivalent of it) within the `Stepping` loop, after
checking and skipping the step if `irA == irB`. Which steps meet this condition
depends on `Pos%Rz`, which is receiver information and has nothing to do with
the ray/beam itself. If `q` goes from negative to positive on one step, the
caustic will go unnoticed until the next step where `irA != irB`. However, if
there are TWO caustics during this time (i.e. at two or more sequential steps
where `irA != irB`), `q` will have gone from negative to positive back to
negative and the phase change will be completely missed. This means field
results at one receiver depend on where other receivers are, which is
non-physical.

mbp had a comment suggesting this be precomputed, which is how it is done in
the corresponding 3D influence functions, and is the solution which has been
applied here.

### Curvature correction rotation divergence

In `CurvatureCorrection2` in `Reflect3DMod`, even if there is no curvature
correction (`R1`, `R2`, and `R3` are all zero), rotating `p` into and out of the
ray's coordinate system introduces slight floating-point error, which
accumulates over every step and can trigger edge case differences later. This
has been improved by keeping `p` constant when there is no curvature correction.

### Intentional asymmetry between X and Y crossing thresholds

In `ReduceStep3D` and `StepToBdry3D`, the threshold for crossing segments
in X was `EPSILON( h1 )` but in Y was `1000.0 * EPSILON( h1 )`. mbp had a
comment questioning this. The 1000 has been removed.

### Memory leaks of temporary arrays

`TopGlobalx`, `TopGlobaly`, `BotGlobalx`, and `BotGlobaly` were never
deallocated. This has been fixed.

### Other

There have been other minor fixes which should not affect results unless you're
doing something strange. See commit logs.


# Features which are probably bugs

Some features of `BELLHOP` / `BELLHOP3D` are probably not intended, but they
were not difficult to emulate in the C++/CUDA version, so their behavior was
left alone. Or, there is code which has no effect and was removed in the
C++/CUDA version.

### Three issues with phase and caustics

The code updating `phase` and computing the final value of `phaseInt` in the
various influence functions has been moved to the functions `IncPhaseIfCaustic`
and `FinalPhase` respectively. There are three differences in how this code
operates between versions, which do not seem to have any mathematical basis and
therefore may be bugs.

First, the step of the ray whose phase is read--`X` for `ray2D( X )%Phase` (or
the same for 3D)--is `iS - 1` in all cases, except for 2D Gaussian Cartesian
where it is `iS`. Given that the code then updates the phase based on whether
the caustic has yet occurred during the step, it seems that starting from the
initial phase of the step (rather than the final phase of the step) would be
correct. So, the 2D Gaussian Cartesian version may be a mistake.

Second, in all the 2D influence functions, but none of the 3D ones, if the
caustic has occurred, the ray phase is discarded, and only the accumulated
phase from previous caustics (plus the new pi/2) is used. This seems like it may
have been a typo which was copied-and-pasted to all the 2D functions, but then
fixed in the 3D functions.

Finally, when checking for caustics, the position of the equals sign in the
greater than / greater-than-or-equal-to / etc. symbols differs. In 2D SGB, the
equal signs are on `qOld`, whereas in all the other cases they are on `q`. Since
this is just detecting a change from positive to negative or vice versa, neither
choice is more correct than the other. But it's not clear why this should differ
between different influence functions.

### Strange behavior in `r` / `ir` initialization in 3D Cartesian

- `rA` is computed as the distance between the start of the ray and the source,
  but the ray always starts at the source so this is always zero.
- `ir` is then initialized to the minimum index of receivers at `r > rA`. This
  will be 1 if there are not receivers at 0, or 2 if there are.
- On the first step only, if the step is not a duplicate, there is a check
  for whether the "ray is moving outward in range", which it always is because
  it always starts from range 0 (`rA` is 0, `rB` is a norm so it must be `>= 0`,
  and if it was 0 it would be a duplicate step and not get there). This sets
  `ir` to 1.
- Therefore, if the first step is a duplicate, receivers at 0 will be skipped,
  but if the first step is not a duplicate, receivers at 0 will be included.

This is unlikely to be desired behavior. Just initialize `rA = 0.0` and
`ir = 1` and remove the other code.

### Dead code in 3D ray-centered

Both hat and Gaussian 3D ray-centered store `m` and `delta` from `B` to `A`
after every step, as if the last step's value was going to be reused. But the
last step's value is recomputed every step anyway, as is needed because steps
can be skipped, making the old `A` value stale.

Gaussian 3D ray-centered also initializes `deltaA`, `irA`, and `mA` before the
stepping loop, but their values are immediately overwritten. It also sets `irA`
to zero with a comment that this is a flag, but there is no code to check this
condition.

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
