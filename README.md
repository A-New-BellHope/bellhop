# BELLHOP
A mirror of the original FORTRAN BELLHOP underwater acoustics simulator, with
some bugs fixed.

# Impressum

Copyright (C) 2021 A-New-BellHope (Jules Jaffe team)
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

# Bugs fixed

### iSegz, iSegr initialization

In BELLHOP, `iSegz` and `iSegr` are initialized (to 1) only once globally. This
means their initial state for each ray is their final state from the previous
ray. *Normally, this shouldn't matter, because the segment indices are
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

The SGB (Bucker's Simple Gaussian Beams) Total Level influence algorithm has a
bug, acknowledged by mbp in code comments, where every step is assumed to be of
the full, maximum step size, which is often 500 m or larger. Thus, one extra
infinitesimally small step is actually one extra full-size step, which
substantially impacts the field results from this beam. Again, this change
happens "randomly" based on the final state of the previous ray/beam traced.

`iSegz` and `iSegr` have been re-initialized to 1 before every ray/beam is
traced, removing this "random" inconsistent behavior. The SGB bug still stands,
but at least it produces reproducible results now.


# Other information

See index.htm for further information.

Code retrieved 12/17/21 from http://oalib.hlsresearch.com/AcousticsToolbox/
claimed to have been updated 11/4/20.

Files pertaining to the other simulators (Krakel, Kraken, KrakenField, Scooter)
have been removed.

The Makefile has been set up to build with gfortran by default (Linux or
Windows with mingw etc.)
