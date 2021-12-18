#Bellhop Beam tracing for ocean acoustics

   Copyright (C) 2009 Michael B. Porter

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



12/94
Added option for reading tabulated reflection coefficient.

 1/95
Added eigenray option which produces a ray plot but only for those rays which connect the source and receiver(s).

 3/95
Sign changed on pressure field to make Bellhop consistent with the other models.

 9/97
Option to propagate time series directly in Bellhop has been removed and replaced with the option of writing a arrivals file containing delays and amplitudes for each arrival. Phase changes have been pulled from the delay vector into a complex amplitude. Matlab script stackarr.m written to stack up the arrivals taking into account all phase changes.

11/97
The ray file is now written in ASCII rather than unformatted binary. This makes it easier for external programs to read ...

12/97
The option of including bottom bathymetry has been added.

 5/98
Small boo-boo fixed. The / option for the sound speed profile is supposed to cause a default to the last used value. That didn't work for the lower halfspace ...

11/98
If you have a very reflective bottom, rays may end up propagating back towards the source. Bellhop has been modified to allow it to search backwards through bathymetry points to properly treat bottom reflections in this case. Also, the code now searches through as many subsequent bottom sections as necessary to find the proper patch of bottom. Thus, if you have a finely sampled bottom the code will be sure and advance far enough to get over the right patch even if your ray step length is large.

3/00
The arrivals file, which gives arrival times and amplitudes, now also gives the ray-angle.

1/01
The code previously took for granted that the first bathymetry point was given at range, r=0.0, i.e. directly below the source. It has now been extended to start by searching the bathymetry table to find the segment below the source. If no such segment exists, Bellhop generates an error message and terminates.

2/01
A phase error was corrected in the treatment of reflection by a halfspace (due to a sign error in calculating the angle of incidence). This problem was not present when halfspaces were first implemented but crept in during a subsequent modification.

2/01
A user doing an 'arrivals' run for 5000 ranges and 32 depths found that the maximum number of arrivals allowed per rcvr was 0. The associated matrices have been increased in size; a threshold has been set to always make the matrices big enough for at least 10 arrivals; a facility has been added such that when there are too many arrivals, the program automatically keeps just the largest ones. Finally, tests have been added throughout where dynamic allocation is done, to generate an informative error message if the allocation fails due to insufficient (virtual) memory.

May 2001
The option for using a top reflection coefficient had an error. Basically it produced a negative angle of incidence so that the lookup of the reflection coefficient always returned the value at 0 degrees incidence. This did not effect bottom reflection coefficient tables.

May 2001
Some fairly sizeable changes ... hopefully with no new introductions of bugs. The handling of non-smooth SSP's has been changed so that the code automatically adjusts the step size during the ray tracing so that it steps precisely onto those discontinuities (of cz). It also makes the requisite curvature change in the dynamic ray equations. With these changes, the code is now safe to use for piecewise linear profiles (option 'C' or 'N' ). There is no effect in the spline option which automatically produces smooth profiles. I now suggest using the piecewise linear option for most cases, since if you don't take care with the spline option, it will wiggle significantly (even though it always passes through the given SSP points). Finally, I have modified the treatment of bottom reflections so that the code is more accurate for steeply sloping bottoms. The code used to handle the curvature corrections assuming the bottom was locally horizontal.

August 2001
Option added to output arrivals file in binary format.

December 2001
Correction to treatment of reflection coefficients. The original code treated the coefficient as being a function of the vertical angle. The current one treats it as a function of the angle of incidence (referenced to the normal to the bottom). This had previously caused an error in treating the surface reflection coefficient because the angles of incidence were negative.

The option of having a non-flat surface has been added. The new version optionally reads an AltImetrY file (ATIFIL) with that info. Also a bug was fixed that caused Bellhop to confuse the last SSP point with the bottom depth. These are the same when a bathymetry file is not used, but with a bathymetry file the use could specify a bathymetry point deeper that the last SSP point. (The bug only caused the code to run slower.)

Bellhop has been changed throughout to use vector quantities for the ray coordinate, ray tangent, and bathymetry. This makes the code much more concise.

May 2002
More complete ray history information has been added to the arrivals file 'ARRFIL'. The file now includes the ray angle at the source and receiver as well as the number of surface and bottom bounces.

A small modification has been made to ensure that the sound speed on the ray path is calculated at the terminal point of the ray. This info is used later in the TL calculation and causes a glitch for the last range. (Thanks to Lu Lian-Gang for bringing this to my attention.)

An option has been added to give attenuation in terms of a 'loss parameter'.

The arrivals file has been modified to write the amplitude and phase of the arrival rather than a complex amplitude. This facilitates interpolation of arrival data between environments in cases where there are many phase rolls between the environments. (The phase is 'unwrapped' and therefore keeps track of the accumulated phase including all 360 degree loops of the origin.)

June 2002
With sloping bottoms, Bellhop can generate rays with angles of incidence greater than 90 degrees. (Imagine a vertically launched ray travelling downward to a bottom sloping upwards. It can also happen with rays that reverse direction due to upslope angle steepening.) Typically people only provide reflection coefficients to 90 degrees. The new version of Bellhop generates a warning message when it needs a reflection coefficient outside the tabulated domain. It also assumes the reflection coefficient is symmetric about 90 degrees. Thus an angle of incidence of 95 degrees is treated like 85 degrees.

October 2002
Added an option to include a source beam pattern. This is done by reading an array of values from a separate source beam-pattern file, SBPFIL.

November 2002
Bellhop was generating NaN's when sources were placed in the sub-bottom. The code was never set up to handle this; however, in some applications it was convenient to provide for this. The new version does so by terminating the ray trace immediately and producing a vanishing pressure contribution for such cases.

December 2002
The Thorpe volume attenuation was not being applied in the case where an arrivals output file was selected (runtype 'A' or 'a'). Separately the rough surface option (altimetry file) was being invoked from option letter 4 rather than 5 as specified in the documentation.

January 2003
The option for automatically calculating the number of beams to use has been implemented. This was in the documentation but somehow had never been put in. The formula used is often going to be too conservative, i.e. it tends to use more beams than you really need.

April 2004
When the number of receiver ranges was reduced (to 1 in the test case), a bug was causing an error in the computed TL. The problem was traced to an assumption in the code that the receiver sampling would be fine enough to catch phase changes resulting from caustics. Fixed ...

December 2005
Added a line source option.

January 2005
Improved the beam interpolation. Shortcomings were obvious in using Bellhop for nearly vertical propagation. Also fixed a problem in which the source coordinate was being counted as a caustic leading to an extra caustic phase change. Added an option for 'irregular' (as opposed to 'rectilinear') receiver coordinates. The purpose was to allow receivers following the bottom terrain, which in turn is needed for reverberation calculations. Invoking the option causes Bellhop to interpret receiver depths and ranges as pairs rather than defining a rectilinear grid.

October 2007
Added an option for curvilinear interpolation of the boundaries. This provides a smoother interpolation than the piecewise-linear option that had been the default. Smooth interpolation is important when the geometric ray option is invoked as geometric ray theory is very sensitive to interpolation issues. This option also requires a change in the format of the altimetry and bathymetry files--- an option letter has been added that selects piecewise linear ('L') or curvilinear ('C'). The changes have been implemented in both the Fortran and Matlab versions of Bellhop.

November 2007
Bug in the treament of altimetry/bathymetry files: Bellhop is supposed to extend boundaries to infinity in a constant fashion, using the values at the endpoints. This was not working properly and caused a problem when the user traced rays to a range beyond where the boundaries were defined.

December 2007
The Bellhop integrator is a polygon method that uses an Euler first step, followed by a midpoint correction. (The polygon or midpoint method is a second-order member of the Runge-Kutta family of codes.) A change was made to have Bellhop recheck the allowed step size in the second or midpoint phase. This provides a more accurate treatment of the boundary reflection.

A new routine was added that checks all the limits on the ray step. The ray step is limited by top and bottom reflections, interface crossings (points where the SSP is sampled), and segment crossings (points where the bathymetry is sampled). With this addition, one does not need to worry about using a small range step in problems with fine sampling in the boundaries, since the code will carefully reduce the step size to land on such boundaries. This allows larger range steps, and therefore some increase speed.

A section of code has been added to the code for geometric and Gaussian beams that does a coarse pre-check of the limits of receivers that may feel the influence of a beam. This speeds things up considerably since certain beam calculations can then be omitted for receivers that have no chance of receiving a beam contribution.

An option was added for range-dependent sound speed profiles, which are specified on a rectangular grid and supplied in an extra SSP file. Beam curvature formulas have been updated to consider the gradients of the SSP in range.

The logic has been implemented to calculate the TL associated with back-travelling rays. Until now, Bellhop could trace rays that went backwards, but could not calculate the field due to such rays.

April 2008
Fixed logic that was supposed to detect the last SSP point by an approximate match to the value of Depth. The test was too stringent, and didn't always find a match because of roundoff.

April 2009
Subroutine for reading tabulated reflection coefficients was not closing files and deallocating variables. Caused problems in KRAKEN runs with multiple profiles. that used such tabulated reflection coefficients.

December 2009
A trap was added in Bellhop to catch rays that escape the domain of the waveguide. This can happen when the curvilinear option is used and then significant CPU time is used to trace a ray to the storage limit, which does not contribute to the field. Converted vectors defining the ray to a structure in an attempt to improve readability and reduce cache failures. Converted some remaining variables to double precision, imagining that single precision saves little on run time, may even slow things down if conversions are required or compiler optimizations are blocked, and frequently turns out to be inadequate precision.

January 2010
It was pointed out to me that the Quad SSP interpolation was doing linear extrapolation for values outside the grid, but that piecewise constant extrapolation would be better. My fix was supposed to limit the interpolation parameter, s, to the interval [ 0, 1 ] suppressing any extrapolation. However, min and max were switched so that s was always set to 0, effectively forcing piecewise constant interpolation/extrapolation in every segment of the grid. Fixed.

August 2010
There was a bug in the new, vectorized, Fortran 2003 version in the case of a single receiver range. (The code was looking through the full length of the receiver vector, which includes a couple of special values used to produce a linearly spaced vector of receiver ranges. Fixed.

June 2011
Two small bugs fixed: the boundary curvature was not being set for the last segment of the bottom or top bathymetry. A dummy variable for the halfspace depth was not being set if the user used '/' for the completely halfspace specification. (Thanks to J. Peterson for reporting this.)

January 2012
A problem (bug) occurred if you didn't have the sound speed start a zero depth, because Bellhop was calculating the depth of the top boundary using the minval of the vector of depths where the sound speed was sampled. That included the whole length of the dimensioned array, rather than just the part used to store the sound speed, so it was using a depth of zero for the top boundary. Actually some parts of the code had the zero depth, and some parts depended on the shallowest depth for the SSP, which is what caused problems. This has been fixed to always use the first depth of the SSP. For ocean acoustic problems, this normally implies that you should make sure you have an SSP point at zero depth, since the rays will reflect at the depth of that first SSP point. 

February 2013
A range-dependent bottom (or top) type is now allowed in Bellhop. The information is added in to the bathymetry (or altimetry) file. Major changes have been made throughout to use more structure variables. The code has also been re-structured to facilitate writing a Mex wrapper so that the Fortran Bellhop can be called from Matlab.

February 2013
When a source is placed on an SSP point, which is also a local minimum then a horizontally launched ray is a pathological situation. If its launch angle or position were slighly perturbed then it would generate a ray that cycles repeatedly in tiny loops. The logic for dynamically adjusting the step size along a ray has always been confused by this and Bellhop would take infinitessimally small steps until running out of storage space for the ray. The new version detects this pathological case and uses the standard stepsize. Hopefully this will not have other bad side effects.

March 2013
A term involving the rz derivative of the SSP had been neglected when the Quad option was used to interpolate a 2D SSP. Similarly, curvature changes in crossing range-segments in the SSP had been ignored. Fixed ... Logic was also added so that the ray stepsize would be automatically adjusted to land on the ranges defining the SSP. Previously we had required the user to control that through the bathymetry or altimetry file.

July 2013
Changes made to make it easier to call Bellhop as a subroutine from Matlab. In particular all the variables can easily be set from the main routine or by reading them from an envfil.
Fixed an error in the formula that handles changes in beam curvature during boundary reflection. The earlier formula implicity assumed the sound speed gradient was perpendicular to the boundary. This was also fixed in Bellhop3D for both Nx2D and 3D options. The results are clearly improved when there is a complicated bathymetry together with ocean SSP gradients.

February 2014
There is a further curvature correction to beams when they reflect from a curved boundary. (Boundaries in Bellhop can either be curved or piecewise linear.) The curvature correction had the wrong sign for the curvature.

May 2014
The arrivals file had inaccurate receive angles in the case when Gaussian beams were used. Bellhop was storing the receiver angle associated with the last beam contributing to a given receiver location. This was not very accurate. Bellhop has been modified to use the average angle of all the contributing beams (in a group of beams associated with a certain eigenray). Bellhop3D has also been modified accordingly.

April 2015
The logic in the routine that calculates the beam influence for Geometric Gaussian beams had been rearranged for greater efficiency. In the process a factor of 2 was dropped in calculating the intensity with the incoherent field option. Thus incoherent calculations with Geometric Gaussian beams were coming out 3 dB too low. This error affected the 2015 release of Bellhop.

June 2015
A user wanted to see ray traces for multiple source depths in the same plot. A minor modification was done to permit this.

November 2015
Bellhop had used an approximate treatement of ocean volume attenuation based directly on the cylindrical range to the receiver. The volume attenuation should really be based on the path length of the ray. This has been a standard approximation in ray models; however, it does not always work well. For instance, for HF acoustic modems in deep water being used over a vertical path, this approximation yields 0 loss (since the cylindrical range to the receiver is 0), thus ignoring significant loss in the depth direction. Furthermore, Bellhop had no provision for a depth-dependent volume attenuation. Both Bellhop and Bellhop3D have been modified to allow a depth-dependent attenuation and to incorporate losses in terms of the length of the ray path. The changes to do this were simple but numerous. All the routines for interpolating the soundspeed profile were modified to return a complex sound speed. The option to write an arrivals file was modified to write a complex delay where the imaginary part represents the attenuation along the eigenray path. Matlab outines that read the arrivals file were also modified. All the beam influence routines were modifed to incorporate the complex delay.  A new test case ('VolAtt') was added to exercise the new capability. With so many changes, I assume some bugs will show up, despite the extensive testing.

April 2016
The incoherent calculations for Geometric Gaussian beams were still about 1.6 dB too low. I had assumed loosely that two adjacent beams were contributing. I modified the formulat to use the sum of an infinite number of beams (with increasing distance form the receiver). This was actually the process that had always been used for the *coherent* sum, but I had been less careful with the incoherent sum. The updated formula now increases to a fracion of a dB with what had been in the Geometric Hat beams. Further refinements are possible for refracted media.

July 2016
Fixed an error in Bellhop3D using Geometric Hat beams in Cartesian coordinates: the old version would generally have a beam make a contribution to a radial and its continuation in the other half plane. For instance, if the beam made a contribution to the 45-degree radial, it would also show an erroneus contribution to the 225-degree radial.

August 2016
Added an option for geometric hat beams in ray-centered coordinates (rather than Cartesian coordinates). The geometric beams in Cartesian coordinates have been a standard. Only the Cerveny-style beams had been implemented in ray-centered coordinates. The interpolation of the beams works better with ray-centered coordinates, especialy for beams that are nearly vertical. However, the Cartesian beams require less computation.

September 2016
Fixed a bug in Bellhop3D in which the 2D option was not clearing the pressure field for a radial. As a result, the field on new radials had the field in storage from the previous radial in areas that would have been shadows.

Modified the Gaussian beams in ray-centered coordinates to use the same algorithm as that used for hat beams in ray-centered coordinates. The new code is faster by about 4x. (That factor depends on lots of particulars of the problem and the set up.)

Bellhop reads the SSP data to a user-specified depth. That depth was being converted to single precision and then in a rare case did not produce a match at the last SSP point due to the lower precision. Fixed.

Added (but did not test) curvilinear altimetry option based on the analogous bathymetry option.

For Geometric, (Gaussian and hat beams) in ray-centered coordinates the 'Arrivals' option was incorrectly using a ray matrix from the 2D case (ray2D vs. ray3D). Fixed.

October 2016
The Bellhop3D altimetry option had probably never been tested. I created a directory 'ParaBot' with a parabolic altimetry and bathymetry file based on the same case in the Bellhop tests. It turned out this had never worked propertly. The logic is, of course, almost identical to the bathymetry case but it wasn't quite right. Based on the ParaBot test it looks like it is now working correctly. The pressure field is symmetric about the horizontal axis, indicating the altimetry and bathymetry files are being treated the same way. the match with the normal Bellhop is not quite perfect and deserves another look. However, it is awfully close.

Bellhop and Bellhop3D have always ignored the density in the volume assuming it was seawater and could be treated as a constant. That constant was taken as 1.0 gm/cm^2, which is very slighly off. In running an aeroacoustic problem for the air that density is about 1000x too high which led to an incorrect boundary reflection coefficient. Both models have been modified to use the density as specified in the environmental file including its depth variation. This effects about 6 different types of interpolation options in both models; however, I have just tested the c-linear option and only in Bellhop.

July 2017
The logic in Bellhop3D had always been sloppy in terms of how it would handle a ray that exited the domain where the bathymetry or altimetry was defined. This could result in accesses to undefined variables. (Bellhop handles this by extending the domain as a constant.) Bellhop3D has been modified to terminate the ray trace for the particular take-off angle involved when this happens. The ray trace then continues to the next launch angle.

Often problems in the ray trace occur because of user errors in defining the bathymetry. For instance, NaNs are often present. Bellhop3D now checks for such NaNs when they're first read and generates a warning message. It also tries to detect when a ray has entered a patch where the bathymetry is improperly defined (because of the Nan) and then terminates the trace of that particular ray.

A problem was also present in that Bellhop3D did not always check for each new ray which bathymetry/altimetry segment the source was over/under. This has been fixed; however, I need to review this since Bellhop uses a different process where it *always* searches for the active segment. This is less efficient but safer. In any case, Bellhop and Bellhop3D should likely be set up the same way.

The Bellhop3D ascii arrivals file has been modified so that it is stylistically consistent with the Bellhop version. (Previously the delay information was written as a complex number by Bellhop3D vs. two reals in Bellhop.)

December 2017
The arrivals option for the case of an Nx2D run had not been implemented. This feature was added so that an arrivals calculation can be done with either Nx2D or 3D options.

May 2018
The geometric Gaussian beam options in Bellhop and Bellhop3D use a 'stent' which sets a lower limit on how narrow a beam can be. That lower limit is based on arclength from the source so that it allows narrow beams close to the source. Otherwise the beams add up to too high a level. Bellhop and Bellhop3D had been using just a coarse estimate of arclength based on the number of steps times the stepsize. However, that estimate was way off in cases where the default stepsize is chosen to be really large. That is normally OK because the codes adjust the stepsize down automatically as the ray steps along. However, the number of steps times the stepsize then overestimates the arclength in a big way. Both codes have been modified to estimate the arclength using the travel time delay, which is tracked as part of the ray trace and incorporates the dynamic changes in steplength.

The semi-coherent option in Bellhop was over-riding any source beam pattern if that was provided. Fixed.

The Gaussian beam option in Bellhop3D (in 3D mode) now can be done in the original ray-centered coordinates or in Cartesian coordinates. The latter is vastly faster.

June 2018
The treatment of volume attenuation in the Matlab version of Bellhop has been upgraded to be consistent with the Fortran version. (See note dated November 2015.)
The fix of Feb. 2013 was also incorporated in BellhopM

July 2018
The Matlab version of Bellhop (BellhopM) has the wrong sign on the phase changes that happen when a ray crosses a caustic (because I used ' instead of .' on the associated vector, inadvertently conjugating it). Fixed. The arclength correction from May 2018 was also implemented in BellhopM.

A problem showed up in Bellhop3D when rays were traced over a bathymetry given in Eastings and Northings. Those values are in the millions, but described (in this case) a small patch. The ray steps are around 1 m and when written in single precion to the ray file, the ray coordinates in Easting and Northings don't have enough precision. The rays then look jittery. The fix was simply to write the ray file in double precision. This was changed in the regulard 2D Bellhop as well, although it's less likely one would create a Bellhop file with that kind of problem.

There were lines in WriteArrivals for the 3D case that factored in 1 / sqrt( r ) for a point source and something else for a line source. Neither one makes sense in a 3D run so they were removed. The 3D code does the equivalent of a sqrt( r) through the actual spreading of the beams in each direction; a line source would require a different treatment because of the 3D.

The format for ARRFiles generated by Bellhop3D has been modified to write the number of bearings (as well as receiver depths and ranges). Azimuthal angles at the source and receiver were also added to the ARRFile.

When Bellhop3D was run in the Nx2D mode it was neglecting to make the normals to the boundaries have unit length. (They were projections of the 3D normal vectors onto the 2D propagation plane; that process changed the length of the normal vectors.) As a result the angles of reflection were slightly off producing errors in the results.

August 2018
The previous version created pyramidal beams with a rectangular base. Thus rotations of the q_tilde, q_hat coordinates were ignored. The new version takes these into account yielding pyramidal beams with a parallelogram base. The resulting beams then produce what could be called a tiling, with essentially no gaps or overlaps between adjacent beams. A lot of work was also done to handle boundary curvature. When the curvilinear boundary option is selected it now automatically calculates the boundary curvature in each direction and applies jumps to the q functions (related to beam curvature) when the beams reflect off those boundaries. The previous version only incorporated curvature in calculating the angle of the reflected rays.

November 2018
The option for a bottom type that varied with range was not working properly. A confusion in the order of the real and imaginary parts of the P-wave and S-wave speeds was causing the sound speed attenuation to be ignored. An additional bug had been introduced when the broadband option was introduced into the acoustic toolbox-- the routine for converting attenuation to a complex sound speed was missing a newly introduced variable. Fixed ...

December 2018
The density interpolation with the hexahedral option was not interpolating density correctly. Since a varying density is rarely used in ocean acoustics this problem would likely not have shown up in any practical case. A modification has been made so that it is now interpolated ln linearly in depth. However, the interpolation is piecewise constant in x and y.

January 2019
PCHIP, Spline, and N2-linear SSP options added to the Matlab version of BELLHOP.

February 2019
2D or 3D soundspeeds are read in from an auxilliary SSPFILE. In that case, a depth-dependent volume attenuation can be input in the envfil. However, the attenuation was being interpolated incorrectly due to a confusion between depth into a layer and the proportional depth in the layer. This has been fixed but not tested. Note that the depths for the SSP values must match between the ENVFILE and the SSPFILE. There is currently no check to enforce that.

The Francois-Garrison formulas for volume attenuation have been added as an alternative to the older Thorp formula. In the process of doing this I found that I had volume attenuation incorporated twice in the Matlab version of BELLHOP. This error was introduced in June 2018 (see above) when I switched from handling attenuation by a simple formula based on cylindrical range to a more accurate one based on attenuation along the arc of the ray path. I failed to remove the older range formula when I added the new one.

May 2019
The incoherent option with paraxial beams (Cerveny) was producing levels that were way too low. The code was summing the absolute values of the pressures, rather than the intensities (pressure squared). With that corrected for both beams in Cartesian coordinates and ray-centered beams, the levels are better, but they still seem a bit off. However, the idea of incoherently summing beams that are in the same family does not really make sense.

In the case of a full 360 degree sweep of rays, BELLHOP and BELLHOP3D were supposed to automatically remove a duplicate beam at, e.g. -180 and 180 degrees. The coding was done lazily, testing for exact equality and with round-off that test was not always satisified. The code has been modified to test for approximate equality. This modification is of little importance, but it does eliminate an artifact in at/tests/BeamPattern/omni.

BELLHOP3D was not including the volume attenuation when the incoherent/semi-coherent option was selected.

June 2019
Elasticity in the halfspaces has been included. Note that BOUNCE can be used to generate reflection coefficients for elastic halfspaces and even more complicated layered sub-bottoms. Then the resulting reflection coefficient can be used in BELLHOP. However, this modification simplifies the process for halfspaces.

July 2020
Thre is a hidden option in BELLHOP to include a beam shift on boundary reflection. (It is hidden or undocumented because it is not considered to be reliable.) That option was being ignored when a ray trace run was selected. As a result, the ray plots did not reflect the true trajectories that were using in a TL calculation. Fixed ... This beamshift option is available in both BELLHOP and BELLHOP3D, but will not produce reliable results with complicated bottoms. There is no logic in the software to deal with the fact that a displaced ray/beam may be moved to the next bathymetric domain. However, this change at least ensures the ray plots will represent the underlying behavior for the TL.

August 2020
In the test case at/tests/BeamPattern/omni.env there are homogeneous lossless halfspaces for the top and bottom. Calculating the boundary impedance uses a square root to get a vertical wavenumber. In this special case it it taking the square root of a negative real. It turns out the gfortran and Intel compilers do different things with those negative reals (which lie on the the branch cut of the square root). They both give a purely imaginary value, but with different signs. The Fortran standard does not lock this down and the behavior can also be a function of the machine architecture. The effect was to generate a divide by zero with the Intel compiler. This was trapped with all the traps enabled. Logic has been added to force the desired sign and fix this inconsistency.

This error could occur with elastic layers and/or with other codes in the Acoustics Toolbox. However, one would rarely use a lossless halspace in those models. Also, SCOOTER, for instance always works a little off the real axis and wouldn't experience this. This also might be an issue in  the Matlab version of BELLHOP depending on the architecture.

September 2020
An error was showing up in the Fortran BELLHOP under gfortran 9.2, but not earlier versions of gfortran and not the Intel Fortran. The error caused the ray trace to take a large number of small steps. I found that the stepsize, h, had been declared as INTENT( OUT ), when it should have been INTENT( INOUT ) since a maximum step was being passed, that was supposed to be reduced by ReduceStep. So the step size is both an input value and an output value. Fixed. This error also was present in BELLHOP3D although there had been no reports of problems. It has been fixed there also.

October 2020
When the number of beams selected was 1, BELLHOP and BELLHOP3D interpreted this as a full 360 degree sweep (since the difference between the first and last beam MOD( 360 ) is 0. In such cases it decrements the number of beams resulting in 0 beams. Using just one declination angle would be a rare choice in BELLHOP; however, in BELLHOP3D it can make sense to trace just one bearing angle for an Nx2D run. In addition, a small change was made to ensure that the vectors for these angles are always allocated to a minimum length of 3. This is necessary because of the coding to automatically generate linearly interpolated values, which is flagged in the third element.
