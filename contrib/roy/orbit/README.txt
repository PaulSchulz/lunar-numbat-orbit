###
### OVERVIEW
###
### This program implements a simple-minded model of Numbat motion,
### in the vicinity of the Earth and Moon.
###
### Included are the effects of gravity from both the earth, moon and sun.
### An on-board engine can be fired to accelerate the craft.
###
### The position of the moon can be approximated using internal functions,
### or positional data of DE405-level accuracy can be used by contacting
### JPL's Horizons ephemeris service (the program will do this automatically
### if requested via the appropriate command-line option).
###
### The state vector of the Numbat is converted to orbital elements if the
### appropriate command-line option is given.
###
### The program now also includes selenographic position and velocity
### information, which may prove useful as the Numbat approaches the
### moon.
###
###
###
###
### --------------------------------------------------------------------------
### 
### DETAILS
###
### All this chaos takes place in something approximating a geocentric
### inertial frame, using equinox J2000.0.  We use two coordinate
### systems herein.  First, Ecliptic coordinates (longitude, latitude
### and distance).  Second, a right-handed, cartesian system defined
### such that the "positive x" axis is in the direction of Ecliptic
### longitude = 0, latitude = 0, and the "positive y" axis points
### towards longitude = pi / 2 rad, latitude 0.  The "positive z" axis
### is oriented such that positive Ecliptic latitudes correspond to z
### values > 0.
###
### Included are the effects of gravity from both the earth, moon and sun.
###
### Both the moon and sun appear in the code as point-like masses, so
### we're effectively including these bodies as perfect spheres of
### uniform mass density.  Exactly how good an approximation this is
### when we get close to the moon is unclear to me at the moment.  The
### Earth is modelled as an oblate spheroid; that is, a point-like
### mass (i.e., gravitational potential = -GM/r) but also including
### the zonal harmonic associatend with the "n = 2" Legendre
### polynomial.  As the magnitude of this n=2 term is approximately
### 1000 times that of the next largest term, only n=2 has been
### included here.
###
### I have chosen the classic, second-order integrator affectionately
### known as "leapfrog".  Google or read any book on numerical methods
### for more info on this.  The moon_internal(), sun() and anomaly()
### functions were adapted from the book "Astronomy with your Personal
### Computer" by P. Duffett-Smith Cambridge University Press.  The
### function libration() was adapted from "Practical astronomy with
### your calculator", also by Duffett-Smith.  The position of the Moon
### is accurate to maybe 50 km using moon_internal().  You can do
### better by specifying the "-horizons" option, as this will obtain
### data from the JPL Horizons web interface.  This gives the program
### access to moon position and libration data at a level of accuracy
### comparable to JPL's DE405.
###
### When requesting data from JPL's Horizons service, data are
### obtained at five-minute intervals, and cubic splines are used
### within this program to interpolate values between these points.
### Data from Horizons are saved in a disk cache file, so they can be
### reused in subsequent runs without having to contact Horizons
### again.
###
### 
###
###
###
### COMMAND LINE OPTIONS
###
### The "-horizons" option will ask the program to use Moon position data
### from the JPL Horizons web interface.  These data will be cached in a
### file on local disk.  This gives the program access to moon position
### data at "DE405 accuracy".
###
### The "-elements" option will output instantaneous orbital elements
### of Numbat if appropriate.  Maybe it's just because I'm an "orbit guy"
### but I find the ever-changing numbers of a state vector quite
### difficult to understand.  If in orbit about the Earth, this option
### will output instantaneous orbital elements at every reporting line.
### Two sets of elements are printed: one is refered to the Ecliptic
### plane and the other to the Earth's Equatorial plane.  If in orbit
### about the moon, elements are printed refered to the Ecliptic plane
### only at the moment.
###
### 
### 
### 
### Data reported during the run is as follows:
### 
### UNIX ctime        ctime in the simulated world.  This is useful to adjust
###                   engine timings, as they are set by ctime.
### 
### UTC               The time in the simulated world, UTC.
###
### Elapsed           Elapsed time in the simulated world,
###                   since the run was started.
###
### [xyz]             Cartesian coordinates of Numbat (m)
### 
### d[xyz]/dt         Velocity components of Numbat (m/s)
### 
### d2[xyz]/dt2       Acceleration components of Numbat (m/s**2)
### 
### Mass              Mass of the Numbat (kg).
###                   May change as the engine burns fuel.
### 
### Eng               Engine state.  May be either 'on' or 'off'.
### 
### RA and Dec        J2000.0 Right Ascension and Declination.
###                   These Equatorial coordinates may be directly
###                   plotted on the sky,
###                   e.g., using Google Earth.
###
### Earth Dist        Distance of Numbat from the centre of Earth.
###                   Reported both in metres and in Earth radii.
### 
### Moon  Dist        Distance of Numbat from the centre of the Moon.
###                   Reported both in metres and in Moon radii.
### 
### Speed wrt Earth   Speed of Numbat (m/s) with respect to the
###                   centre of the Earth.
### 
### Speed wrt EMB     Speed of Numbat (m/s) with respect to the
###                   Earth-Moon barycentre.  This is a dynamical
###                   point which is the centre of mass of the
###                   Earth-Moon system.
### 
### Speed wrt Moon    Speed of Numbat (m/s) with respect to the
###                   centre of the Moon.
###
### EnergyE           Total energy of the Numbat with respect to Earth.
###                   That is, kinetic energy with respect to the Earth
###                   summed with the potential energy in the Earth's
###                   gravitational potential.  If Numbat is near the
###                   Earth, this is a convenient measure of whether
###                   you're in orbit.  If < 0, you're in orbit.
###
### EnergyM           Total energy of the Numbat with respect to the Moon.
###                   That is, kinetic energy with respect to the Moon
###                   summed with the potential energy in the Moon's
###                   gravitational potential.  If Numbat is near the
###                   Moon, this is a convenient measure of whether
###                   you're in orbit.  For example, if you burn the engine
###                   to drop in to lunar orbit and this term remains > 0,
###                   you didn't do it right.  :)
### 
### Moon [xyz]        Cartesian coordinates of the Moon (m).
###                   These are useful for plotting purposes.
###
###   For example, to make a 3D plot of Numbat and the Moon,
###   you could do this
###
###   [roy@localhost ~]$ ./orbit >output.out
###   ^C
###   [roy@localhost ~]$ cat output.out  | grep -v Note | grep -v UTC \
###   | awk '{print $5 " " $6 " " $7 "\n" $27 " " $28 " " $29}' >orbit.xyz
###   [roy@localhost ~]$ gnuplot
###   gnuplot> set term post eps
###   gnuplot> set size ratio 1
###   gnuplot> set style data dots
###   gnuplot> set output "orbit.ps"
###   gnuplot> splot "orbit.xyz"
###   gnuplot> quit
###   [roy@localhost ~]$ 
###
###
### The following data are probably only of much use as the Numbat
### approaches the Moon.
###
### Selenographic     The selenographic longitude and latitude (degrees)
###                   of the "sub-Numbat" point.  That is, the lunar
###                   longitude and latitude of the point on the moon's
###                   surface which lies directly underneath the craft.
###
### Altitude          The distance (m) of the craft above the "nominal"
###                   surface of the moon.  More precisely, the radius
###                   of the moon (the variable $Rm) is used for this value.
###
### VerticalSp        Speed at which the craft is moving towards the lunar
###                   surface (m/s).
###                   Values > 0 mean that altitude is increasing.
###
### HorizSp           Speed of Numbat in the plane which lies tangential to
###                   the moon's surface at the "sub-spacecraft" point.
###                   Basically the "horizontal speed" referred to the
###                   surface of the moon.
###
### Angle             The position-angle of the horizontal velocity vector
###                   (degrees).
###                   A value of 0 means the craft is moving towards the
###                   north pole.
###                   A value of 90 degrees means the craft is moving due East
###                   (i.e., selenographic longitude increasing,
###                   latitude invariant).
###
###
###
### 
### 
### ORBITAL ELEMENT DATA
###
### If you're using the "-elements" option, each line of reported data
### will be accompanied by orbital data, including orbital elements.
### These are as follows:
### 
###      a = length of semi-major axis (m)
###      e = orbital eccentricity
###      i = orbital inclination (degrees)
###  Omega = longitude of the ascending node (degrees)
###     LP = longitude of periapsis (degrees)
###  omega = argument of periapsis (degrees)
###      M = mean anomaly (degrees)
###      T = true anomaly (degrees)
###     PD = periapsis distance (m)
###     AD = apoapsis distance (m)
###      P = orbital period (s)
### 
### 
### 
### 
### LICENSE
### 
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### Please see <http://www.gnu.org/licenses/> for more information.
### You know you want to.
### 
### 
### 
### HISTORY
###
### [17-Feb-2009 Roy Duncan]    Initial (poor) coding.
### [02-Mar-2009 Roy Duncan]    Added GNU GPL.
### [XX-Apr-2009 Roy Duncan]    Added differential solar gravitation
###                             Added orbital elements
###                             Added Earth's equatorial bulge.
###                             Added code to access JPL Horizons interface.
### [09-May-2009 Roy Duncan]    Fixed up documentation.
### [26-May-2009 Roy Duncan]    Added senenographic handling.
### [29-May-2009 Roy Duncan]    Added vertical and horizontal velocity data.
###                             Clean up (some of) the mess.
###
###
###
###
### VERSION
###
### This is version 2009-May-29.
###
###
