#!/usr/bin/perl -w

use strict;
use Math::Trig;
use POSIX qw(strftime floor);
use Getopt::Long;
use Math::Spline qw(spline binsearch);
use Math::Derivative qw(Derivative2);
use LWP::UserAgent;
use Date::Parse;

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

###############################################################################
######################## Constants and assumed values #########################
###############################################################################
# Universal Gravitational Constant, SI units
my $G = 6.672E-11;

# Mass of Earth (kg)
my $Me = 5.9742E+24;

# Mass of Moon (kg)
my $Mm = 7.36E+22;

# Mass of Sun (kg)
my $Ms = 1.989E+30;

# Equatorial radius of the Earth (m)
my $Re = 6.378136E+06;

# Radius of the Moon (m)
my $Rm = 1.7374E+06;

# Speed of light (m/s)
my $c = 299792458;

# Earth axial tilt (rad)
my $earth_axial_tilt = deg2rad(23.43928);

# Moon axial tilt (rad)
my $moon_axial_tilt = deg2rad(1.543);

# Precession stuff: J2000.0 epoch in ctime
my $J2k_epoch = 946727936;

# Precession stuff: rate of movement of zero point of ecliptic
# longitude (rad/sec)
my $luni_solar_precession_rate = 7.68156698238282e-12;





###############################################################################
########################## Parameters for this run  ###########################
###############################################################################
# Start simulation at this ctime
my $initial_ctime = 1235259090;

# Time step, seconds.  If you make this bigger, be careful!
my $dt = 1;

# Vehicle initial mass, kg.
my $initial_mass = 400;

# Reporting step.  Print a line of status data every $report_every
# seconds of "simulated time".
my $report_every = 300;

# Print a data header every $header_every lines of output.
# Set to 0 to disable, but you'll still get a header as the first
# line.
my $header_every = 20;

# When engine running, you get this many Newtons of force.
my $engine_force = 10000;

# When engine running, vehicle loses this mass of propellant at this
# rate (kg per second).
my $mass_loss = 0.002;

# Fuel mass (kg).  This is so we know when or if we've run out.
my $fuel_mass = 1;

# Initial total energy of the vehicle wrt Earth (Joules)
# E for a circular orbit of radius r is -GmM/2r or in the parlance of
# this script -0.5 * $G * $Me * $initial_mass / $initial_radius
my $initial_energy = -7.5E+09;

# Initial distance from Earth Centre of Mass (m).  cf mean Earth radius of 6.3781E+06 m.
my $initial_radius = 6.628E+06;

# Call it quits if distance from Earth exceeds this value (m)
my $max_earth_distance = 1.0E+09;

# Call it quits if distance from Earth falls below this value (m).
# This is approximately one earth radius + 50 km.
my $min_earth_distance = 6.428E+06;

# Okay, let's set up some contrived initial conditions.

# Orbiting in the Ecliptic plane, initially setting y = 0 and dy such
# that values of y will be initially positive.  This is the same
# "sense" in which the moon orbits.  We might want to drop all of this
# stuff and instead specify the initial conditions as orbital
# elements.  First, calculate initial speed (m/s), appropriate to set
# vehicle energy as specified.

my $initial_speed = 
    sqrt(2 * 
	 (($initial_energy / $initial_mass) + ($G * $Me / $initial_radius)));
my $numbat = { # Radius vector
               'x'  => $initial_radius,
               'y'  => 0,
               'z'  => 0,
               # Velocity vector
               'dx' => 0,
               'dy' => $initial_speed,
               'dz' => 0,
               # Vehicle mass.
               'm'   => $initial_mass,
               # Engine state: 
	       #   1 = running, applying force parallel to the direction
               #       of Numbat motion (wrt Earth),
               #   0 = not running,
               #  -1 = running, applying force antiparallel to the
               #       direction of Numbat motion (wrt Earth).
               'engine' => 0
             };
               

### Engine start and stop times.  This is a set of ctimes for engine
### state changes, and is kinda clunky.  A '1' or a '-1' at the end
### means turn on the engine at the specified ctime, and a '0' at the
### end means turn the engine off at the specified ctime.  The force
### of the engine is directed parallel to the Numbat's direction of
### motion, with respect to Earth, if '1' is used.  This will act to
### increase the speed of the craft.  If a '-1' is used, the force is
### directed antiparallel to the direction of motion, and acts to
### decrease speed.  This needs to be more flexible.  For instance,
### accelerating the Numbat in the direction of motion wrt Earth isn't
### terribly useful when shaping a lunar orbit.  Because the moon is
### in motion wrt Earth means that there would always be a component
### of the force directed tangential to the moon's orbit, which isn't
### terribly useful.

my @engine_states = ( # Burn to enter lunar orbit: 155 seconds
                      #'1235355005|0',
                      #'1235354850|-1',
                      # Burn to leave earth orbit: 180 seconds
                      '1235275860|0',
                      '1235275680|1'
                    );
               


###############################################################################
# Here we initialise some Horizons-related stuff and
# process command line options.
# First, the cache file we use to store Horizons data we've obtained.
my $horizons_moondata_file = "./orbit_moondata_cache.txt";

# Initialise data structures for Horizons data.
my %horizons = ();
my %spline = ();
my %months = ( 'Jan' => '01',
               'Feb' => '02',
               'Mar' => '03',
               'Apr' => '04',
               'May' => '05',
               'Jun' => '06',
               'Jul' => '07',
               'Aug' => '08',
               'Sep' => '09',
               'Oct' => '10',
               'Nov' => '11',
               'Dec' => '12');

my $print_orbit_element_data = 0;
my $use_jpl_horizons_data = 0;

GetOptions( 'elements'=>\$print_orbit_element_data,
            'horizons'=>\$use_jpl_horizons_data);
 
my $ua = 
    ($use_jpl_horizons_data) ? LWP::UserAgent->new(env_proxy => 0, keep_alive => 0, timeout => 20) : 0;



###############################################################################
#############################  Subroutines  ###################################
###############################################################################

# Print output column header
sub print_header {

    print "UNIX ctime       UTC                     Elapsed      x             y             z               dx/dt         dy/dt         dz/dt           d2x/dt2       d2y/dt2       d2z/dt2          Mass      Eng    RA       Dec        Earth Distance        Moon Distance           Speed wrt Earth    Speed wrt EMB      Speed wrt Moon   EarthE    MoonE       Moon x        Moon y        Moon z          Selenographic   Altitude      VerticalSp    HorizSp      Heading\n";
    return(0);
}

###############################################################################
# This subroutine reads moon positional data from disk.
# These data are sourced from the Horizons service and are
# stored on disk for when the model is run again.
#
sub read_moon_data_cache {

    my $error = 0;
    my $line_count = 0;

    # Init data structures if needed
    unless (defined $horizons{'moon'}->{'timestamps'}) {
        my @horizons_moondata_timestamps;
        $horizons{'moon'}->{'timestamps'} = \@horizons_moondata_timestamps;
    }
    unless (defined $horizons{'moon'}->{'longitudes'}) {
        my @horizons_moondata_longitudes;
        $horizons{'moon'}->{'longitudes'} = \@horizons_moondata_longitudes;
    }
    unless (defined $horizons{'moon'}->{'latitudes'}) {
        my @horizons_moondata_latitudes;
        $horizons{'moon'}->{'latitudes'} = \@horizons_moondata_latitudes;
    }
    unless (defined $horizons{'moon'}->{'distances'}) {
        my @horizons_moondata_distances;
        $horizons{'moon'}->{'distances'} = \@horizons_moondata_distances;
    }
    unless (defined $horizons{'moon'}->{'sublong'}) {
        my @horizons_moondata_sublong;
        $horizons{'moon'}->{'sublong'} = \@horizons_moondata_sublong;
    }
    unless (defined $horizons{'moon'}->{'sublat'}) {
        my @horizons_moondata_sublat;
        $horizons{'moon'}->{'sublat'} = \@horizons_moondata_sublat;
    }

    # Read in
    if (open(DATA, "<$horizons_moondata_file")) {
        while(<DATA>) {
            chomp;
            my @bits = split(/\|/, $_);
            if ($#bits == 5) {
                push(@{$horizons{'moon'}->{'timestamps'}}, $bits[0]);
                push(@{$horizons{'moon'}->{'longitudes'}}, $bits[1]);
                push(@{$horizons{'moon'}->{'latitudes'}},  $bits[2]);
                push(@{$horizons{'moon'}->{'distances'}},  $bits[3]);
                push(@{$horizons{'moon'}->{'sublong'}},    $bits[4]);
                push(@{$horizons{'moon'}->{'sublat'}},     $bits[5]);
                $line_count++;
            } else {
                $error = $line_count unless ($error);
            }
        }
        close(DATA);
    }
    if ($error) {
        print "An error was encountered in reading the cache file, named \"$horizons_moondata_file\"\n";
        print "The error was encountered on line $error.\n";
        return(1);
    }
    print "A total of $line_count lines were read from the cache file.\n";
    return(0);
}
###############################################################################
# This subroutine makes sure that the Horizons data held in
# memory is in good order.  For example, 300s between array entries.
sub clean_moon_data_cache {

    my @entries;
    for (my $i = 0; $i <= $#{$horizons{'moon'}->{'timestamps'}}; $i++) {
        push(@entries, sprintf("%10d|%f|%f|%f|%f|%f", $horizons{'moon'}->{'timestamps'}[$i], $horizons{'moon'}->{'longitudes'}[$i], $horizons{'moon'}->{'latitudes'}[$i], $horizons{'moon'}->{'distances'}[$i], $horizons{'moon'}->{'sublong'}[$i], $horizons{'moon'}->{'sublat'}[$i]));
    }
    @{$horizons{'moon'}->{'timestamps'}} = ();
    @{$horizons{'moon'}->{'longitudes'}} = ();
    @{$horizons{'moon'}->{'latitudes'}}  = ();
    @{$horizons{'moon'}->{'distances'}}  = ();
    @{$horizons{'moon'}->{'sublong'}}  = ();
    @{$horizons{'moon'}->{'sublat'}}  = ();
    $horizons{'moon'}->{'latest_ctime'} = 0;
    foreach my $entry (sort @entries) {
        my @bits = split(/\|/, $entry);
        # Only 5-minute mark entries, thanks
        next unless (abs(($bits[0] / 300) - int($bits[0] / 300)) < 1.0E-10);
        next unless ($bits[0] > $horizons{'moon'}->{'latest_ctime'});
        push(@{$horizons{'moon'}->{'timestamps'}}, $bits[0]);
        push(@{$horizons{'moon'}->{'longitudes'}}, $bits[1]);
        push(@{$horizons{'moon'}->{'latitudes'}},  $bits[2]);
        push(@{$horizons{'moon'}->{'distances'}},  $bits[3]);
        # Make sure sublong values oscillate around 0; i.e., no 360 -> 0 discontinuities.
        $bits[4] -= 360 if ($bits[4] > 180);
        push(@{$horizons{'moon'}->{'sublong'}},    $bits[4]);
        push(@{$horizons{'moon'}->{'sublat'}},     $bits[5]);
        $horizons{'moon'}->{'latest_ctime'} = $bits[0];
    }
    return(0);
}
###############################################################################
# This subroutine writes moon positional data to disk.
# These data are sourced from the Horizons service and are
# stored on disk for when the model is run again.
sub write_moon_data_cache {

    my $error = 0;
    my $line_count = 0;
    unlink "$horizons_moondata_file.new" if (-f "$horizons_moondata_file.new");
    if (open(DATA, ">$horizons_moondata_file.new")) {
        for (my $i = 0; $i <= $#{$horizons{'moon'}->{'timestamps'}}; $i++) {
            print DATA $horizons{'moon'}->{'timestamps'}[$i], "|", $horizons{'moon'}->{'longitudes'}[$i], "|", $horizons{'moon'}->{'latitudes'}[$i], "|", $horizons{'moon'}->{'distances'}[$i], "|", $horizons{'moon'}->{'sublong'}[$i], "|", $horizons{'moon'}->{'sublat'}[$i], "\n";
            $line_count++;
        }
        close(DATA);
    } else {
        $error = 1;
    }
    unless ($error) {
        rename "$horizons_moondata_file.new", "$horizons_moondata_file"
	    if (-f "$horizons_moondata_file.new");
    }
    if ($error) {
        print "An error was encountered in writing a new cache file, ",
	    "named \"$horizons_moondata_file\"\n";
        return(1);
    }
    print "A total of $line_count lines were written to the cache file.\n";
    return(0);
}
###############################################################################
# This subroutine contacts JPL's Horizons system and retrieves data
# for moon position.  This could benefit from a rewrite with WWW::Mechanize.
# Note that we obtain data for more than just moon position, but that is
# not currently utilised.
#
sub get_horizons_moon_data {

    my ($ctime) = @_;
    # Form start time string for Horizons
    my @bits = gmtime($ctime);
    my $start_time =
	sprintf("%4d-%02d-%02d", $bits[5] + 1900, $bits[4] + 1, $bits[3]);
    # Form stop time string for Horizons - we get 72 hrs at a time
    @bits = gmtime($ctime + 3 * 86400);
    my $stop_time =
	sprintf("%4d-%02d-%02d", $bits[5] + 1900, $bits[4] + 1, $bits[3]);

    my $uri = "http://ssd.jpl.nasa.gov/horizons.cgi";
    my $pause = 1;

    print "Contacting JPL Horizons service...\n";
    my $response = $ua->get("$uri");
    my $cookie_string = '';
    my $error = 0;
    if ($response->is_success) {
        my $content = $response->content;

        # Check out cookies
        foreach my $line (split(/\n/, $response->headers->as_string)) {
            next unless ($line =~ m/Set-Cookie/i);
            (undef, $cookie_string, undef) = split(/[:;]/, $line, 3);
            $cookie_string =~ s/^ +//g;
        }
    } else {
        $error++;
        print $response->status_line, "\n";
    }
    $error++ unless ($cookie_string =~ m/CGI/);

    # Set target body to Moon
    # Use "MB:10" for the Sun...
    unless ($error) {
        sleep $pause;
        print "Setting Horizons target body to Moon.\n";
        $response = $ua->get("$uri?body=MB:301&select_body=1", Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Set ephemeris type to OBSERVER
    unless ($error) {
        sleep $pause;
        print "Setting Horizons ephemeris type.\n";
        $response =
	    $ua->get("$uri?table_type=OBSERVER&set_table_type=1",
		     Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Set origin as geocentric
    unless ($error) {
        sleep $pause;
        print "Setting Horizons origin to geocentric.\n";
        $response =
	  $ua->get("$uri?s_lookup=Search&l_str=500",
		   Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Set time span and reporting interval appropriately
    unless ($error) {
        sleep $pause;
        print "Setting Horizons report time span.\n";
        $response =
	  $ua->get("$uri?set_time_span=1&start_time=$start_time&stop_time=$stop_time&step_size=5&interval_mode=m", Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Set output type as TEXT
    unless ($error) {
        sleep $pause;
        print "Setting Horizons data output type.\n";
        $response =
	  $ua->get("$uri?set_display=1&display=TEXT",
		   Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Set required quantities and units
    unless ($error) {
        sleep $pause;
        print "Setting Horizons required quantities.\n";
        $response = $ua->get("$uri?oq_14=1&oq_17=1&oq_20=1&oq_31=1&set_table=1&time_fmt=CAL&time_digits=SECONDS&ang_format=DEG&output_units=KM-S&range_units=KM&apparent=AIRLESS&airmass=&elev_max=&solar_elong_min=0&solar_elong_max=180&rts_flag=NO&ref_system=J2000&set_table_settings=1", Cookie => $cookie_string);
        $error++ unless ($response->is_success);
    }

    # Report the results of these changes
    unless ($error) {
        my $content = $response->content;
        foreach my $line (split(/\n/, $content)) {
            # Check ephemeris type is correctly set
            if (($line =~ m/s_type=1/) && ($line !~ m/OBSERVER/)) {
                print "Could not set ephemeris type.\n";
                $error++;
            }
            # Check target body is correctly set
            if (($line =~ m/s_body=1/) && ($line !~ m/Moon.*\[301\]/)) {
                print "Could not set target body.\n";
                $error++;
            }
            # Check observer location is correctly set
            if (($line =~ m/s_loc=1/) && ($line !~ m/Geocentric.*\[500\]/)) {
                print "Could not set observer location.\n";
                $error++;
            }
            # Check time span is correctly set
            if (($line =~ m/s_time=1/) && ($line !~ m/$start_time.*$stop_time/)) {
                print "Could not set ephemeris time span.\n";
                $error++;
            }
            # Check quantities are correctly set
            if (($line =~ m/s_tset=1/)
		&& ($line !~ m/QUANTITIES.*SECONDS.*DEG.*KM/)) {
                print "Could not set required ephemeris data.\n";
                $error++;
            }
            # Check output format is correctly set
            if (($line =~ m/s_disp=1/) && ($line !~ m/plain text/)) {
                print "Could not set required output format.\n";
                $error++;
            }
        }
    }
    if ($error) {
        print "Weird!\n";
        exit(1);
    }
    # GO!
    sleep $pause;
    print "Retrieving ephemeris now.\n";
    $response = $ua->get("$uri?go=1", Cookie => $cookie_string);
    $error++ unless ($response->is_success);

    # Report the ephemeris
    unless ($error) {
        my $content = $response->content;
        my $print = 0;
        foreach my $line (split(/\n/, $content)) {
            $print = 0 if ($line =~ m/^\$\$EOE$/);;
            if ($print) {
                $line =~ s/^ +| +$//g;
                my @data = split(/ +/, $line);
                my @date_bits = split(/\-/, $data[0]);
                my $timestring = sprintf("%4d-%02d-%02dT%8s +0000", $date_bits[0], $months{$date_bits[1]}, $date_bits[2], $data[1]);
                my $data_timestamp = str2time($timestring);
                my $obsrv_long = $data[2];
                my $obsrv_lat  = $data[3];
                my $np_ang     = $data[4];
                my $np_dist    = $data[5];
                my $delta      = $data[6] * 1000;
                my $deldot     = $data[7] * 1000;
                my $ec_long    = $data[8];
                my $ec_lat     = $data[9];
                push(@{$horizons{'moon'}->{'timestamps'}}, $data_timestamp);
                push(@{$horizons{'moon'}->{'longitudes'}}, $ec_long);
                push(@{$horizons{'moon'}->{'latitudes'}},  $ec_lat);
                push(@{$horizons{'moon'}->{'distances'}},  $delta);
                push(@{$horizons{'moon'}->{'sublong'}},    $obsrv_long);
                push(@{$horizons{'moon'}->{'sublat'}},     $obsrv_lat);
            }
            $print = 1 if ($line =~ m/^\$\$SOE$/);;
        }
    }
}
###############################################################################
# Gateway for moon position calculations. 
# Use internal approximations or DE405-level accuracy (external) ?
sub moon {

    my ($ctime, $retardation_time) = @_;

    unless ($use_jpl_horizons_data) {
        # Use internal approximation to lunar position
        my ($longitude, $latitude, $distance) = moon_internal($ctime, $retardation_time);
        my ($sublong, $sublat) = libration($ctime, $longitude, $latitude);
        return($longitude, $latitude, $distance, $sublong, $sublat);
    } else {
        my ($longitude, $latitude, $distance, $sublong, $sublat) = moon_external($ctime, $retardation_time);
        return($longitude, $latitude, $distance, $sublong, $sublat);
    }
}
###############################################################################
# Calculate the Ecliptic longitude (rad), latitude (rad) and distance (m)
# of the moon for the given ctime.  Coordinates are equinox of date.
sub moon_external {

    my ($ctime, $retardation_time) = @_;
    # Get data from Horizons cache
    # First, do we need to access more data?  Yes if $ctime > the latest
    # ctime we have in the cache.  It's probably not really necessary
    # to have this here, because regeneration of the spline "looks ahead"
    # several minutes and as such the code which deals with spline
    # regeneration (below) should require more Horizons data first.
    if ($ctime >= $horizons{'moon'}->{'latest_ctime'}) {
        get_horizons_moon_data($ctime);
        clean_moon_data_cache();
        write_moon_data_cache();
    }

    # Second, do we need to regenerate the splines?
    # Check data structures are initialised first.
    unless (defined $spline{'moon'}->{'timestamps'}) {
        my @spline_moondata_timestamps;
        $spline{'moon'}->{'timestamps'} = \@spline_moondata_timestamps;
        $spline{'moon'}->{'rebuild'} = 1;
    }
    unless (defined $spline{'moon'}->{'longitudes'}) {
        my @spline_moondata_longitudes;
        $spline{'moon'}->{'longitudes'} = \@spline_moondata_longitudes;
        $spline{'moon'}->{'rebuild'} = 1;
    }
    unless (defined $spline{'moon'}->{'latitudes'}) {
        my @spline_moondata_latitudes;
        $spline{'moon'}->{'latitudes'} = \@spline_moondata_latitudes;
        $spline{'moon'}->{'rebuild'} = 1;
    }
    unless (defined $spline{'moon'}->{'distances'}) {
        my @spline_moondata_distances;
        $spline{'moon'}->{'distances'} = \@spline_moondata_distances;
        $spline{'moon'}->{'rebuild'} = 1;
    }
    unless (defined $spline{'moon'}->{'sublong'}) {
        my @spline_moondata_sublong;
        $spline{'moon'}->{'sublong'} = \@spline_moondata_sublong;
        $spline{'moon'}->{'rebuild'} = 1;
    }
    unless (defined $spline{'moon'}->{'sublat'}) {
        my @spline_moondata_sublat;
        $spline{'moon'}->{'sublat'} = \@spline_moondata_sublat;
        $spline{'moon'}->{'rebuild'} = 1;
    }

    # Regenerate the splines here.
    if ($spline{'moon'}->{'rebuild'}) {
        # Search for appropriate index in the larger data arrays.
        my $new_index = -1;
        for (my $i = 0; $i <= $#{$horizons{'moon'}->{'timestamps'}}; $i++) {
            if (${$horizons{'moon'}->{'timestamps'}}[$i] > $ctime) {
                $new_index = $i;
                last;
            }
        }
        # What happened?
        # Do we have data for the "end" of the spline?
        my $need_more_data = 0;
        $need_more_data++ if ($new_index == -1);
        $need_more_data++ unless (defined ${$horizons{'moon'}->{'timestamps'}}[$new_index - 2]);
        $need_more_data++ unless (defined ${$horizons{'moon'}->{'timestamps'}}[$new_index - 1]);
        $need_more_data++ unless (defined ${$horizons{'moon'}->{'timestamps'}}[$new_index]);
        $need_more_data++ unless (defined ${$horizons{'moon'}->{'timestamps'}}[$new_index + 1]);
        $need_more_data++ unless (defined ${$horizons{'moon'}->{'timestamps'}}[$new_index + 2]);
        unless ($need_more_data) {
            $need_more_data++ unless ((${$horizons{'moon'}->{'timestamps'}}[$new_index - 1] - ${$horizons{'moon'}->{'timestamps'}}[$new_index - 2]) == 300);
            $need_more_data++ unless ((${$horizons{'moon'}->{'timestamps'}}[$new_index] - ${$horizons{'moon'}->{'timestamps'}}[$new_index - 1]) == 300);
            $need_more_data++ unless ((${$horizons{'moon'}->{'timestamps'}}[$new_index + 1] - ${$horizons{'moon'}->{'timestamps'}}[$new_index]) == 300);
            $need_more_data++ unless ((${$horizons{'moon'}->{'timestamps'}}[$new_index + 2] - ${$horizons{'moon'}->{'timestamps'}}[$new_index + 1]) == 300);
        }
        if ($need_more_data) {
            get_horizons_moon_data($ctime);
            clean_moon_data_cache();
            write_moon_data_cache();
            # Search for appropriate index... again!
            $new_index = -1;
            for (my $i = 0; $i <= $#{$horizons{'moon'}->{'timestamps'}}; $i++) {
                if (${$horizons{'moon'}->{'timestamps'}}[$i] > $ctime) {
                    $new_index = $i;
                    last;
                }
            }
        }
        if ($new_index < 0) {
            # Oh, my!  Still no valid data at the end?  Abort in confusion.
            print "Aborting in confusion.  Cound not obtain appropriate lunar data from Horizons for ctime \"$ctime\".\n";
            exit(1);
        }
        my $start_index = $new_index - 2;
        my $end_index   = $new_index + 1;
        $start_index = 0 if ($start_index < 0);

        # Pull appropriate data into $spline
        $spline{'moon'}->{'timestamps'} = ();
        $spline{'moon'}->{'longitudes'} = ();
        $spline{'moon'}->{'latitudes'} = ();
        $spline{'moon'}->{'distances'} = ();
        $spline{'moon'}->{'sublong'} = ();
        $spline{'moon'}->{'sublat'} = ();
        foreach my $i (0 .. 3) {
            push(@{$spline{'moon'}->{'timestamps'}}, ${$horizons{'moon'}->{'timestamps'}}[$start_index + $i]);
            push(@{$spline{'moon'}->{'longitudes'}}, ${$horizons{'moon'}->{'longitudes'}}[$start_index + $i]);
            push(@{$spline{'moon'}->{'latitudes'}},  ${$horizons{'moon'}->{'latitudes'}}[$start_index + $i]);
            push(@{$spline{'moon'}->{'distances'}},  ${$horizons{'moon'}->{'distances'}}[$start_index + $i]);
            push(@{$spline{'moon'}->{'sublong'}},  ${$horizons{'moon'}->{'sublong'}}[$start_index + $i]);
            push(@{$spline{'moon'}->{'sublat'}},  ${$horizons{'moon'}->{'sublat'}}[$start_index + $i]);
        }
        # Catch splines in which 360 -> 0 wrapping in moon ecliptic longitude values takes place.  Tends to upset the derivatives.
        foreach my $i (0 .. 2) {
            if (${$spline{'moon'}->{'longitudes'}}[$i + 1] < ${$spline{'moon'}->{'longitudes'}}[$i]) {
                ${$spline{'moon'}->{'longitudes'}}[$i + 1] += 360;
            }
        }
        # Make second-order derivatives
        my @moon_longitudes2 = Derivative2($spline{'moon'}->{'timestamps'}, $horizons{'moon'}->{'longitudes'});
        $spline{'moon'}->{'longitudes2'} = \@moon_longitudes2;
        my @moon_latitudes2  = Derivative2($spline{'moon'}->{'timestamps'}, $horizons{'moon'}->{'latitudes'});
        $spline{'moon'}->{'latitudes2'} = \@moon_latitudes2;
        my @moon_distances2  = Derivative2($spline{'moon'}->{'timestamps'}, $horizons{'moon'}->{'distances'});
        $spline{'moon'}->{'distances2'} = \@moon_distances2;
        my @moon_sublong2  = Derivative2($spline{'moon'}->{'timestamps'}, $horizons{'moon'}->{'sublong'});
        $spline{'moon'}->{'sublong2'} = \@moon_sublong2;
        my @moon_sublat2  = Derivative2($spline{'moon'}->{'timestamps'}, $horizons{'moon'}->{'sublat'});
        $spline{'moon'}->{'sublat2'} = \@moon_sublat2;

        $spline{'moon'}->{'rebuild'} = 0;
    }

    # Third, interpolate splines for data.
    my $index = binsearch($spline{'moon'}->{'timestamps'}, $ctime);
    my $dist = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'distances'}, $spline{'moon'}->{'distances2'}, $index, $ctime);
    $ctime += $dist / $c;
    $ctime -= $retardation_time;
    my $long = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'longitudes'}, $spline{'moon'}->{'longitudes2'}, $index, $ctime);
    my $lat = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'latitudes'}, $spline{'moon'}->{'latitudes2'}, $index, $ctime);
    $dist = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'distances'}, $spline{'moon'}->{'distances2'}, $index, $ctime);
    my $sublong = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'sublong'}, $spline{'moon'}->{'sublong2'}, $index, $ctime);
    my $sublat = spline($spline{'moon'}->{'timestamps'}, $spline{'moon'}->{'sublat'}, $spline{'moon'}->{'sublat2'}, $index, $ctime);
    return(deg2rad($long), deg2rad($lat), $dist, $sublong, $sublat);
}
###############################################################################
# Calculate the Ecliptic longitude (rad), latitude (rad) and distance (m)
# of the moon for the given ctime.  Coordinates are equinox of date.
sub moon_internal {

    my ($ctime, $retardation_time) = @_;
    $ctime -= $retardation_time;
    my @bits = gmtime($ctime);
    my $year = $bits[5] + 1900;
    my $month = $bits[4] + 1;
    my $day = $bits[3] + $bits[2] / 24 + $bits[1] / 1440 + $bits[0] / 86400;
    if ($month < 3) {
         $month += 12;
         $year--;
    }
    my $a = int($year / 100.0);
    my $b = 2.0 - $a + int($a / 4.0);
    my $c = int(365.25 * $year)  - 694025.0;
    my $d = int(30.6001 * ($month + 1));
    my $djd = $b + $c + $d + $day - 0.5;
    my $t = $djd / 36525.0;
    my $t2 = $t * $t;

    my $m1 = 2.732158213E+01;
    my $m2 = 3.652596407E+02;
    my $m3 = 2.755455094E+01;
    my $m4 = 2.953058868E+01;
    my $m5 = 2.721222039E+01;
    my $m6 = 6.798363307E+03;
    my $q = $djd;
    $m1 = $q / $m1;
    $m2 = $q / $m2;
    $m3 = $q / $m3;
    $m4 = $q / $m4;
    $m5 = $q / $m5;
    $m6 = $q / $m6;
    $m1 = 360 * ($m1 - int($m1));
    $m2 = 360 * ($m2 - int($m2));
    $m3 = 360 * ($m3 - int($m3));
    $m4 = 360 * ($m4 - int($m4));
    $m5 = 360 * ($m5 - int($m5));
    $m6 = 360 * ($m6 - int($m6));

    my $ml = 2.70434164E+02 + $m1 - (1.133E-03 - 1.9E-06 * $t) * $t2;
    my $ms = 3.58475833E+02 + $m2 - (1.5E-04 + 3.3E-06 * $t) * $t2;
    my $md = 2.96104608E+02 + $m3 + (9.192E-03 + 1.44E-05 * $t) * $t2;
    my $me = 3.50737486E+02 + $m4 - (1.436E-03 - 1.9E-06 * $t) * $t2;
    my $mf = 11.250889 + $m5 - (3.211E-03 + 3.0E-07 * $t) * $t2;
    my $na = 2.59183275E+02 - $m6 + (2.078E-03 + 2.2E-06 * $t) * $t2;
    $a = deg2rad(51.2 + 20.2 * $t);
    my $s1 = sin($a);
    my $s2 = sin(deg2rad($na));
    $b = 346.56 + (132.87 - 9.1731E-03 * $t) * $t;
    my $s3 = 3.964E-03 * sin(deg2rad($b));
    $c = deg2rad($na + 275.05 - 2.3 * $t);
    my $s4 = sin($c);
    $ml = $ml + 2.33E-04 * $s1 + $s3 + 1.964E-03 * $s2;
    $ms = $ms - 1.778E-03 * $s1;
    $md = $md + 8.17E-04 * $s1 + $s3 + 2.541E-03 * $s2;
    $mf = $mf + $s3 - 2.4691E-02 * $s2 - 4.328E-03 * $s4;
    $me = $me + 2.011E-03 * $s1 + $s3 + 1.964E-03 * $s2;
    my $e = 1.0 - (2.495E-03 + 7.52E-6 * $t) * $t;
    my $e2 = $e * $e;
    $ml = deg2rad($ml);
    $ms = deg2rad($ms);
    $na = deg2rad($na);
    $me = deg2rad($me);
    $mf = deg2rad($mf);
    $md = deg2rad($md);

    my $l = 6.28875 * sin($md) + 1.274018 * sin(2 * $me - $md);
    $l = $l + 6.58309E-01 * sin(2 * $me) + 2.13616E-01 * sin(2 * $md);
    $l = $l - $e * 1.85596E-01 * sin($ms) - 1.14336E-01 * sin(2 * $mf);
    $l = $l + 5.8793E-02 * sin(2 * ($me - $md));
    $l = $l + 5.7212E-02 * $e * sin(2 * $me - $ms - $md) + 5.332E-02 * sin(2 * $me + $md);
    $l = $l + 4.5874E-02 * $e * sin(2 * $me - $ms) + 4.1024E-02 * $e * sin($md - $ms);
    $l = $l - 3.4718E-02 * sin($me) - $e * 3.0465E-02 * sin($ms + $md);
    $l = $l + 1.5326E-02 * sin(2 * ($me - $mf)) - 1.2528E-02 * sin(2 * $mf + $md);
    $l = $l - 1.098E-02 * sin(2 * $mf - $md) + 1.0674E-02 * sin(4 * $me - $md);
    $l = $l + 1.0034E-02 * sin(3 * $md) + 8.548E-03 * sin(4 * $me - 2 * $md);
    $l = $l - $e * 7.91E-03 * sin($ms - $md + 2 * $me) - $e * 6.783E-03 * sin(2 * $me + $ms);
    $l = $l + 5.162E-03 * sin($md - $me) + $e * 5E-03 * sin($ms + $me);
    $l = $l + 3.862E-03 * sin(4 * $me) + $e * 4.049E-03 * sin($md - $ms + 2 * $me);
    $l = $l + 3.996E-03 * sin(2 * ($md + $me)) + 3.665E-03 * sin(2 * $me - 3 * $md);
    $l = $l + $e * 2.695E-03 * sin(2 * $md - $ms) + 2.602E-03 * sin($md - 2 * ($mf + $me));
    $l = $l + $e * 2.396E-03 * sin(2 * ($me - $md) - $ms) - 2.349E-03 * sin($md + $me);
    $l = $l + $e * 2.249E-03 * sin(2 * ($me - $ms)) - $e * 2.125E-03 * sin(2 * $md + $ms);
    $l = $l - $e * 2.079E-03 * sin(2 * $ms) + $e2 * 2.059E-03 * sin(2 * ($me - $ms) - $md);
    $l = $l - 1.773E-03 * sin($md + 2 * ($me - $mf)) - 1.595E-03 * sin(2 * ($mf + $me));
    $l = $l + $e * 1.22E-03 * sin(4 * $me - $ms - $md) - 1.11E-03 * sin(2 * ($md + $mf));
    $l = $l + 8.92E-04 * sin($md - 3 * $me) - $e * 8.11E-04 * sin($ms + $md + 2 * $me);
    $l = $l + $e * 7.61E-04 * sin(4 * $me - $ms - 2 * $md);
    $l = $l + $e2 * 7.04E-04 * sin($md - 2 * ($ms + $me));
    $l = $l + $e * 6.93E-04 * sin($ms - 2 * ($md - $me));
    $l = $l + $e * 5.98E-04 * sin(2 * ($me - $mf) - $ms);
    $l = $l + 5.5E-04 * sin($md + 4 * $me) + 5.38E-04 * sin(4 * $md);
    $l = $l + $e * 5.21E-04 * sin(4 * $me - $ms) + 4.86E-04 * sin(2 * $md - $me);
    $l = $l + $e2 * 7.17E-04 * sin($md - 2 * $ms);
    my $mm = $ml + deg2rad($l);
    while ($mm < 0) {
        $mm += 2.0 * pi;
    }
    while ($mm >= (2.0 * pi)) {
        $mm -= 2.0 * pi;
    }

    my $g = 5.128189 * sin($mf) + 2.80606E-01 * sin($md + $mf);
    $g = $g + 2.77693E-01 * sin($md - $mf) + 1.73238E-01 * sin(2 * $me - $mf);
    $g = $g + 5.5413E-02 * sin(2 * $me + $mf - $md) + 4.6272E-02 * sin(2 * $me - $mf - $md);
    $g = $g + 3.2573E-02 * sin(2 * $me + $mf) + 1.7198E-02 * sin(2 * $md + $mf);
    $g = $g + 9.267E-03 * sin(2 * $me + $md - $mf) + 8.823E-03 * sin(2 * $md - $mf);
    $g = $g + $e * 8.247E-03 * sin(2 * $me - $ms - $mf) + 4.323E-03 * sin(2 * ($me - $md) - $mf);
    $g = $g + 4.2E-03 * sin(2 * $me + $mf + $md) + $e * 3.372E-03 * sin($mf - $ms - 2 * $me);
    $g = $g + $e * 2.472E-03 * sin(2 * $me + $mf - $ms - $md);
    $g = $g + $e * 2.222E-03 * sin(2 * $me + $mf - $ms);
    $g = $g + $e * 2.072E-03 * sin(2 * $me - $mf - $ms - $md);
    $g = $g + $e * 1.877E-03 * sin($mf - $ms + $md) + 1.828E-03 * sin(4 * $me - $mf - $md);
    $g = $g - $e * 1.803E-03 * sin($mf + $ms) - 1.75E-03 * sin(3 * $mf);
    $g = $g + $e * 1.57E-03 * sin($md - $ms - $mf) - 1.487E-03 * sin($mf + $me);
    $g = $g - $e * 1.481E-03 * sin($mf + $ms + $md) + $e * 1.417E-03 * sin($mf - $ms - $md);
    $g = $g + $e * 1.35E-03 * sin($mf - $ms) + 1.33E-03 * sin($mf - $me);
    $g = $g + 1.106E-03 * sin($mf + 3 * $md) + 1.02E-03 * sin(4 * $me - $mf);
    $g = $g + 8.33E-04 * sin($mf + 4 * $me - $md) + 7.81E-04 * sin($md - 3 * $mf);
    $g = $g + 6.7E-04 * sin($mf + 4 * $me - 2 * $md) + 6.06E-04 * sin(2 * $me - 3 * $mf);
    $g = $g + 5.97E-04 * sin(2 * ($me + $md) - $mf);
    $g = $g + $e * 4.92E-04 * sin(2 * $me + $md - $ms - $mf) + 4.5E-04 * sin(2 * ($md - $me) - $mf);
    $g = $g + 4.39E-04 * sin(3 * $md - $mf) + 4.23E-04 * sin($mf + 2 * ($me + $md));
    $g = $g + 4.22E-04 * sin(2 * $me - $mf - 3 * $md) - $e * 3.67E-04 * sin($ms + $mf + 2 * $me - $md);
    $g = $g - $e * 3.53E-04 * sin($ms + $mf + 2 * $me) + 3.31E-04 * sin($mf + 4 * $me);
    $g = $g + $e * 3.17E-04 * sin(2 * $me + $mf - $ms + $md);
    $g = $g + $e2 * 3.06E-04 * sin(2 * ($me - $ms) - $mf) - 2.83E-04 * sin($md + 3 * $mf);
    my $w1 = 4.664E-04 * cos($na);
    my $w2 = 7.54E-05 * cos($c);
    my $bm = deg2rad($g) * (1.0 - $w1 - $w2);

    my $pm = 9.50724E-01 + 5.1818E-02 * cos($md) + 9.531E-03 * cos(2 * $me - $md);
    $pm = $pm + 7.843E-03 * cos(2 * $me) + 2.824E-03 * cos(2 * $md);
    $pm = $pm + 8.57E-04 * cos(2 * $me + $md) + $e * 5.33E-04 * cos(2 * $me - $ms);
    $pm = $pm + $e * 4.01E-04 * cos(2 * $me - $md - $ms);
    $pm = $pm + $e * 3.2E-04 * cos($md - $ms) - 2.71E-04 * cos($me);
    $pm = $pm - $e * 2.64E-04 * cos($ms + $md) - 1.98E-04 * cos(2 * $mf - $md);
    $pm = $pm + 1.73E-04 * cos(3 * $md) + 1.67E-04 * cos(4 * $me - $md);
    $pm = $pm - $e * 1.11E-04 * cos($ms) + 1.03E-04 * cos(4 * $me - 2 * $md);
    $pm = $pm - 8.4E-05 * cos(2 * $md - 2 * $me) - $e * 8.3E-05 * cos(2 * $me + $ms);
    $pm = $pm + 7.9E-05 * cos(2 * $me + 2 * $md) + 7.2E-05 * cos(4 * $me);
    $pm = $pm + $e * 6.4E-05 * cos(2 * $me - $ms + $md) - $e * 6.3E-05 * cos(2 * $me + $ms - $md);
    $pm = $pm + $e * 4.1E-05 * cos($ms + $me) + $e * 3.5E-05 * cos(2 * $md - $ms);
    $pm = $pm - 3.3E-05 * cos(3 * $md - 2 * $me) - 3.0E-05 * cos($md + $me);
    $pm = $pm - 2.9E-05 * cos(2 * ($mf - $me)) - $e * 2.9E-05 * cos(2 * $md + $ms);
    $pm = $pm + $e2 * 2.6E-05 * cos(2 * ($me - $ms)) - 2.3E-05 * cos(2 * ($mf - $me) + $md);
    $pm = $pm + $e * 1.9E-05 * cos(4 * $me - $ms - $md);
    $pm = deg2rad($pm);
    my $dist = 6378140 / sin($pm);

    return($mm, $bm, $dist);
}
###############################################################################
sub anomaly {

    my ( $mean_anomaly, $eccentricity ) = @_;
    my $tpi = 2 * pi;
    my $m = $mean_anomaly - ($tpi * int($mean_anomaly / $tpi));
    my $eccentric_anomaly = $m;
    my $dla = $eccentric_anomaly - ($eccentricity * sin($eccentric_anomaly)) - $m;
    my $loops = 0;

    while ((abs($dla) > 1.0E-08) && ($loops < 100)) {
        $dla /= (1 - ($eccentricity * cos($eccentric_anomaly)));
        $eccentric_anomaly -= $dla;
        $dla = $eccentric_anomaly - ($eccentricity * sin($eccentric_anomaly)) - $m;
        $loops++;
    }

    my $value = sqrt((1 + $eccentricity) / (1 - $eccentricity)) * tan($eccentric_anomaly / 2);
    my $true_anomaly = 2 * atan($value);
    return($true_anomaly, $eccentric_anomaly);
}
###############################################################################
# Calculate the Ecliptic longitude (rad) and distance (m)
# of the sun for the given ctime.  Coordinates are equinox of date.
sub sun {

    my ($ctime, $retardation_time) = @_;
    my @bits = gmtime($ctime - $retardation_time);
    my $year = $bits[5] + 1900;
    my $month = $bits[4] + 1;
    my $day = $bits[3] + $bits[2] / 24 + $bits[1] / 1440 + $bits[0] / 86400;
    if ($month < 3) {
         $month += 12;
         $year--;
    }
    my $a = int($year / 100.0);
    my $b = 2.0 - $a + int($a / 4.0);
    my $c = int(365.25 * $year)  - 694025.0;
    my $d = int(30.6001 * ($month + 1));
    my $djd = $b + $c + $d + $day - 0.5;
    my $t = $djd / 36525.0;
    my $t2 = $t * $t;

    $a = 1.000021359E+02 * $t;
    $b = 360.0 * ($a - int($a));
    my $l = 2.7969668E+02 + 3.025E-04 * $t2 + $b;
    $a = 9.999736042E+01 * $t;
    $b = 360.0 * ($a - int($a));
    my $m1 = 3.5847583E+02 - (1.5E-04 + 3.3E-06 * $t) * $t2 + $b;
    my $ec = 1.675104E-02 - 4.18E-05 * $t - 1.26E-07 * $t2;

    my $am = deg2rad($m1);
    my ($at, $ae) = &anomaly($am, $ec);

    $a = 6.255209472E+01 * $t;
    $b = 360.0 * ($a - int($a));
    my $a1 = deg2rad(153.23 + $b);
    $a = 1.251041894E+02 * $t;
    $b = 360.0 * ($a - int($a));
    my $b1 = deg2rad(216.57 + $b);
    $a = 9.156766028E+01 * $t;
    $b = 360.0 * ($a - int($a));
    my $c1 = deg2rad(312.69 + $b);
    $a = 1.236853095E+03 * $t;
    $b = 360.0 * ($a - int($a));
    my $d1 = deg2rad(350.74 - 1.44E-03 * $t2 + $b);
    my $e1 = deg2rad(231.19 + 20.2 * $t);
    $a = 1.831353208E+02 * $t;
    $b = 360.0 * ($a - int($a));
    my $h1 = deg2rad(353.4 + $b);

    my $d2 = 1.34E-03 * cos($a1) + 1.54E-03 * cos($b1) + 2.0E-03 * cos($c1) + 1.79E-03 * sin($d1) + 1.78E-03 * sin($e1);
    my $d3 = 5.43E-06 * sin($a1) + 1.575E-05 * sin($b1) + 1.627E-05 * sin($c1) + 3.076E-05 * cos($d1) + 9.27E-06 * sin($h1);

    my $sr = $at + deg2rad ($l - $m1 + $d2);
    my $rr = 1.49598E+11 * (1.0000002 * (1.0 - $ec * cos($ae)) + $d3);
    $sr += 2 * pi if ($sr < 0);
    $sr -= 2 * pi if ($sr > 2 * pi);

    return($sr, $rr);
}
###############################################################################
# Lunar libration stuff: reserved for future expansion
sub libration {

    ###
    ### Returns the Lunar libration in longitude and latitude.
    ###
    my ($ctime, $long, $lat) = @_;
    my @bits = gmtime($ctime);
    my $year = $bits[5] + 1900;
    my $month = $bits[4] + 1;
    my $day = $bits[3] + $bits[2] / 24 + $bits[1] / 1440 + $bits[0] / 86400;
    if ($month < 3) {
         $month += 12;
         $year--;
    }
    my $a = int($year / 100.0);
    my $b = 2.0 - $a + int($a / 4.0);
    my $c = int(365.25 * $year)  - 694025.0;
    my $d = int(30.6001 * ($month + 1));
    my $djd = $b + $c + $d + $day - 36525.5;
    my $t = $djd / 36525.0;

    my $f = deg2rad(93.2720993 + 483202.0175273 * $t - .0034029 * $t * $t);
    $b = $f / (2 * pi);
    $a = 2 * pi * ($b - floor($b));
    $a += (2 * pi) if ($a < 0);
    $f = $a;

    my $omega = deg2rad(125.044555 - 1934.1361849 * $t + .0020762 * $t * $t);
    $b = $omega / (2 * pi);
    $a = 2 * pi * ($b - floor($b));
    $a += (2 * pi) if ($a < 0);
    $omega = $a;

    #my $w = $long - $omega;
    my $w = $omega - $long;
    my $y = -1 * sin($w) * cos($lat) * cos($moon_axial_tilt) - sin($lat) * sin($moon_axial_tilt);
    my $x = cos($w) * cos($lat);
    $a = atan2($y, $x);
    my $ell = $a - $f;

    # kludge to catch cases of 'round the back' angles
    $ell += (2 * pi) if ($ell < (pi / -2));
    $ell -= (2 * pi) if ($ell > (pi / 2));

    my $bee = asin(sin($w) * cos($lat) * sin($moon_axial_tilt) - sin($lat) * cos($moon_axial_tilt));

    return(rad2deg($ell), rad2deg($bee));
}
###############################################################################
# Convert cartesian coordinates (positions or velocities) between
# planes.  This is currently used to convert the state vector of
# Numbat from a coordinate system aligned with the ecliptic to one
# aligned equatorially.  This is so we can convert the state vector to
# orbital elements relative to both planes, as the elements relative
# to both planes are of interest.  Inputs are (x, y, z) in "ecliptic
# space" and outputs are (x, y, z) in "equatorial" space.  This is
# implemented using a straight-out transformation matrix, of the form
#
#                 | 1         0          0   |
#  R_equatorial = | 0       cos T     -sin T | * R_ecliptic
#                 | 0       sin T      cos T |
#
# where T is the obliquity of the ecliptic for epoch J2000.0.
sub planechange {

    my ($x, $y, $z, $T) = @_;

    my $new_x = $x;
    my $new_y = cos($T) * $y - sin($T) * $z;
    my $new_z = sin($T) * $y + cos($T) * $z;
    return ($new_x, $new_y, $new_z);
}
###############################################################################
# Calculate cartesian components of gravitational acceleration
# by the earth at the given coordinates.
sub earth_gravity {

    my ($x_position, $y_position, $z_position) = @_;

    # First, shift the origin of the cartesian coordinates to Numbat's position.
    # This allows us to get cartesian coordinates of Earth as seen by Numbat.
    my $delta_x = 0 - $x_position;
    my $delta_y = 0 - $y_position;
    my $delta_z = 0 - $z_position;
    # Transform these cartesian coordinates to Ecliptic.
    # This gives Ecliptic coordinates of the Earth as seen from Numbat.
    my $lambda = atan2($delta_y, $delta_x);
    my $beta   = atan($delta_z / sqrt($delta_x * $delta_x + $delta_y * $delta_y));
    my $earth_distance2 = $delta_x * $delta_x + $delta_y * $delta_y + $delta_z * $delta_z;
    # Calculate magnitude of the acceleration vector
    my $j2 = 1.0827E-03;
    my $acceleration = $G * $Me / $earth_distance2;

    # Quadrupole factor
    my ($x, $y, $z) = planechange($x_position, $y_position, $z_position, $earth_axial_tilt);
    my $sdec = sin(atan($z / sqrt($x * $x + $y * $y)));
    my $factor = 1 - (3 * $j2 * $Re * $Re * (1.5 * ($sdec * $sdec) - 0.5) / $earth_distance2);
    # Use magnitude and ecliptic coordinates (direction) to resolve this acceleration
    # in to cartesian components.
    my $z_component = $acceleration * $factor * sin($beta);
    my $y_component = $acceleration * $factor * cos($beta) * sin($lambda);
    my $x_component = $acceleration * $factor * cos($beta) * cos($lambda);
    return($x_component, $y_component, $z_component);

}
###############################################################################
# Calculate cartesian components of gravitational acceleration
# by the moon at the given coordinates and ctime.  Also returns
# the light travel time to the moon so we can retard by this amount
# in the next loop.  Using a retarded potential is probably overkill,
# as the influence of the Moon is greatest when the craft is nearest,
# at which point the retardation time will be small.  Anyway...
sub moon_gravity {

    my ($x_position, $y_position, $z_position, $time, $retardation_time) = @_;

    my ($longitude, $latitude, $distance, undef, undef) = moon($time, $retardation_time);
    # Correct for luni-solar precession.  This is a small < 0 adjustment to Ecliptic longitude.
    # Ignore retardation time here
    $longitude += ($J2k_epoch - $time) * $luni_solar_precession_rate;
    # Convert the moon's Ecliptic coordinates into cartesian coordinates
    my $z = $distance * sin($latitude);
    my $y = $distance * cos($latitude) * sin($longitude);
    my $x = $distance * cos($latitude) * cos($longitude);
    # Now shift to get the displacement of the moon from Numbat.
    # This allows us to get cartesian coordinates of the Moon as seen by Numbat.
    my $delta_x = $x - $x_position;
    my $delta_y = $y - $y_position;
    my $delta_z = $z - $z_position;
    # Transform these cartesian coordinates to Ecliptic.
    # This gives Ecliptic coordinates of the Moon as seen from Numbat.
    my $lambda = atan2($delta_y, $delta_x);
    my $beta   = atan($delta_z / sqrt($delta_x * $delta_x + $delta_y * $delta_y));
    my $moon_distance2 = $delta_x * $delta_x + $delta_y * $delta_y + $delta_z * $delta_z;
    my $acceleration = $G * $Mm / $moon_distance2;
    # Use magnitude and ecliptic coordinates (direction) to resolve this acceleration
    # in to cartesian components.
    my $z_component = $acceleration * sin($beta);
    my $y_component = $acceleration * cos($beta) * sin($lambda);
    my $x_component = $acceleration * cos($beta) * cos($lambda);
    my $new_retardation_time = sqrt($moon_distance2) / $c;
    return($x_component, $y_component, $z_component, $new_retardation_time);
}
###############################################################################
# Calculate cartesian components of gravitational acceleration
# by the sun at the given coordinates and ctime.  Also returns
# the light travel time to the sun so we can retard by this amount
# in the next loop.
sub sun_gravity {

    my ($x_position, $y_position, $z_position, $time, $retardation_time) = @_;

    my ($longitude, $distance) = sun($time, $retardation_time);
    # Correct for luni-solar precession.  This is a small < 0 adjustment to Ecliptic longitude.
    # Ignore retardation time here
    $longitude += ($J2k_epoch - $time) * $luni_solar_precession_rate;
    # Also correct for aberration, pretty much a constant 20.5 arcsec for the sun.
    $longitude -= 9.93E-05;
    # Convert the sun's Ecliptic coordinates into cartesian coordinates
    my $z = 0;
    my $y = $distance * sin($longitude);
    my $x = $distance * cos($longitude);
    # Now shift to get the displacement of the sun from Numbat.
    # This allows us to get cartesian coordinates of the sun as seen by Numbat.
    my $delta_x = $x - $x_position;
    my $delta_y = $y - $y_position;
    my $delta_z = $z - $z_position;
    # Transform these cartesian coordinates to Ecliptic.
    # This gives Ecliptic coordinates of the Sun as seen from Numbat.
    my $lambda = atan2($delta_y, $delta_x);
    my $beta   = atan($delta_z / sqrt($delta_x * $delta_x + $delta_y * $delta_y));
    my $sun_distance2 = $delta_x * $delta_x + $delta_y * $delta_y + $delta_z * $delta_z;
    my $acceleration = $G * $Ms / $sun_distance2;
    # Use magnitude and ecliptic coordinates (direction) to resolve this acceleration
    # in to cartesian components.
    my $z_component = $acceleration * sin($beta);
    my $y_component = $acceleration * cos($beta) * sin($lambda);
    my $x_component = $acceleration * cos($beta) * cos($lambda);
    my $new_retardation_time = sqrt($sun_distance2) / $c;
    return($x_component, $y_component, $z_component, $new_retardation_time);
}
###############################################################################
# Calculate cartesian components of acceleration due to on-board engine.
# The ctime is used to determine engine states (from the @engine_states array).
# The acceleration is oriented in the direction of the Numbat's instantaneous
# velocity with respect to Earth, or in the direction exactly opposite to this.
# This will need to be made more flexible if this model is to be used to
# adjust the lunar orbit.
sub engine {

    my ($time) = @_;

    # Determine engine state for this time $t
    my $new_engine_state = 0;
    my $ran_out_of_fuel = 0;
    foreach my $entry (@engine_states) {
        my ($ctime, $state) = split(/\|/, $entry);
        $new_engine_state = $state;
        last if ($ctime < $time);
    }

    # Your mower is out of petrol, and you were mowing at the time
    if (($$numbat{'m'} <= ($initial_mass - $fuel_mass)) && ($$numbat{'engine'} != 0)) {
        $new_engine_state = 0;
        $ran_out_of_fuel = 1;
    }
    # Your mower is out of petrol, but you weren't mowing anyway
    $new_engine_state = 0 if ($$numbat{'m'} <= ($initial_mass - $fuel_mass));

    my $x_component = my $y_component = my $z_component = 0;
    if ($new_engine_state) {
        # Speed of the craft, relative to the Earth.
        my $dx = $$numbat{'dx'};
        my $dy = $$numbat{'dy'};
        my $dz = $$numbat{'dz'};
        my $speed = sqrt($dx * $dx + $dy * $dy + $dz * $dz);
        my $mass = $$numbat{'m'};

        $x_component = ($engine_force / $mass) * ($dx / $speed) * $new_engine_state;
        $y_component = ($engine_force / $mass) * ($dy / $speed) * $new_engine_state;
        $z_component = ($engine_force / $mass) * ($dz / $speed) * $new_engine_state;
    }
    return($new_engine_state, $ran_out_of_fuel, $x_component, $y_component, $z_component);
}
###############################################################################
# Calculate instantaneous, osculating orbital elements from Numbat state vector
# relative to another state vector.
# Returns: Semi-major axis (m),
#          Eccentriciy (dimensionless),
#          Inclination (degrees),
#          Longitude of the ascending node (degrees),
#          Longitude of periapsis (degrees),
#          Argument of periapsis (degrees),
#          Mean anomaly (degrees).
#          True anomaly (degrees).
#          Periapsis distance (m)
#          Apoapsis distance (m)
#          Orbital period (s)
#
sub state2elements {

    my ($x, $y, $z, $dx, $dy, $dz, $mass, $type) = @_;

    # Get position components relative to r
    $x *= -1;
    $y *= -1;
    $z *= -1;

    # Get velocity components relative to dr
    $dx *= -1;
    $dy *= -1;
    $dz *= -1;

    # Rotate coordinate system if required
    if ($type eq 'equatorial') {
        ($x, $y, $z) = planechange($x, $y, $z, $earth_axial_tilt);
        ($dx, $dy, $dz) = planechange($dx, $dy, $dz, $earth_axial_tilt);
    }

    # Dot products
    my $r = sqrt($x * $x + $y * $y + $z * $z);
    my $dr = sqrt($dx * $dx + $dy * $dy + $dz * $dz);

    # Form h (r x dr)
    my $hx = $y * $dz - $z * $dy;
    my $hy = $z * $dx - $x * $dz;
    my $hz = $x * $dy - $y * $dx;

    # Inclination
    my $i = atan2(sqrt($hx * $hx + $hy * $hy) , $hz);

    # Longitude of the Ascending Node
    my $lan = atan2($hx, (-1 * $hy));
    $lan += 2 * pi if ($lan < 0);

    # Argument of latitude
    my $h = sqrt($hx * $hx + $hy * $hy + $hz * $hz);
    my $u = atan2($z * $h, ($y * $hx - $x * $hy)) + pi;
    $u -= 2 * pi if ($u >= 2 * pi);

    # Length of the semi-major axis
    my $a = 1 / ((2 / $r) - ($dr * $dr / ($G * $mass)));

    # Eccentricity
    my $ecE = 1 - $r / $a;
    my $esE = ($x * $dx + $y * $dy + $z * $dz) / sqrt($G * $mass * $a);
    my $e = sqrt($ecE * $ecE + $esE * $esE);

    # Eccentric anomaly
    my $E = atan2($esE, $ecE);

    # Mean anomaly
    my $M = $E - $esE;

    # True anomaly
    my $T = atan2((sqrt(1 - $e * $e) * $esE), ($ecE - $e * $e));

    # Argument of periapsis
    my $ap = $u - $T;
    $ap += 2 * pi if ($ap < 0);

    # Longitude of periapsis
    my $lp = $ap + $lan;
    $lp -= 2 * pi if ($lp >= 2 * pi);

    $T += 2 * pi if ($T < 0);
    $M += 2 * pi if ($M < 0);

    my $periapsis_distance = $a * (1 - $e);
    my $apoapsis_distance = $a * (1 + $e);
    my $orbital_period = 2 * pi * sqrt($a * $a * $a / ($G * $mass));

    return($a, $e, rad2deg($i), rad2deg($lan), rad2deg($lp), rad2deg($ap), rad2deg($M), rad2deg($T), $periapsis_distance, $apoapsis_distance, $orbital_period);


}
###############################################################################







###############################################################################
#############################  Main program  ##################################
###############################################################################

# No buffering on STDOUT
$| = 1;



# Start with a header
print_header();

# Define some global variables
my $engine_delta_v = my $engine_run_time = my $earth_distance = my $moon_distance = 0;
my $header_print_counter = 0;
my $last_report = 0;
my $moon_retardation_time = 1;
my $sun_retardation_time = 500;
my $t = $initial_ctime;
my $moon_node_longitude = -999;
my $max_moon_node_longitude_rate = 2.0E-04;

# Make sure default engine state is 'off'
push(@engine_states, '0|0');

# Initialise Horizons array if required
if ($use_jpl_horizons_data) {
    my $rv = read_moon_data_cache();
    if ($rv) {
        # Revert to internal sources of position data on error
        $use_jpl_horizons_data = 0;
        print "No Horizons - sorry.\n";
    }
    # We've changed the content of the cache arrays, so let's sanity-check just to be safe.
    clean_moon_data_cache();
    # If we're using Horizons data, we need to fit a spline to it.
    $spline{'moon'}->{'rebuild'} = 1 unless ($rv);
}


# Start loopin'
while (1) {

    unless (defined $$numbat{'d2x'}) {
        # No data for iteration "i + 1" passed to us via the hash from the last loop.
        # Oh, well, better calculate it here then.

        ###############################################
        # Calculate acceleration due to Earth gravity #
        ###############################################
        my ($ax, $ay, $az) = earth_gravity($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'});
        $$numbat{'d2x'} = $ax;
        $$numbat{'d2y'} = $ay;
        $$numbat{'d2z'} = $az;

        ###############################################
        # Calculate acceleration due to lunar gravity #
        ###############################################
        ($ax, $ay, $az, $moon_retardation_time) = moon_gravity($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, $t, $moon_retardation_time);
        $$numbat{'d2x'} += $ax;
        $$numbat{'d2y'} += $ay;
        $$numbat{'d2z'} += $az;

        ###############################################
        #### Apply any force from a Numbat engine #####
        ###############################################
        (undef, undef, $ax, $ay, $az) = engine($t);
        $$numbat{'d2x'} += $ax;
        $$numbat{'d2y'} += $ay;
        $$numbat{'d2z'} += $az;
    }



    ###
    ### Handle engine turn on or turn off
    ###
    my ($new_engine_state, $out_of_fuel, undef, undef, undef) = engine($t);

    # Your mower is out of petrol, and you were mowing at the time
    if ($out_of_fuel) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        printf("%20s     Note: you are out of fuel!  Sorry this sounds so much like a video game.\n", $date);
    }

    # Handle turning on the engine
    if (($new_engine_state != 0) && ($$numbat{'engine'} == 0)) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        $engine_delta_v = 0;
        $engine_run_time = 0;
        printf("%20s     Note: engine state has changed to ON.\n", $date);
    }
    # Handle turning off the engine
    if (($new_engine_state == 0) && ($$numbat{'engine'} != 0)) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        printf("%20s     Note: engine state has changed to OFF.\n", $date);
        printf("%20s     Note: burn duration was %4.2f seconds.\n", $date, $engine_run_time);
        printf("%20s     Note: change in speed from this burn was %4.2f m/s.\n", $date, $engine_delta_v);
    }
    $$numbat{'engine'} = $new_engine_state;
    if ($new_engine_state) {
        $engine_delta_v += $dt * $new_engine_state * $engine_force / $$numbat{'m'};
        $engine_run_time += $dt;
    }





    ###############################################
    ##############   REPORT LINE  #################
    ###############################################
    if (($t - $last_report) >= $report_every) {
        # Make simulation "run time"
        my $duration = $t - $initial_ctime;
        my $simulated_time = '';
        if ($duration < 3600) {
            $simulated_time = sprintf("%3.1f mins", $duration / 60);
        } elsif ($duration < 86400) {
            $simulated_time = sprintf("%4.2f hrs", $duration / 3600);
        } else {
            $simulated_time = sprintf("%5.3f days", $duration / 86400);
        }

        # Make "current time" of the simulation
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));

        # Make Equatorial coordinates of Numbat (rad)
        my ($x, $y, $z) = planechange($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, $earth_axial_tilt);
        my $ra = atan2($y, $x);
        $ra += 2 * pi if ($ra < 0);
        my $dec = atan($z / sqrt($x * $x + $y * $y));

        # Earth distance
        $earth_distance = sqrt($$numbat{'x'} * $$numbat{'x'} + $$numbat{'y'} * $$numbat{'y'} + $$numbat{'z'} * $$numbat{'z'});

        # Moon distance
        my ($longitude, $latitude, $distance, $sublong, $sublat) = moon($t, $moon_retardation_time);
        $longitude += ($J2k_epoch - $t) * $luni_solar_precession_rate;
        my $moon_z = $distance * sin($latitude);
        my $moon_y = $distance * cos($latitude) * sin($longitude);
        my $moon_x = $distance * cos($latitude) * cos($longitude);
        my $delta_x = $$numbat{'x'} - $moon_x;
        my $delta_y = $$numbat{'y'} - $moon_y;
        my $delta_z = $$numbat{'z'} - $moon_z;
        $moon_distance = sqrt($delta_x * $delta_x + $delta_y * $delta_y + $delta_z * $delta_z);

        # Ecliptic longitude and latitude of Earth as seen from Moon.
        my $earth_long_as_seen_from_moon = $longitude + pi;
        my $earth_lat_as_seen_from_moon  = -1 * $latitude;

        # Ecliptic longitude and latitude of Numbat as seen from Moon
        my $numbat_long_as_seen_from_moon = atan2($delta_y, $delta_x);
        my $numbat_lat_as_seen_from_moon  = atan($delta_z / sqrt($delta_x * $delta_x + $delta_y * $delta_y));

        # Speed relative to Earth
        my $earth_speed = sqrt($$numbat{'dx'} * $$numbat{'dx'} + $$numbat{'dy'} * $$numbat{'dy'} + $$numbat{'dz'} * $$numbat{'dz'});

        # Speed of Numbat relative to moon
        # Make moon cartesian coords at time (t - 1) and (t + 1), and use this to determine approximate velocity of moon
        ($longitude, $latitude, $distance, undef, undef) = moon($t - 1, $moon_retardation_time);
        $longitude += ($J2k_epoch - $t) * $luni_solar_precession_rate;
        my $moon_z_at_t_minus_1 = $distance * sin($latitude);
        my $moon_y_at_t_minus_1 = $distance * cos($latitude) * sin($longitude);
        my $moon_x_at_t_minus_1 = $distance * cos($latitude) * cos($longitude);
        ($longitude, $latitude, $distance, undef, undef) = moon($t + 1, $moon_retardation_time);
        $longitude += ($J2k_epoch - $t) * $luni_solar_precession_rate;
        my $moon_z_at_t_plus_1 = $distance * sin($latitude);
        my $moon_y_at_t_plus_1 = $distance * cos($latitude) * sin($longitude);
        my $moon_x_at_t_plus_1 = $distance * cos($latitude) * cos($longitude);
        my $moon_dz = ($moon_z_at_t_plus_1 - $moon_z_at_t_minus_1) / 2;
        my $moon_dy = ($moon_y_at_t_plus_1 - $moon_y_at_t_minus_1) / 2;
        my $moon_dx = ($moon_x_at_t_plus_1 - $moon_x_at_t_minus_1) / 2;
        # Get velocity components relative to moon
        my $dx = $moon_dx - $$numbat{'dx'};
        my $dy = $moon_dy - $$numbat{'dy'};
        my $dz = $moon_dz - $$numbat{'dz'};
        my $moon_speed = sqrt($dx * $dx + $dy * $dy + $dz * $dz);

        # Orbital elements of Moon
        my @elements = state2elements($moon_x, $moon_y, $moon_z, $moon_dx, $moon_dy, $moon_dz, $Me, 'ecliptic');
        # Calculate cartesian components of Longitude of Ascending Node direction (z = 0 here because the ascending node lies in the Ecliptic plane).
        my $moon_lan_x = cos(deg2rad($elements[3]));
        my $moon_lan_y = sin(deg2rad($elements[3]));
        # Calculate cartesian components of Longitude of Perihelion direction
        my $moon_lp_z = sin(deg2rad($elements[2]));
        my $moon_lp_y = cos(deg2rad($elements[2])) * sin(deg2rad($elements[4]));
        my $moon_lp_x = cos(deg2rad($elements[2])) * cos(deg2rad($elements[4]));
        # Form cross product component in ecliptic plane.  This gives us the ecliptic longitude of the normal to the orbital plane.
        my $moon_orbit_plane_x = $moon_lan_y * $moon_lp_z;
        my $moon_orbit_plane_y = -1 * $moon_lan_x * $moon_lp_z;
        my $moon_orbit_plane_angle = atan2($moon_orbit_plane_y, $moon_orbit_plane_x);
        # To get the ecliptic longitude of the moon's polar axis, juat add pi radians!
        # This is because the Moon is in a "Cassini state".  Google it.
        my $moon_polar_angle = $moon_orbit_plane_angle + pi;

        # Moon polar angle correction to latlong.
        # Here we define "adjusted longitude" as the zero point of longitude.  This is analogous to the
        # Earth's longitude zero point, where the Ecliptic and equatorial planes cross.
        # We go through these hoops so we can properly account for the inclination of the moon's polar axis.
        # That is, circling the moon in the cliptic plane would see your selenographic latitude
        # change sinusoidally from ~-1.5 degrees to ~+1.5 degrees.  We obviously want to account for this
        # and hopefully this code does that to a reasonable approximation.
        my $adjusted_longitude = $earth_long_as_seen_from_moon - $moon_polar_angle + pi / 2;
        my $moon_latitude_at_ecliptic  = $sublat  + rad2deg($latitude) * cos($moon_axial_tilt * cos($adjusted_longitude));
        my $moon_longitude_at_ecliptic = $sublong - rad2deg($latitude) * sin($moon_axial_tilt * cos($adjusted_longitude));

        # Ecliptic longitude of Numbat, relative to Moon's pole
        $adjusted_longitude = $numbat_long_as_seen_from_moon - ($moon_polar_angle - pi / 2);
        $z = $moon_distance * sin($numbat_lat_as_seen_from_moon);
        $y = $moon_distance * cos($numbat_lat_as_seen_from_moon) * sin($adjusted_longitude);
        $x = $moon_distance * cos($numbat_lat_as_seen_from_moon) * cos($adjusted_longitude);
        ($x, $y, $z) = planechange($x, $y, $z, $moon_axial_tilt);
        my $new_longitude = atan2($y, $x);
        my $new_latitude  = atan($z / sqrt($x * $x + $y * $y));
        my $numbat_selenographic_latitude  = rad2deg($new_latitude);
        my $numbat_selenographic_longitude = $moon_longitude_at_ecliptic  + (rad2deg($moon_polar_angle) - 90) - rad2deg($earth_long_as_seen_from_moon) + rad2deg($new_longitude);
        while ($numbat_selenographic_longitude >= 360) {
            $numbat_selenographic_longitude -= 360;
        }
        # Altitude above the (nominal) lunar surface
        my $altitude = $moon_distance - $Rm;

        # Speed of approach to the lunar surface
        my $vv = -1 * ($delta_x * $dx + $delta_y * $dy + $delta_z * $dz) / $moon_distance;
        # Get vector components of motion tangential to surface...
        # ...in plane parallel to the ecliptic...
        my $hv_x = -1 * $dx - $vv * $delta_x / $moon_distance;
        my $hv_y = -1 * $dy - $vv * $delta_y / $moon_distance;
        my $hv_xy = sqrt($hv_x * $hv_x + $hv_y * $hv_y);
        # ...and perpendicular to the ecliptic.
        my $hv_z = -1 * $dz - $vv * $delta_z / $moon_distance;
        # What direction is xy?  Do a cross product to find out!
        # We're interested in the sign of the z-component of the product.
        my $sign = $delta_y * $dx - $delta_x * $dy;
        $hv_xy *= -1 if ($sign < 0);
        # Add a small, negative velocity to the hv_xy component here
        # to account for the moon's rotation.  Here we use an (equatorial)
        # rotation speed of 4.627 m/s.  As the moon rotates from west to east,
        # this is equivalent to a westward speed iin m/s of 4.627 * cos(latitude).
        $hv_xy -= 4.627 * cos(deg2rad($numbat_selenographic_latitude));
        my $hv =  sqrt($hv_z * $hv_z + $hv_xy * $hv_xy);
        my $angle = atan2($hv_xy, $hv_z);
        # Adjust for inclination of moon's polar axis to ecliptic.
        $angle -= $moon_axial_tilt * cos($adjusted_longitude);

        # Velocity of EMB relative to Earth
        my $emb_dx = $moon_dx * $Mm / ($Me + $Mm);
        my $emb_dy = $moon_dy * $Mm / ($Me + $Mm);
        my $emb_dz = $moon_dx * $Mm / ($Me + $Mm);

        # Speed of Numbat relative to EMB
        $dx = $emb_dx - $$numbat{'dx'};
        $dy = $emb_dy - $$numbat{'dy'};
        $dz = $emb_dz - $$numbat{'dz'};
        my $emb_speed = sqrt($dx * $dx + $dy * $dy + $dz * $dz);

        # Total energy of Numbat, ignoring gravitational potential of the Moon and with kinetic energy relative to Earth
        my $kinetic_energy = 0.5 * $$numbat{'m'} * ($$numbat{'dx'} * $$numbat{'dx'} + $$numbat{'dy'} * $$numbat{'dy'} + $$numbat{'dz'} * $$numbat{'dz'});
        my $potential_energy = -1 * $G * $$numbat{'m'} * $Me / $earth_distance;
        my $total_energy_wrt_earth = $kinetic_energy + $potential_energy;
        # Total energy of Numbat, ignoring gravitational potential of the Earth and with kinetic energy relative to Moon
        $kinetic_energy = 0.5 * $$numbat{'m'} * $moon_speed * $moon_speed;
        $potential_energy = -1 * $G * $$numbat{'m'} * $Mm / $moon_distance;
        my $total_energy_wrt_moon = $kinetic_energy + $potential_energy;

        # Print column headers if needed
        if (($header_print_counter >= $header_every) && ($header_every)) {
            print_header();
            $header_print_counter = 0;
        }
        $header_print_counter++ if ($header_every);

        # Print state information
        printf("%14.3f   %20s   %10s   %13.6e %13.6e %13.6e   %13.6e %13.6e %13.6e   %13.6e %13.6e %13.6e    %8.3f   %3s   %8.5f %+8.4f   %13.6e %5.3f   %13.6e %5.3f   %13.6e      %13.6e      %13.6e    %9.2e %9.2e   %13.6e %13.6e %13.6e   %7.3f %+7.3f %13.6e %13.6e %13.6e %+8.4f\n",
            $t, $date, $simulated_time, $$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, $$numbat{'dx'}, $$numbat{'dy'}, $$numbat{'dz'},
            $$numbat{'d2x'}, $$numbat{'d2y'}, $$numbat{'d2z'}, $$numbat{'m'}, ($$numbat{'engine'} ? 'on' : 'off'), rad2deg($ra / 15), rad2deg($dec),
            $earth_distance, $earth_distance / $Re, $moon_distance, $moon_distance / $Rm, $earth_speed, $emb_speed, $moon_speed,
            $total_energy_wrt_earth, $total_energy_wrt_moon, $moon_x, $moon_y, $moon_z,
            $numbat_selenographic_longitude, $numbat_selenographic_latitude, $altitude, $vv, $hv, rad2deg($angle));

        # Calculate orbital element data if appropriate
        if ($print_orbit_element_data) {
            if (($total_energy_wrt_earth < 0) && ($total_energy_wrt_moon > 0)) {
                # Looks like we're in Earth orbit.  Calculate orbital elements around Earth, relative to both Equatorial and Ecliptic planes.
                my @elements = state2elements($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, $$numbat{'dx'}, $$numbat{'dy'}, $$numbat{'dz'}, $Me, 'ecliptic');
                printf("%14.3f   %20s   %10s   Orbital elements, referred to Earth ecliptic: a=%13.6e e=%1.5f i=%6.2f Omega=%6.2f LP=%6.2f omega=%6.2f M=%6.2f T=%6.2f  PD=%13.6e AD=%13.6e P=%8.1f\n",
                    $t, $date, $simulated_time, @elements);
                @elements = state2elements($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, $$numbat{'dx'}, $$numbat{'dy'}, $$numbat{'dz'}, $Me, 'equatorial');
                printf("%14.3f   %20s   %10s   Orbital elements, referred to Earth equator:  a=%13.6e e=%1.5f i=%6.2f Omega=%6.2f LP=%6.2f omega=%6.2f M=%6.2f T=%6.2f  PD=%13.6e AD=%13.6e P=%8.1f\n",
                    $t, $date, $simulated_time, @elements);
            } elsif ($total_energy_wrt_moon < 0) {
                # Looks like we're in lunar orbit.  Calculate orbital elements around the Moon, relative to the Ecliptic plane.
                my @elements = state2elements(($$numbat{'x'} - $moon_x), ($$numbat{'y'} - $moon_y), ($$numbat{'z'} - $moon_z), ($$numbat{'dx'} - $moon_dx), ($$numbat{'dy'} - $moon_dy), ($$numbat{'dz'} - $moon_dz), $Mm, 'ecliptic');
                printf("%14.3f   %20s   %10s   Orbital elements, referred to Moon  ecliptic:  a=%13.6e e=%1.5f i=%6.2f Omega=%6.2f LP=%6.2f omega=%6.2f M=%6.2f T=%6.2f  PD=%13.6e AD=%13.6e P=%8.1f\n",
                    $t, $date, $simulated_time, @elements);
            }
        }
        $last_report = $t;
    }




    ###############################################
    ##############   INTEGRATOR   #################
    ###############################################
    # At this stage we have position, velocity and acceleration for time $t
    # Now we must create position at time $t+$dt so we can calculate
    # acceleration for time $t+$dt for leapfrog
    # So, first update the position vector thus
    # (you remember your old "s = u t + 0.5 a t**2" from school, don't you?)
    $$numbat{'x'} += ($$numbat{'dx'} * $dt) + (0.5 * $$numbat{'d2x'} * $dt * $dt);
    $$numbat{'y'} += ($$numbat{'dy'} * $dt) + (0.5 * $$numbat{'d2y'} * $dt * $dt);
    $$numbat{'z'} += ($$numbat{'dz'} * $dt) + (0.5 * $$numbat{'d2z'} * $dt * $dt);

    # Now use new position vector to make acceleration at time $t + $dt
    my ($mx, $my, $mz);
    my ($next_d2x, $next_d2y, $next_d2z) = earth_gravity($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'});
    ($mx, $my, $mz, $moon_retardation_time) = moon_gravity($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, ($t + $dt), $moon_retardation_time);
    $next_d2x += $mx;
    $next_d2y += $my;
    $next_d2z += $mz;

    # Differential effects of solar gravity
    my ($longitude, $latitude, $distance, undef, undef) = moon(($t + $dt), $moon_retardation_time);
    $longitude += ($J2k_epoch - $t) * $luni_solar_precession_rate;
    my $moon_z = $distance * sin($latitude);
    my $moon_y = $distance * cos($latitude) * sin($longitude);
    my $moon_x = $distance * cos($latitude) * cos($longitude);
    # Cartesian components of distance from Earth to EMB
    my $emb_x = $moon_x * $Mm / ($Me + $Mm);
    my $emb_y = $moon_y * $Mm / ($Me + $Mm);
    my $emb_z = $moon_z * $Mm / ($Me + $Mm);
    my ($nx, $ny, $nz, $ex, $ey, $ez);
    # This is how we add the differential acceleration from the sun:
    # we add the gravitational acceleration of the sun as seen by Numbat,
    # and subtract the gravitational acceleration of the sun at the Earth-Moon barycentre.
    ($nx, $ny, $nz, $sun_retardation_time) = sun_gravity($$numbat{'x'}, $$numbat{'y'}, $$numbat{'z'}, ($t + $dt), $sun_retardation_time);
    ($ex, $ey, $ez, $sun_retardation_time) = sun_gravity($emb_x, $emb_y, $emb_z, ($t + $dt), $sun_retardation_time);
    $next_d2x += $nx - $ex;
    $next_d2y += $ny - $ey;
    $next_d2z += $nz - $ez;

    # Acceleration due to engine force
    (undef, undef, $mx, $my, $mz) = engine($t + $dt);
    $next_d2x += $mx;
    $next_d2y += $my;
    $next_d2z += $mz;

    # Now use these data to update the velocity.  Think "v = u + a t" although here a is a mean of a(t) and a(t + dt).
    $$numbat{'dx'} += 0.5 * $dt * ($$numbat{'d2x'} + $next_d2x);
    $$numbat{'dy'} += 0.5 * $dt * ($$numbat{'d2y'} + $next_d2y);
    $$numbat{'dz'} += 0.5 * $dt * ($$numbat{'d2z'} + $next_d2z);

    # Modify mass if engine is on
    $$numbat{'m'} -= $mass_loss * $dt * abs($$numbat{'engine'});

    # Pass calculated acceleration components for time $t + $dt through to next loop,
    # so that we don't have to calculate them again when we get there.
    $$numbat{'d2x'} = $next_d2x;
    $$numbat{'d2y'} = $next_d2y;
    $$numbat{'d2z'} = $next_d2z;

    # Next time step
    $t += $dt;

    # See if spline needs regeneration
    if ($use_jpl_horizons_data) {
        $spline{'moon'}->{'rebuild'} = 1 if (binsearch($spline{'moon'}->{'timestamps'}, $t) > 1);
    }

    # Limits exceeded?
    if ($earth_distance > $max_earth_distance) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        printf("%20s     Note: Earth distance %d has exceeded maximum allowed value (%d).  Stopping now.\n", $date, $earth_distance, $max_earth_distance);
        exit(0);
    } elsif ($earth_distance < $min_earth_distance) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        printf("%20s     Note: I think you've hit the Earth.  Earth distance %d has fallen below allowed value (%d).  Stopping now.\n", $date, $earth_distance, $min_earth_distance);
        exit(0);
    } elsif ($moon_distance < $Rm) {
        my $date = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($t));
        printf("%20s     Note: I think you've hit the Moon.  Moon distance %d has fallen below allowed value (%d).  Stopping now.\n", $date, $moon_distance, $Rm);
        exit(0);
    }
}
