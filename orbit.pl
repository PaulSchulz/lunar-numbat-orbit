#!/usr/bin/perl
# -*- cperl -*-
#
# Copyright 2009 Paul Schulz <paul@mawsonlakes.org>
#
# This file is part of Orbit.
#
# Orbit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Orbit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Orbit.  If not, see <http://www.gnu.org/licenses/>.

$|=1;

use Orbit;

# Reference: http://en.wikipedia.org/wiki/Earth
my $earth_radius = 6371.0;        # km             - Equatorial radius
my $earth_mass   = 5.9736e24;     # kg
my $earth_ev     = 11.186;        # km/s

# Reference: http://en.wikipedia.org/wiki/Moon
my $moon_radius  = 1737.10;       # km             - Mean radius
my $moon_mass    = 7.3477e22;     # kg
my $moon_perigee = 363104;        # km
my $moon_apogee  = 405696;        # km
my $moon_avgspeed = 1.022;        # km/s

# Reference: http://en.wikipedia.org/wiki/Gravitational_constant
my $G            = 6.67428e-11;   # m^3 kg^-1 s^-1 - Gravitational constant

my $x            = $earth_radius;
my $y            = 0.0;

my $vx           = 0.0;
my $vy           = 0.0;

my $t            = 0.0;
my $dt           = 1.0;           # s

my $ship_mass    = 100;           # kg             - Example

my $scenario = [
		{
		 'description' => 'Escape velocity',
		 'rtc' => '1234567890',
		 't' => 0.0,
		 'x' => $earth_radius + 100, # Edge of space
		 'y' => 0.0,
		 'vx' => 0.0,
		 'vy' => $earth_ev,
		},
		{
		 'description' => 'In orbit with the moon',
		 'rtc' => '1234567890',
		 't' => 0.0,
		 'x' => $moon_perigee,
		 'y' => 0.0,
		 'vx' => 0.0,
		 'vy' => $moon_avgspeed,
		},
	       ];


sub initialise ( $ ) {
  my ($scen) = @_;

  $rtc = $scenario->[$scen]->{'rtc'};
  $t = $scenario->[$scen]->{'t'};
  $x = $scenario->[$scen]->{'x'};
  $y = $scenario->[$scen]->{'y'};
  $vx = $scenario->[$scen]->{'vx'};
  $vy = $scenario->[$scen]->{'vy'};

}

sub telemetry ( ) {

  printf( "%s:%s\n","rtc", $rtc );
  printf( "%s:%s\n","t", $t );
  printf( "%s:%s\n","x", $x );
  printf( "%s:%s\n","y", $y );
  printf( "%s:%s\n","vx", $vx );
  printf( "%s:%s\n","vy", $vy );
  printf( "\n" );

};

# Calculated in kg * m / s^2
sub gravity ( $ $ ) {

  my ( $x, $y ) = @_;
  my ( $radius, $theta ) = rec2pol( $x, $y );

  $radius = $radius * 1000;  # Convert km to m;

  my $constant = $G * $earth_mass * $ship_mass;
  my $fmag = -1.0 * $constant / $radius**2 ;

  my ( $fx, $fy ) = pol2rec( $fmag, $theta );

  return ( $fx, $fy );
}

initialise(1);

# The following assumes the time interval of 1 sec.
while ( 1 ) {

  my ($fx,$fy) = gravity( $x, $y );

  # a = F /m
  $vx += $fx / $ship_mass / 1000;
  $vy += $fy / $ship_mass / 1000;

  $x += $vx;
  $y += $vy;

  telemetry();

  $t += $dt;

  sleep(1);
}
