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

# Reference: http://en.wikipedia.org/wiki/Earth
my $earth_radius = 6371.0;        # km             - Equatorial radius
my $earth_mass   = 5.9736e24;     # kg

# Reference: http://en.wikipedia.org/wiki/Gravitational_constant
my $G            = 6.67428e-11;   # m^3 kg^-1 s^-1 - Gravitational constant

my $x            = $earth_radius;
my $y            = 0.0;

my $vx           = 0.0;
my $vy           = 0.0;

my $t            = 0.0;
my $dt           = 1.0;           # s

my $ev = 11.186;  # km/h - Escape velocity


my $scenario = [
		{
		 'rtc' => '1234529684',
		 't' => 0.0,
		 'x' => $earth_radius + 100, # Edge of space
		 'y' => 0.0,
		 'vx' => 0.0,
		 'vy' => $ev,
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


sub gravity ( $ $ ) {

  my ( $x, $y ) = @_;

  my $constant = $G * $earth_mass;

#  $fx = $constant  / ( $x * $x + $y * $y );
#  $fy = $constant  / ( $x * $x + $y * $y );

  $fx = 0.0;
  $fy = 0.0;

  return ( $fx, $fy );
}

initialise(0);

while ( 1 ) {

  my ($fx,$fy) = gravity($x,$y);

  $x += $vx;
  $y += $vy;

  telemetry();

  $t += $dt;

  sleep(1);
}
