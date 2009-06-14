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

use Curses;
use Orbit;

initscr();
clear();
noecho();
cbreak();

sub template ( ) {

  move(2,2);
  printw( "  time(s)   :" );

  move(3,2);
  printw( "  x   (km)  :" );

  move(4,2);
  printw( "  y   (km)  :" );

  move(5,2);
  printw( "  v_x (km/s):" );

  move(6,2);
  printw( "  v_y (km/s):" );

  # Derived measurements
  move(8,2);
  printw( "r_mag   (km):" );

  move(9,2);
  printw( "r_theta(rad):" );

  move(10,2);
  printw( "v_mag (km/s):" );

  move(11,2);
  printw( "v_theta(rad):" );

  refresh();

}

# Display parameter
my $dx1 = 30;
my $dx2 = 79;
my $dx =  int($dx1 + ($dx2 - $dx1) / 2);

my $dy1 = 1;
my $dy2 = 23;
my $dy =  int($dy1 + ($dy2 - $dy1) / 2);

sub near_earth_orbit_display ( ) {

# This image would need to be adjusted if the display parameters
# are changed.

  move($dy-9,$dx-21);
  printw("              ,------------.              ");
  move($dy-8,$dx-21);
  printw("          .--'              `--.          ");
  move($dy-7,$dx-21);
  printw("       .-'                      `-.       ");
  move($dy-6,$dx-21);
  printw("     .'                            `.     ");
  move($dy-5,$dx-21);
  printw("   .'                                `.   ");
  move($dy-4,$dx-21);
  printw("  '                                    `  ");
  move($dy-3,$dx-21);
  printw(" :                                      : ");
  move($dy-2,$dx-21);
  printw(" :                                      : ");
  move($dy-1,$dx-21);
  printw(":                                        :");
  move($dy,$dx-21);
  printw(":                                        :");
  move($dy+1,$dx-21);
  printw(" :                                      : ");
  move($dy+2,$dx-21);
  printw(" :                                      : ");
  move($dy+3,$dx-21);
  printw("  `                                    '  ");
  move($dy+4,$dx-21);
  printw("   `.                                ,'   ");
  move($dy+5,$dx-21);
  printw("     `.                            ,'     ");
  move($dy+6,$dx-21);
  printw("       `-.                      ,-'       ");
  move($dy+7,$dx-21);
  printw("          `--.              .--'          ");
  move($dy+8,$dx-21);
  printw("              `------------'              ");

  refresh();

  return $earth_radius+1000;

}

sub earth_orbit_display ( ) {

  move($dy-1,$dx-2);
  printw(" .-.");
  move($dy,$dx-2);
  printw("(   )");
  move($dy+1,$dx-2);
  printw(" `-'");

  #    move(23,79);
  #    printw(".");

  refresh();

  return $earth_radius+6000;

}

sub earth_moon_display ( ) {

  move($dy,$dx);
  printw("O");

  #    move(23,79);
  #    printw(".");

  refresh();

  return $moon_apogee;

}

template();

my %data;

my $scale = near_earth_orbit_display();

my $X=0;
my $Y=0;

while ( $line = <> ) {
  chomp $line;
  if ( $line ne '' ){
    my ($key,$value) = split(':',$line,2);

    $data{$key} = $value;

  } else {

    move(2,16);
    printw("%s", sprintf("%11s",sec2dur(int($data{'t'}))) );

    move(3,16);
    printw("%s", sprintf("%11.3f",$data{'x'}) );

    move(4,16);
    printw("%s", sprintf("%11.3f",$data{'y'}) );

    move(5,16);
    printw("%s", sprintf("%11.3f",$data{'vx'}) );

    move(6,16);
    printw("%s", sprintf("%11.3f",$data{'vy'}) );

    # Derived values
    my $x = $data{'x'};
    my $y = $data{'y'};

    my ( $r_mag, $r_theta ) = rec2pol ( $x, $y );

    move(8,16);
    printw("%s", sprintf("%11.3f",$r_mag) );

    move(9,16);
    printw("%s", sprintf("%11.3f",$r_theta) );

    my $vx = $data{'vx'};
    my $vy = $data{'vy'};

    my $v_mag   = sqrt ( $vx*$vx + $vy*$vy );
    my $v_theta = atan2( $vy, $vx );

    move(10,16);
    printw("%s", sprintf("%11.3f",$v_mag) );

    move(11,16);
    printw("%s", sprintf("%11.3f",$v_theta) );

    move($Y,$X);
    printw(" ");

    $X = $x / $scale * ( $dx2 - $dx1 )/2.0 + $dx;
    $Y = $y / $scale * ( $dy2 - $dy1 )/2.0 + $dy;

    move( $Y, $X );
    printw("+");

    refresh();
  }


}

# sub quit __END__ {
#   endwin();
# }
