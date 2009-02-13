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

initscr();
clear();
noecho();
cbreak();

sub template ( ) {

  move(1,2);
  printw( "time(s)   :" );

  move(2,2);
  printw( "x   (km)  :" );

  move(3,2);
  printw( "y   (km)  :" );

  move(4,2);
  printw( "v_x (km/h):" );

  move(5,2);
  printw( "v_y (km/h):" );

  refresh();

}

template();

my %data;

while ( $line = <> ) {
  chomp $line;
  if ( $line ne '' ){
    my ($key,$value) = split(':',$line,2);

    $data{$key} = $value;

  } else {

    move(1,14);
    printw("%s", sprintf("%11.3f",$data{'t'}) );

    move(2,14);
    printw("%s", sprintf("%11.3f",$data{'x'}) );

    move(3,14);
    printw("%s", sprintf("%11.3f",$data{'y'}) );

    move(4,14);
    printw("%s", sprintf("%11.3f",$data{'vx'}) );

    move(5,14);
    printw("%s", sprintf("%11.3f",$data{'vy'}) );

    refresh();
  }


}

# sub quit __END__ {
#   endwin();
# }
