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

my %data;

while ( $line = <> ) {
  chomp $line;
  if ( $line ne '' ){
    my ($key,$value) = split(':',$line,2);

    $data{$key} = $value;

  } else {
    $key='t';  printf( "%s:%-12s ", $key,$data{$key} );
    $key='x';  printf( "%s:%-12.3f ", $key,$data{$key} );
    $key='y';  printf( "%s:%-12.3f ", $key,$data{$key} );
    $key='vx'; printf( "%s:%-12.3f ", $key,$data{$key} );
    $key='vy'; printf( "%s:%-12.3f ", $key,$data{$key} );

    print( "\n" );
  }


}

