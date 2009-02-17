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

# Convert rectalinear coordinates (x,y) to polar coordinates
# (r,theta).
sub rec2pol ( $ $ ) {
  my ( $x, $y ) = @_;

  my $radius = sqrt( $x * $x + $y * $y );
  my $theta = atan2 ( $y, $x );

  return ( $radius, $theta );
}

# Convert polar coordinates (r,theta) to rectalinear coordinates
# (x,y).
sub pol2rec ( $ $ ) {
  my ( $radius, $theta ) = @_;

  my $x = $radius * cos( $theta );
  my $y = $radius * sin( $theta );

  return ( $x, $y );
}

# Convert seconds into a duration string (eg. 1h4m20s)
# Assumes that the duration is positive.
sub sec2dur ( $ ) {
  my ($sec) =@_;
  my $dur = '';

  my $d = -1;
  my $h = -1;
  my $m = -1;
  my $s = -1;

  my $units=0;

  if ( $sec >= 24*60*60 )
    {
      $d = int( $sec / (24*60*60) );
      $sec = $sec-($d*24*60*60);

      $units++;
    }
  if ( $sec >= 60*60 )
    {
      $h = int( $sec / (60*60) );
      $sec = $sec-($h*60*60);

      $units++;
    }
  if ( $sec >= 60 )
    {
      $m = int( $sec / (60) );
      $sec = $sec-($m*60);

      $units++;
    }
  $s = $sec;

  if ( $d != -1 ){
    $dur .= $d."d";
  }
  if ( $h != -1 ){
    $dur .= $h."h";
  }
  if ( $m != -1 ){
    $dur .= $m."m";
  }
  if ( $s != -1 ){
    $dur .= $s."s";
  }

  if ( $dur eq '' ){
    $dur = '-';
  }

  return $dur;
}


1;
