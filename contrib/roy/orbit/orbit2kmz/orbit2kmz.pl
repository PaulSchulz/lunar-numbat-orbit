#!/usr/bin/perl -w

###
### Convert orbit.pl output to a KMZ archive, suitable for
### display in Google Earth.
###
### Usage is something like:
###
###     ./orbit.pl > output.txt
###     cd orbit2kmz
###     cat ../output.txt | ./orbit2kmz.pl
###
### This will make a KMZ archive named "orbit.kmz", which should contain
### both Numbat and Moon positional data suitable for viewing in Google Earth.
###
### More doco coming...
###
###
### HISTORY
###
### [01-Jun-2009 Roy] Initial coding.
###
###
###
###
use strict;
use Math::Trig;


my $zip = 'zip';
my $kmlfile = "doc.kml";
my $zipfile = "orbit.kmz";

# Earth axial tilt (rad)
my $earth_axial_tilt = deg2rad(23.43928);

# Need planechange() to make RA and Dec of the Moon.
sub planechange {

    my ($x, $y, $z, $T) = @_;

    my $new_x = $x;
    my $new_y = cos($T) * $y - sin($T) * $z;
    my $new_z = sin($T) * $y + cos($T) * $z;
    return ($new_x, $new_y, $new_z);
}


# Open the KML output file
unlink "$kmlfile.new" if (-e "$kmlfile.new");
open(KML, ">$kmlfile.new") or die "Can't open output file $kmlfile.new : $!";

# date for KML
my $kml_datestring_full = scalar gmtime(time());


print KML <<EOF;
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated on $kml_datestring_full UTC by lunarnumbat.org       -->
<kml xmlns="http://earth.google.com/kml/2.2" hint="target=sky"
     xmlns:atom="http://www.w3.org/2005/Atom">
<Document>	
        <name>Position of Lunar Numbat</name>
	<atom:author>      
		<atom:name>lunarnumbat.org</atom:name>
	</atom:author>
	<atom:link href="http://www.lunarnumbat.org/" />

	<open>0</open>
	<StyleMap id="NumbatStyle">
		<Pair>
			<key>normal</key>
			<styleUrl>#NumbatNormal</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#NumbatHighlight</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="NumbatHighlight">
		<IconStyle>
			<scale>1.5</scale>
			<Icon>
				<href>files/numbat.png</href>
			</Icon>
		</IconStyle>
		<BalloonStyle>
			<text>&lt;center&gt;&lt;font face=&apos;Verdana&apos; size=&apos;+2&apos; color=&apos;#64634c&apos;&gt;&lt;b&gt;\$[name]&lt;/b&gt;&lt;/font&gt;&lt;/center&gt;&lt;br/&gt;\$[description]</text>
		</BalloonStyle>
	</Style>
	<Style id="NumbatNormal">
		<IconStyle>
			<scale>1.5</scale>
			<Icon>
				<href>files/numbat.png</href>
			</Icon>
		</IconStyle>
		<BalloonStyle>
			<text>&lt;center&gt;&lt;font face=&apos;Verdana&apos; size=&apos;+2&apos; color=&apos;#64634c&apos;&gt;&lt;b&gt;\$[name]&lt;/b&gt;&lt;/font&gt;&lt;/center&gt;&lt;br/&gt;\$[description]</text>
		</BalloonStyle>
	</Style>
	<StyleMap id="MoonStyle">
		<Pair>
			<key>normal</key>
			<styleUrl>#MoonNormal</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#MoonHighlight</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="MoonHighlight">
		<IconStyle>
			<scale>1.5</scale>
			<Icon>
				<href>files/moon.png</href>
			</Icon>
		</IconStyle>
		<BalloonStyle>
			<text>&lt;center&gt;&lt;font face=&apos;Verdana&apos; size=&apos;+2&apos; color=&apos;#64634c&apos;&gt;&lt;b&gt;\$[name]&lt;/b&gt;&lt;/font&gt;&lt;/center&gt;&lt;br/&gt;\$[description]</text>
		</BalloonStyle>
	</Style>
	<Style id="MoonNormal">
		<IconStyle>
			<scale>1.5</scale>
			<Icon>
				<href>files/moon.png</href>
			</Icon>
		</IconStyle>
		<BalloonStyle>
			<text>&lt;center&gt;&lt;font face=&apos;Verdana&apos; size=&apos;+2&apos; color=&apos;#64634c&apos;&gt;&lt;b&gt;\$[name]&lt;/b&gt;&lt;/font&gt;&lt;/center&gt;&lt;br/&gt;\$[description]</text>
		</BalloonStyle>
	</Style>
	<Style id="check-hide-children">        <!-- define the style for the Document -->
		<ListStyle>
			<listItemType>checkHideChildren</listItemType>
		</ListStyle>
	</Style>
	<styleUrl>#check-hide-children</styleUrl> <!-- add the style to the Document -->
EOF


# Begin looping
my @data;
my $counter = 0;
while(<STDIN>) {
    next unless ($_ =~ m/^1/);
    my @bits = split(/ +/, $_);
    my $ctime = $bits[0];
    my $ra = $bits[15];
    my $dec = $bits[16];
    my $ra_gs = $ra * 15 - 180;
    my $distance_from_earth = $bits[17];
    my $distance_from_moon  = $bits[19];
    my $speed_wrt_earth     = $bits[21];
    my $speed_wrt_moon      = $bits[23];

    # Make Moon RA and Dec
    my ($x, $y, $z) = planechange($bits[26], $bits[27], $bits[28], $earth_axial_tilt);
    my $moon_ra = atan2($y, $x);
    $moon_ra += 2 * pi if ($moon_ra < 0);
    my $moon_dec = atan($z / sqrt($x * $x + $y * $y));
    my $moon_ra_gs = rad2deg($moon_ra) - 180;
    $moon_ra = rad2deg($moon_ra / 15);
    $moon_dec = rad2deg($moon_dec);
    push(@data, "$ctime|$ra|$dec|$ra_gs|$distance_from_earth|$distance_from_moon|$speed_wrt_earth|$speed_wrt_moon|$moon_ra|$moon_dec|$moon_ra_gs");
}

###
### Loop over data
###
for (my $i = 0; $i < $#data; $i++) {
    my ($ctime, $ra, $dec, $ra_gs, $distance_from_earth, $distance_from_moon, $speed_wrt_earth, $speed_wrt_moon, $moon_ra, $moon_dec, $moon_ra_gs) = split(/\|/, $data[$i]);
    my ($next_ctime, undef) = split(/\|/, $data[$i + 1], 2);
    my ($seconds_start, $minutes_start, $hours_start, $day_start, $month_start, $year_start, undef, undef, undef) = gmtime($ctime);
    my ($seconds_finish, $minutes_finish, $hours_finish, $day_finish, $month_finish, $year_finish, undef, undef, undef) = gmtime($next_ctime);
    $year_start += 1900;
    $year_finish += 1900;
    $month_start++;
    $month_finish++;


    ###
    ### Make RA and Dec strings - there's probably a more elegant way to do this.
    ###
    my $deg = int($ra);
    my $min = int(($ra - $deg) * 60);
    my $sec = ($ra - $deg - ($min / 60)) * 3600;
    my $ra_string = sprintf("%02d<sup>h</sup> %02d<sup>m</sup> %04.1f<sup>s</sup>", $deg, $min, $sec);
    # Dec string
    my $sign = "+";
    if ($dec < 0) {
        $sign = "-";
        $dec = abs($dec);
    }
    $deg = int($dec);
    $min = int(($dec - $deg) * 60);
    $sec = int(($dec - $deg - ($min / 60)) * 3600);
    my $dec_string = sprintf("%1s%02d&#176; %02d' %02d''", $sign, $deg, $min, $sec);
    $dec *= -1.0 if ($sign eq '-');
    # For the Moon
    $deg = int($moon_ra);
    $min = int(($moon_ra - $deg) * 60);
    $sec = ($moon_ra - $deg - ($min / 60)) * 3600;
    my $moon_ra_string = sprintf("%02d<sup>h</sup> %02d<sup>m</sup> %04.1f<sup>s</sup>", $deg, $min, $sec);
    # Moon Dec string
    $sign = "+";
    if ($moon_dec < 0) {
        $sign = "-";
        $moon_dec = abs($moon_dec);
    }
    $deg = int($moon_dec);
    $min = int(($moon_dec - $deg) * 60);
    $sec = int(($moon_dec - $deg - ($min / 60)) * 3600);
    my $moon_dec_string = sprintf("%1s%02d&#176; %02d' %02d''", $sign, $deg, $min, $sec);
    $moon_dec *= -1.0 if ($sign eq '-');


    ##
    ## Write KML for Numbat and the Moon at this point.
    ##
    print KML "<Placemark>\n";
    print KML "   <TimeSpan>\n";
    printf(KML "      <begin>%4d-%02d-%02dT%02d:%02d:%02dZ</begin>\n", $year_start, $month_start, $day_start, $hours_start, $minutes_start, $seconds_start);
    printf(KML "      <end>%4d-%02d-%02dT%02d:%02d:%02dZ</end>\n", $year_finish, $month_finish, $day_finish, $hours_finish, $minutes_finish, $seconds_finish);
    print KML "   </TimeSpan>\n";
    print KML "    <name>Lunar Numbat</name>\n";
    print KML "    <description>\n";
    printf(KML "        <![CDATA[<p><table><tr><td colspan=\"2\">%04d-%02d-%02dT%02d:%02d:%02d UTC:</td></tr><tr><td>RA (J2000.0):</td><td>%s</td></tr><tr><td>Declination:</td><td>%s</td></tr><tr><td>Earth distance:</td><td>%3.1f km</td></tr><tr><td>Speed wrt Earth:</td><td>%4.2f km/s</td></tr><tr><td>Moon distance:</td><td>%4.2f km</td></tr><tr><td>Speed wrt Moon:</td><td>%4.2f km/s</td></tr></table></p>]]>\n", $year_start, $month_start, $day_start, $hours_start, $minutes_start, $seconds_start, $ra_string, $dec_string, $distance_from_earth / 1000, $speed_wrt_earth / 1000, $distance_from_moon / 1000, $speed_wrt_moon / 1000);
    print KML "    </description>\n";
    print KML "    <Point>\n";
    printf(KML "       <coordinates>%7.5f,%6.4f,0</coordinates>\n", $ra_gs, $dec);
    print KML "    </Point>\n";
    print KML "    <styleUrl>#NumbatStyle</styleUrl>\n";
    print KML "</Placemark>\n";

    # Do MOON's KML
    print KML "<Placemark>\n";
    print KML "   <TimeSpan>\n";
    printf(KML "      <begin>%4d-%02d-%02dT%02d:%02d:%02dZ</begin>\n", $year_start, $month_start, $day_start, $hours_start, $minutes_start, $seconds_start);
    printf(KML "      <end>%4d-%02d-%02dT%02d:%02d:%02dZ</end>\n", $year_finish, $month_finish, $day_finish, $hours_finish, $minutes_finish, $seconds_finish);
    print KML "   </TimeSpan>\n";
    print KML "    <name>The Moon</name>\n";
    print KML "    <description>\n";
    printf(KML "        <![CDATA[<p><table><tr><td colspan=\"2\">%04d-%02d-%02dT%02d:%02d:%02d UTC:</td></tr><tr><td>RA (J2000.0):</td><td>%s</td></tr><tr><td>Declination:</td><td>%s</td></tr><tr><td>Numbat distance:</td><td>%4.2f km</td></tr><tr><td>Speed wrt Numbat:</td><td>%4.2f km/s</td></tr></table></p>]]>\n", $year_start, $month_start, $day_start, $hours_start, $minutes_start, $seconds_start, $moon_ra_string, $moon_dec_string, $distance_from_moon / 1000, $speed_wrt_moon / 1000);
    print KML "    </description>\n";
    print KML "    <Point>\n";
    printf(KML "       <coordinates>%7.5f,%6.4f,0</coordinates>\n", $moon_ra_gs, $moon_dec);
    print KML "    </Point>\n";
    print KML "    <styleUrl>#MoonStyle</styleUrl>\n";
    print KML "</Placemark>\n";
}


# Close out KML
print KML "</Document>\n";
print KML "</kml>\n";
close(KML);
rename "$kmlfile.new", "$kmlfile";
# Zip it!  Zip it all!
unlink "$zipfile.new" if (-f "$zipfile.new");
system("$zip -r $zipfile.new doc.kml files");
rename "$zipfile.new", "$zipfile";
# Tidying up.
unlink "$kmlfile.new" if (-e "$kmlfile.new");
unlink "$kmlfile" if (-e "$kmlfile");

# Bye!
exit(0);
