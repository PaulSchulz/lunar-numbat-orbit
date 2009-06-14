#!/usr/bin/perl -w

# Filter to convert outout from Roy Duncan's orbit program into tagged
# telemetry data.

while ( my $line = <> ){
    chomp $line;

    if ( $line =~ /^UNIX ctime/ || $line =~ /Note/ ){
    } else {
	# Split the line data on whitespace.
	@data = split( /\s+/, $line );

	printf( "rtc:%f\n",
		$data[0] );
	printf( "x:%e\ny:%e\nz:%e\n",
		$data[4],
		$data[5],
		$data[6] );
	printf( "vx:%e\nvy:%e\nvz:%e\n",
		$data[7],
		$data[8],
		$data[9] );
	printf( "ax:%e\nay:%e\naz:%e\n",
		$data[10],
		$data[11],
		$data[12] );

	print "\n";
    }
    
}
	       
