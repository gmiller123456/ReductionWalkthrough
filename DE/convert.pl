#!/usr/bin/perl
use strict;

my $f;

print "let de431=[\r\n";

open($f,"ascp02000.431");
while(my $l=<$f>){
	$l=~s/\r*\n*//gis;
	if(length($l)> 20){
		$l=~s/D/E/g;
		if(length($l)<40){
			print "$l,0,0,\r\n";
		} else {
			$l=~s/(\d) /$1\,/gis;
			print "$l,\r\n";
		}
	} else {
		print "//$l\r\n";
	}
}

print "];\r\n";

close($f);