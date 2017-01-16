#!/bin/env perl

$ARGV[0] =~ /l164b(\d+)h(\d+).info/;
$beta = $1/100;
$h = $2/100;
#    print "$beta $h\n";
$ARGV[0] =~ /(l164b\d\d\dh\d\d).info/;
$latstub = $1;


while(<>) {
    if(/gauge.ssplaq = (.*)/) {	$ssplaq = $1; }
    if(/gauge.stplaq = (.*)/) {	$stplaq = $1; }
}

#print $ARGV[0];


$u0 = sqrt(sqrt(0.333333333 * 0.5 * ($ssplaq + $stplaq)));

print "$latstub $beta $h $ssplaq $stplaq $u0\n"; 
