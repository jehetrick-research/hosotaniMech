#!/usr/bin/env perl

# General purpose extract/sum/average columns of an output file

$TAG = shift;
$colx = shift;
$coly = shift;


# $TAG = 'MEASU1';
# $colx = 1;
# $coly = 6;

while(<>) {
    if(/^$TAG /) {
	@line = split();
	$x = $line[$colx];
	$y{$x} += $line[$coly];
	$count{$x}++;
    }
}

# Add std.dev later

foreach $x (sort keys %y) {
    $yave = $y{$x}/$count{$x};
    print "$x $yave\n";
#    print "$yave, ";
}
