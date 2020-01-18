#! /usr/bin/perl -wn
#
# script to convert old scimat files to new scimat files
# old scimat: it just list authors' IDs and papers' IDs in order
# new scimat: parts of authors' IDs and papers' IDs contain not just the
# IDs, but their positions too

BEGIN {
	$index = 0;
}

$line = $_;

if ($line =~ m/^#### .*IDs/) {
	print $line;
	$index = 1;
} elsif ($index > 0) {
	chomp $line;
	print $line, ",", $index, "\n";
	$index++;
} else {
	print $line;
}
