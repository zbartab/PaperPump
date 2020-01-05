#! /usr/bin/perl -w
#
# a script to sort out article entries from the dblp-2019-12-01.xml
# publication database
#
# only those articles are written which have a doi.org entry

my $store = 0;
#my $lineno = 1;
my $record = "";
my $print = 0;

while (my $line = <>) {
	#if ($lineno <= 3) {
	#print $line;
	#$lineno++;
	#}
	if ($line =~ /^<article/) {
		$record = "";
		$store = 1;
	}
	$record .= $line if ($store);
	$print++ if ($store and $line =~ m+<ee>https://doi\.org+);
	if ($line =~ m+^</article+) {
		if ($print > 0) {
			my $ee = "";
			my $year = "";
			my @authors;
			$record =~ s/&(.)[a-zA-Z]+;/$1/g;
			my @entities = split "\n", $record;
			foreach my $rec (@entities) {
				if ($rec =~ m+<ee>https://doi\.org+) {
					$ee = $rec;
					$ee =~ s+<ee>https://(.*)</ee>+$1+;
				} elsif ($rec =~ /<year>/) {
					$year = $rec;
					$year =~ s+<year>(.*)</year>+$1+;
				} elsif ($rec =~ /<author/) {
					my $author = $rec;
					$author =~ s+<author.*>(.*)</author>+$1+;
					push @authors, $author;
				}
				#print $rec, "\n" if ($rec =~ m+</*article+ or $rec =~ /<author/
				#or $rec =~ /<year>/ or $rec =~ m+<ee>\s*https://doi\.org+);
			}
			if ($ee and $year and @authors) {
				print $ee, "\n";
				print join("; ", @authors), "\n";
				print $year, "\n";
			}
			#print $record;
		}
		$store = 0;
		$print = 0;
	}
}

#print "</dblp>\n";
