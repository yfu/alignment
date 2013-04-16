#!/usr/bin/perl
#
# perl alignment.pl [query.txt] [template.txt]
#
# $gapo, a non-negative number, is the NEGATIVE of gap opening penalty
# $gape, a non-negative number, is the NEGATIVE of gap extension penalty
#

$gapo = 10;
$gape = 1;
open(IN, "blosum62.txt") or die;
chomp($line = <IN>);
($a, @aas) = split /\s+/, $line;
while (chomp($line = <IN>)) {
    ($aa, @vals) = split /\s+/, $line;
    for ($i=0; $i<@vals; $i++) {
	$s{$aa}{$aas[$i]} = $vals[$i];
    }
}
close(IN);

open(IN, "$ARGV[0]") or die;
chomp($query = <IN>);
close(IN);
$query =~ s/\s+//g;

open(IN, "$ARGV[1]") or die;
chomp($template = <IN>);
close(IN);
$template =~ s/\s+//g;

#print "$query\n";
#print "$template\n";

@query = split //, $query;
@template = split //, $template;

# global alignment
for ($i=0; $i<@query; $i++) {
    $m[$i][0] = $s{$query[$i]}{$template[0]};
    $pi[$i][0] = 0;
    $pj[$i][0] = 0;
}

for ($j=0; $j<@template; $j++) {
    $m[0][$j] = $s{$query[0]}{$template[$j]};
    $pi[0][$j] = 0;
    $pj[0][$j] = 0;
}

for ($i=1; $i<@query; $i++) {
#    print "row: $i\n";
    for ($j=1; $j<@template; $j++) {
	$maxm = $m[$i-1][$j-1] + $s{$query[$i]}{$template[$j]};
	$maxi = $i-1;
	$maxj = $j-1;
	for ($k=0; $k<$i-1; $k++) {
	    $m = $m[$k][$j-1] + $s{$query[$i]}{$template[$j]};
	    $m -= $gapo + $gape * ($i-$k-2);
	    if ($m > $maxm) {
		$maxm = $m;
		$maxi = $k;
		$maxj = $j-1;
	    }
	}
	for ($k=0; $k<$j-1; $k++) {
	    $m = $m[$i-1][$k] + $s{$query[$i]}{$template[$j]};
	    $m -= $gapo + $gape * ($j-$k-2);
	    if ($m > $maxm) {
		$maxm = $m;
		$maxi = $i-1;
		$maxj = $k;
	    }
	}
	$m[$i][$j] = $maxm;
	$pi[$i][$j] = $maxi;
	$pj[$i][$j] = $maxj;
    }
}

# for ($j=0; $j<@template; $j++) {
#    for ($i=0; $i<@query; $i++) {
#    	printf "%2d(%1d,%1d) ", $m[$i][$j], $pi[$i][$j], $pj[$i][$j];
#    }
#    print "\n";
# }


# identify maximum score
$maxm = $m[@query-1][@template-1];
$pii = @query-1;
$pjj = @template-1;
for ($i=0; $i<@query-1; $i++) {
    if ($m[$i][@template-1] > $maxm) {
	$maxm = $m[$i][@template-1];
	$pii = $i;
	$pjj = @template-1;
    }
}

for ($j=0; $j<@template-1; $j++) {
    if ($m[@query-1][$j] > $maxm) {
	$maxm = $m[@query-1][$j];
	$pii = @query-1;
	$pjj = $j;
    }
}

print "max score: $maxm\n";

# trace back
$ti = 0;
$tj = 0;
for ($i=@query-1; $i>$pii; $i--) {
    $aln1[$ti++] = $query[$i];
    $aln2[$tj++] = '-';
}
for ($j=@template-1; $j>$pjj; $j--) {
    $aln1[$ti++] = '-';
    $aln2[$tj++] = $template[$j];
}

$aln1[$ti++] = $query[$pii];
$aln2[$tj++] = $template[$pjj];
    
while ($pii > 0 && $pjj > 0) {
    for ($i=$pii-1; $i>$pi[$pii][$pjj]; $i--) {
	$aln1[$ti++] = $query[$i];
	$aln2[$tj++] = '-';
    }
    for ($j=$pjj-1; $j>$pj[$pii][$pjj]; $j--) {
	$aln1[$ti++] = '-';
	$aln2[$tj++] = $template[$j];
    }
    $pii2 = $pi[$pii][$pjj];
    $pjj2 = $pj[$pii][$pjj];
    $pii = $pii2;
    $pjj = $pjj2;
    $aln1[$ti++] = $query[$pii];
    $aln2[$tj++] = $template[$pjj];
}

for ($i=$pii-1; $i>=0; $i--) {
    $aln1[$ti++] = $query[$i];
    $aln2[$tj++] = '-';
}
for ($j=$pjj-1; $j>=0; $j--) {
    $aln1[$ti++] = '-';
    $aln2[$tj++] = $template[$j];
}

open(OUT, ">alignment.txt") or die;
for ($i=@aln1-1; $i>=0; $i--) {
    print OUT "$aln1[$i]";
    print "$aln1[$i]";
}
print OUT "\n";
print "\n";
for ($j=@aln2-1; $j>=0; $j--) {
    print OUT "$aln2[$j]";
    print "$aln2[$j]";
}
print OUT "\n";
print "\n";
close(OUT);
