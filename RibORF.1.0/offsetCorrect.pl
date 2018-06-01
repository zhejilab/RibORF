if ($#ARGV < 2) {
	print "## Correct the offset between 5' end of fragments and A-site"."\n";
	print "usage: perl offsetCorrect.pl -f readFile -p offsetParameterFile -o readCorrectedFile"."\n";
	print "-f readFile: read mapping file before offset correction, SAM format;"."\n";
	print "-p offsetParameterFile: parameters for offset correction, 1st column: read length, 2nd column: offset distance;"."\n";
	print "-o readCorrectedFile: output file after offset correction."."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("fpo",\%args);

my $readfile=$args{f};
if (! $readfile) {
	print "No read mapping file"."\n";
    exit;
}
my $offFile=$args{p};
if (! $offFile) {
	print "No parameter file for offset correction"."\n";
    exit;
}
my $outputFile=$args{o};
if (! $outputFile) {
	print "No outpout file"."\n";
    exit;
}

my %dist;
open (AN, "$offFile");
while (<AN>) {
	chomp;
	my @s=split /\t/, $_;
	$dist{$s[0]}=$s[1];
}
close AN;

open (IN, "$readfile");
open (OUT, ">$outputFile");
while (<IN>) {
	chomp;
	if ($_ =~ /^@/) {
		print OUT $_."\n";
	} else {
		my @s1=split /\t/, $_;
		my @s2=split /\D+/, $s1[5];
		my @s3=split /\d+/, $s1[5];
		my $len=0;
		for (my $i=0; $i<$#s3; $i++) {
			if ($s3[$i+1] eq 'M') {
				$len+=$s2[$i];
			}
		}
		if (exists ($dist{$len})) {
			my $asite=$dist{$len};
			my $loc=$s1[3];
			my $ra=0;
			my $ind=1;
			if ($s1[1]==0) {
				for (my $i=0; $i<$#s3 && $ind==1; $i++) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc+=$asite;
							$ind=0;
						} else {
							$loc+=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc+=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc+=$s2[$i];
					} 
				}
			} else {
				for (my $i=0; $i<$#s3; $i++) {
					$loc+=$s2[$i];
				}
				$loc--;
				for (my $i=($#s3-1); $i>=0 && $ind==1; $i--) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc-=$asite;
							$ind=0;
						} else {
							$loc-=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc-=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc-=$s2[$i];
					}
				}
			}
			my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$loc."\t"."50"."\t"."1M"."\t".$s1[6]."\t".$s1[7]."\t".$s1[8]."\t".substr($s1[9], $asite-1, 1)."\t".substr($s1[10], $asite-1, 1);
			for (my $i=11; $i<=$#s1; $i++) {
				$out.="\t".$s1[$i];
			}
			print OUT $out."\n";
		}
	}
}
close IN;
close OUT;


