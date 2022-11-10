if ($#ARGV < 2) {
	print "##merge predicted ORFs"."\n";
	print "usage: perl mergeORF.pl -f readCorrectedFile -c predictedORFFile -o outputDir"."\n";
	print "-f readCorrectedFile: input read mapping file after offset correction;"."\n";
	print "-c predictedORFFile: predicted ORFs in genePred format;"."\n";
	print "-o outputDir: output directory;"."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("fco",\%args);
my $readfile=$args{f};
if (! $readfile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $orfFile=$args{c};
if (! $orfFile) {
	print "No predicted ORF file"."\n";
    exit;
}
my $outputDir=$args{o};
if (! $outputDir) {
	print "No output directory"."\n";
    exit;
}

my @temp=split /~/, $orfFile;
my %tog;
my $snum=0;
foreach my $file (sort @temp) {
	open (IN, "$file");
	$snum++;
	while (<IN>) {
		chomp;
		my @s=split /\t/, $_;
		my $k=$s[1];
		for (my $i=2; $i<=$#s; $i++) {
			$k.="\t".$s[$i];
		}
		$tog{$k}[0]=$s[0];
		$tog{$k}[1].=$snum.":";
	}
	close IN;
}

open (OUT, ">$outputDir/orf.tog.genepred.txt");
foreach my $k (sort keys %tog) {
	my @s=split /:/, $tog{$k}[1];
	my %num;
	for (my $i=0; $i<=$#s; $i++) {
		$num{$s[$i]}++;
	}
	print OUT $tog{$k}[0]."\t".$k."\t".keys(%num)."\n";
}
close OUT;

my %read;
open (IN, "$readfile");
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\D+/, $s1[5];
	my @s3=split /\d+/, $s1[5];
	if ($s1[1] == 0) {
		my $k=$s1[2].":"."+".":".($s1[3]-1);
		$read{$k}++;
	} elsif ($s1[1] == 16) {
		my $loc=$s1[3];
		for (my $i=0; $i<$#s3; $i++) {
			$loc+=$s2[$i];
		}
		my $k=$s1[2].":"."-".":".($loc-2);
		$read{$k}++;
	}
}
close IN;

print "done reading read file"."\n";

open (IN, "$outputDir/orf.tog.genepred.txt");
my %orf;
my %ind;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my $loc1; #start codon
	my $loc2; #stop codon
	if ($s1[2] eq '+') {
		$loc1=$s1[1].":".$s1[2].":".$s1[3];
		$loc2=$s1[1].":".$s1[2].":".$s1[4];
	} else {
		$loc1=$s1[1].":".$s1[2].":".$s1[4];
		$loc2=$s1[1].":".$s1[2].":".$s1[3];
	}
	my $k=$s1[0];
	for (my $i=1; $i<$#s1; $i++) {
		$k.="\t".$s1[$i];
	}
	$orf{$loc1}.=$k."\n";
	$orf{$loc2}.=$k."\n";
	$ind{$k}[0] = 1;
	$ind{$k}[1] = $s1[$#s1];
}
close IN;

foreach my $km (sort keys %orf) {
	my @s1=split /\n/, $orf{$km};
	if ($#s1 > 0 && $orf{$km} =~ /\|ATG/) {
		my $out;
		for (my $i=0; $i<=$#s1; $i++) {
			if ($s1[$i] =~ /\|ATG/) {
				$out.=$s1[$i]."\n";
			} else {
				$ind{$s1[$i]}[0] = 0;
			}
		}
		$orf{$km}=$out;
	}
}

### remove ORFs with unique regions but containing 0 zero reads
foreach my $km (sort keys %orf) {
	my @s1=split /\n/, $orf{$km};
	if ($#s1 > 0) {
		print $km."\n";
		my %loc;
		for (my $l=0; $l<=$#s1; $l++) {
			my @s2=split /\t/, $s1[$l];
			my @s3=split /,/, $s2[8];
			my @s4=split /,/, $s2[9];
			for (my $i=0; $i<=$#s3; $i++) {
				for (my $j=$s3[$i]; $j<$s4[$i]; $j++) {
					if ($j >=$s2[5] && $j < $s2[6]) {
						my $k=$s2[1].":".$s2[2].":".$j;
						$loc{$k}++;
					}	
				}
			}
		}
		my $out;
		for (my $l=0; $l<=$#s1; $l++) {
			my @s2=split /\t/, $s1[$l];
			my @s3=split /,/, $s2[8];
			my @s4=split /,/, $s2[9];
			my $ulen=0;
			my $uread=0;
			for (my $i=0; $i<=$#s3; $i++) {
				for (my $j=$s3[$i]; $j<$s4[$i]; $j++) {
					if ($j >=$s2[5] && $j < $s2[6]) {
						my $k=$s2[1].":".$s2[2].":".$j;
						if ($loc{$k}==1) {
							$ulen++;
							if (exists ($read{$k})) {
								$uread+=$read{$k};
							}
						}
					}
				}
			}
			if ($ulen > 0 && $uread==0) {
				$out.=$s1[$l]."\n";
			}
		}
		my @s=split /\n/, $out;
		if ($#s >=0 && $#s < $#s1) {
			for (my $l=0; $l<=$#s; $l++) {
				$ind{$s[$l]}[0] = 0;
			}
		}
	}
}

########## remove ORFs without unique regions, expect these with the same start and stop codons
foreach my $km (sort keys %orf) {
	my @a1=split /\n/, $orf{$km};
	my $a2;
	for (my $l=0; $l<=$#a1; $l++) {
		if ($ind{$a1[$l]}[0] == 1) {
			$a2.=$a1[$l]."\n";
		}
	}
	my @s1=split /\n/, $a2;
	if ($#s1 > 0) {
		print $km."\n";
		for (my $l=0; $l<=$#s1; $l++) {
			my %loc;
			my @s2=split /\t/, $s1[$l];
			my @s3=split /,/, $s2[8];
			my @s4=split /,/, $s2[9];
			for (my $i=0; $i<=$#s3; $i++) {
				for (my $j=$s3[$i]; $j<$s4[$i]; $j++) {
					if ($j >=$s2[5] && $j < $s2[6]) {
						my $k=$s2[1].":".$s2[2].":".$j;
						$loc{$k}=1;
					}
				}
			}
			for (my $m=0; $m<=$#s1; $m++) {
				if ($m != $l) {
					foreach my $k (keys %loc) {
						$loc{$k}=1;
					}
					my @b2=split /\t/, $s1[$m];
					my @b3=split /,/, $b2[8];
					my @b4=split /,/, $b2[9];
					if (($b2[5] ne $s2[5]) || ($b2[6] ne $s2[6])) {
						my $ulen=0;
						for (my $i=0; $i<=$#b3; $i++) {
							for (my $j=$b3[$i]; $j<$b4[$i]; $j++) {
								if ($j >=$b2[5] && $j < $b2[6]) {
									my $k=$b2[1].":".$b2[2].":".$j;
									if (exists ($loc{$k})) {
										$loc{$k}=0;
									}
								}
							}
						}
						foreach my $k (keys %loc) {
							if ($loc{$k}==1) {
								$ulen++;
							}
						}
						if ($ulen==0) {
							$ind{$s1[$l]}[0] = 0;
							$ind{$s1[$m]}[1] += $ind{$s1[$l]}[1];
						}
					}
				}
			}
			
		}
	}
}

open (OUT1, ">$outputDir/retained.ORF.genepred.txt");
open (OUT2, ">$outputDir/filtered.ORF.genepred.txt");
foreach my $k (sort keys %ind) {
	if ($ind{$k}[0] == 1 && $ind{$k}[1] > 1) {
		print OUT1 $k."\n";
	} else {
		print OUT2 $k."\n";
	}
}
close OUT1;
close OUT2;

