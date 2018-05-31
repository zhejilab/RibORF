if ($#ARGV < 2) {
	print "##Annotate potential ORF"."\n";
	print "usage: perl ORFcandidate.pl -g genomeSequenceFile -t transcriptomeFile -o outputDir [-s startCodon] [-l orfLengthCutoff]"."\n";
	print "-g genomeSequenceFile: input the genome assembly file in fasta format;"."\n";
	print "-t transcriptomeFile: reference transcriptome annotation file in genePred format;"."\n";
	print "-o outputDir: output directory;"."\n";
	print "-s startCodon [optional]: start codon types, default: ATG;CTG;GTG;TTG;ACG;"."\n";
	print "-l orfLengthCutoff [optional]: cutoff of minimum candidate ORF length, default: 6."."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("gtosl",\%args);
my $genome=$args{g}; 
if (! $genome) {
	print "No reference genome file"."\n";
    exit;
}
my $genepred=$args{t}; 
if (! $genepred) {
	print "No reference transcriptome annotation file"."\n";
    exit;
}
my $outdir=$args{o}; 
if (! $outdir) {
	print "No output directory"."\n";
    exit;
}
my $scodon=$args{s}; 
my $lenoff=$args{l};

($scodon="ATG/CTG/GTG/TTG/ACG") if (!$scodon);
($lenoff=6) if (! $lenoff);

my @s=split /\//, $scodon;
my %sc;
for (my $i=0; $i<=$#s; $i++) {
	$sc{$s[$i]}=1;
}

### read genome file
my %seq;
my $id;
open (SE, "$genome");
while (<SE>) {
	chomp;
	if ($_ =~ /^>/) {
		my @a=split /\s+/, (substr ($_, 1));
		$id=$a[0];
	} else {
		$seq{$id}.=$_;
	}
}
close SE;

### output candidate ORF
open (IN, "$genepred");
open (OUT1, ">$outdir/candidateORF.genepred.txt");
open (OUT2, ">$outdir/candidateORF.fa");
my $line;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my @val;
	my $out;
	my $len=0;
	my $loc1=0;
	my $loc2=0;
	my $k=$s1[0].":".$s1[1].":".$s1[2];
	for (my $i=0; $i<=$#s2; $i++) {
		$out.=substr($seq{$s1[1]}, $s2[$i], ($s3[$i]-$s2[$i]));
		$len+=$s3[$i]-$s2[$i];
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			push @val, $j;
			if ($j == $s1[5]) {
				$loc1=$#val+1;
			}
			if ($j == $s1[6]) {
				$loc2=$#val+1;
			}
		}
	}
	$out=uc($out);
	if ($s1[2] eq '-') {
		$out=reverse($out);
		$out =~ tr/ACGT/TGCA/;
		@val = reverse (@val);
		my $a=$loc1;
		$loc1=$len-$loc2+2;
		$loc2=$len-$a+2;
	}
	if ($s1[5] == $s1[6]) {
		$loc1=0;
		$loc2=0;
	}
	
	my $ra=0;
	for (my $i=0; $i<=($len-2); $i++) {
		my $b=substr ($out, $i, 3);
		if (exists ($sc{$b})) {
			my $ind=0; 
			for (my $j=$i; $j<=$len && $ind==0; $j+=3) {
				$d=substr ($out, $j, 3);
				if ($d eq 'TAG' || $d eq 'TAA' || $d eq 'TGA') {
					if (($j-$i+3) >= $lenoff) {
					my $seq=substr($out, $i, $j-$i+3);
					$ra++;
					my $type;
					my $loc3=$i+1;
					my $loc4=$j+4;
					if ($loc1==0) {
						$type="noncoding";
					} elsif ($loc3==$loc1 && $loc4==$loc2) {
						$type="canonical";
					} elsif ($loc3 < $loc1 && $loc4 < $loc1) {
						$type="uORF";
					} elsif ($loc3 < $loc1 && $loc4 >= $loc1 && $loc4 < $loc2) {
						$type="overlap.uORF";
					} elsif ($loc3 >= $loc1 && $loc4 < $loc2) {
						$type="internal";
					} elsif ($loc3 > $loc1 && $loc3 < $loc2 && $loc4 > $loc2) {
						$type="external";
					} elsif ($loc3 > $loc1 && $loc4 == $loc2) {
						$type="truncation";
					} elsif ($loc3 < $loc1 && $loc4 == $loc2) {
						$type="extension";
					} elsif ($loc3 >= $loc2) {
						$type="polycistronic";
					} elsif ($loc3 == $loc1 && $loc4 != $loc2) {
						$type="seqerror";
					} elsif ($loc3 <= $loc1 && $loc4 > $loc2) {
						$type="readthrough";
					} else {
						$type="other";
					}
					my $id=$k."|".$ra."|".$len.":".($i+1).":".($j+4)."|".$type."|".$b;
					print OUT2 ">".$id."\n".$seq."\n";
					my $out2=$id."\t".$s1[1]."\t".$s1[2]."\t".$s1[3]."\t".$s1[4];
					if ($s1[2] eq '+') {
						$out2.="\t".$val[$i]."\t".$val[$j+3];
					} else {
						$out2.="\t".$val[$j+2]."\t".$val[$i-1];
					} 
					$out2.="\t".$s1[7]."\t".$s1[8]."\t".$s1[9];
					print OUT1 $out2."\n";
					}
					$ind=1;					
				}
			}
		}
	}
}
close IN;
close OUT1;
close OUT2;





