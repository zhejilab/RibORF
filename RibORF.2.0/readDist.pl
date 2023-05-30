if ($#ARGV < 2) {
	print "## Calcualte read distribution around CDS start and stop codons"."\n";
	print "usage: perl readDist.pl -f readFile -g geneFile -o outputDir [-d readLength] [-l leftNum] [-r rightNum]"."\n";
	print "-f readFile: read alignments to the transcriptome and reference genome, in SAM format;"."\n";
	print "-g geneFile: ORF annotation file in genepred format;"."\n";
	print "-o outputDir: output directory;"."\n";
	print "-d readLength [optional]: specified RPF length, default: 25,26,27,28,29,30,31,32,33,34,;"."\n";
	print "-l leftNum [optional]: N nucleotides upstream start codon and downstream stop codon, default: 30;"."\n";
	print "-r rightNum [optional]: N nucleotides downstream start codon and upstream stop codon, default: 50."."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("fgodlr",\%args);

my $readfile=$args{f};
if (! $readfile) {
	print "No read mapping file"."\n";
    exit;
}
my $geneFile=$args{g};
if (! $geneFile) {
	print "No ORF annotation file"."\n";
    exit;
}
my $outputDir=$args{o};
if (! $geneFile) {
	print "No output directory"."\n";
    exit;
}
my $lenDist=$args{d};
my $start=$args{l};
my $end=$args{r};
if (! $lenDist) {
	$lenDist="18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,";
}
my %len;
my @s=split /,/, $lenDist;
for (my $i=0; $i<=$#s; $i++) {
	$len{$s[$i]}=1;
}
if (! $start) {
	$start=30;
}
if (! $end) {
	$end=50;
}

my %sta;

foreach my $readLength (sort keys %len) {
open (IN, "$readfile");
my %read;
my $tot=0;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\D+/, $s1[5];
	my @s3=split /\d+/, $s1[5];
	if (length($s1[9]) == $readLength) {
	if ($s1[1] == 0) {
		my $k=$s1[2].":"."+".":".($s1[3]-1);
		$read{$k}++;
		$tot++;
	} elsif ($s1[1] == 16) {
		my $loc=$s1[3];
		for (my $i=0; $i<$#s3; $i++) {
			$loc+=$s2[$i];
		}
		my $k=$s1[2].":"."-".":".($loc-2);
		$read{$k}++;
		$tot++;
	}
	}
}
close IN;
$sta{$readLength}[0]=$tot;

my %count5;
my %count3;
open (AN, "$geneFile");
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[$#s1-1];
	my @s3=split /,/, $s1[$#s1];
	my @val;
	my @ex;
	my $k=-1;
	my $codon5;
	my $codon3;
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			$k++;
			push @val, $j;
			if (exists ($read{$j})) {
				push @ex, $read{$j};
			} else {
				push @ex, "0";
			} 
			if ($s1[2] eq '+') {
				if ($j==$s1[5]) {
					$codon5=$k;
				} elsif ($j==($s1[6]-3)) {
					$codon3=$k;
				}
			} else {
				if ($j==($s1[5]+2)) {
					$codon5=$k;
				} elsif ($j==($s1[6]-1)) {
					$codon3=$k;
				}
			} 
		}
	}
	#if (($codon3-$codon5) > $end) {
	if ($codon5 > $start && ($codon3-$codon5) > $end && ($#val-$codon3) > $start) {
		$count5{"0"}++;
		$count3{"0"}++;
		for (my $i=-$start; $i<=$end; $i++) {
			my $k;
			if ($s1[2] eq '+') {
				$k=$s1[1].":".$s1[2].":".$val[$codon5+$i];
			} else {
				$k=$s1[1].":".$s1[2].":".$val[$codon3-$i];
			}
			if (exists ($read{$k})) {
				$count5{$i+$start+1}+=$read{$k};
			}
		} 
		for (my $i=-$end; $i<=$start; $i++) {
			my $k;
			if ($s1[2] eq '+') {
				$k=$s1[1].":".$s1[2].":".$val[$codon3+$i];
			} else {
				$k=$s1[1].":".$s1[2].":".$val[$codon5-$i];
			}
			if (exists ($read{$k})) {
				$count3{$i+$end+1}+=$read{$k};
			}
		} 
	}
}
close AN;

if (! exists ($count5{"0"})) {
	print "Error: no hits!"."\n";
	exit;
}

open (OUT, ">$outputDir/read.dist.sample.$readLength.txt");
my $out="start.codon";
for (my $i=1; $i<=($end+$start+1); $i++) {
	if (! exists ($count5{$i})) {
		$count5{$i}=0;
	}
	$out.="\t".sprintf("%.3e", $count5{$i}/$count5{"0"}/$tot*1000000);
}
print OUT $out."\n";

my $out="stop.codon";
for (my $i=1; $i<=($end+$start+1); $i++) {
	if (! exists ($count3{$i})) {
		$count3{$i}=0;
	}
	$out.="\t".sprintf("%.3e", $count3{$i}/$count5{"0"}/$tot*1000000);
}
print OUT $out."\n";
close OUT;

for (my $i=($end+$start+1); $i>=$start; $i-=3) {
	my $k1=$count5{$i};
	my $k2=$count5{$i-1};
	my $k3=$count5{$i-2};
	my $tk=$k1+$k2+$k3;
	if ($tk==0) {
		$tk=3;
		$k1=1;
		$k2=1;
		$k3=1;
	}
	$sta{$readLength}[1]+=$k1/$tk;
	$sta{$readLength}[2]+=$k2/$tk;
	$sta{$readLength}[3]+=$k3/$tk;
}

open (RI, ">$outputDir/readDist.plot.$readLength.R");
print RI "A <- read.table (\"$outputDir/read.dist.sample.$readLength.txt\", sep=\"\\t\")"."\n";
print RI "loc1 <- c(-$start:$end)"."\n";
print RI "loc2 <- c(-$end:$start)"."\n";
print RI "B1 <- apply(A[1,2:ncol(A)], 2, sum)"."\n";
print RI "B2 <- apply(A[2,2:ncol(A)], 2, sum)"."\n";
print RI "pdf (file=\"$outputDir/plot.readDist.$readLength.pdf\")"."\n";
print RI "par(mfrow=c(2,2), mar=c(2,2,2,1))"."\n";
print RI "plot(loc1, B1, xlim=c(-20, 30),  type=\"h\", lwd=2, ylim=c(0,max(c(B1,B2))), main=$readLength, xlab=\"Flanking start codon\")"."\n";
print RI "abline (v=-20, col=8, lwd=1, lty=\"dotted\")"."\n";
print RI "abline (v=-10, col=8, lwd=1, lty=\"dotted\")"."\n";
print RI "abline (v=0, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "abline (v=10, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "abline (v=20, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "abline (v=30, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "plot(loc2, B2, xlim=c(-30, 20), type=\"h\", lwd=2, ylim=c(0,max(c(B1,B2))), main=$readLength, xlab=\"Flanking stop codon\")"."\n";
print RI "abline (v=-30, col=8, lwd=1, lty=\"dotted\")"."\n";
print RI "abline (v=-20, col=8, lwd=1, lty=\"dotted\")"."\n";
print RI "abline (v=-10, col=8, lwd=1, lty=\"dotted\")"."\n";
print RI "abline (v=0, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "abline (v=10, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "abline (v=20, col=8, lwd=1,  lty=\"dotted\")"."\n";
print RI "dev.off()"."\n";
close RI;

my $comm="Rscript $outputDir/readDist.plot.$readLength.R";
system ($comm);
}

open (ST, ">$outputDir/sta.read.dist.$lenDist.txt");
print ST "fragment.length"."\t"."read.num"."\t"."frame2"."\t"."frame3"."\t"."frame1"."\n";
foreach my $k (sort keys %sta) {
	my $su=$sta{$k}[1]+$sta{$k}[2]+$sta{$k}[3];
	if ($su==0) {
		$su=1;
	}
	print ST $k."\t".$sta{$k}[0]."\t".sprintf("%.5f", $sta{$k}[1]/$su)."\t".sprintf("%.5f", $sta{$k}[2]/$su)."\t".sprintf("%.5f", $sta{$k}[3]/$su)."\n";
}
close ST;

open (RE, ">$outputDir/sta.readDist.plot.$lenDist.R");
print RE "A <- read.table (\"$outputDir/sta.read.dist.$lenDist.txt\", sep=\"\\t\", header=T)"."\n";
print RE "pdf (file=\"$outputDir/sta.readDist.plot.$lenDist.pdf\")"."\n";
print RE "barplot(A[,2], names.arg = A[,1], beside=T, ylab=\"Read number\")"."\n";
print RE "barplot(as.matrix(t(A[,3:5])), names.arg = A[,1], beside=T, ylim=c(0,1), legend=T)"."\n";
print RE "abline(h=0.5, lty=\"dotted\")"."\n";
print RE "dev.off()"."\n";
close RE;

my $comm="Rscript $outputDir/sta.readDist.plot.$lenDist.R";
system ($comm);

