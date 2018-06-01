if ($#ARGV < 2) {
	print "##ribORF predication"."\n";
	print "usage: perl ribORF.pl -f readCorrectedFile -c candidateORFFile -o outputDir [-l orfLengthCutoff] [-r orfReadCutoff] [-p predictPvalueCutoff]"."\n";
	print "-f readCorrectedFile: input read mapping file after offset correction;"."\n";
	print "-c candidateORFFile: candidate ORFs in genePred format;"."\n";
	print "-o outputDir: output directory, with files reporting testing parameters and predicted translating probability;"."\n";
	print "-l orfLengthCutoff [optional]: cutoff of ORF length (nt), default: 6;"."\n";
	print "-r orfReadCutoff [optional]: cutoff of supported read number, default: ;"."\n";
	print "-p predictPvalueCutoff [optional]: cutoff used to select predicted translated ORF, default: 0.7."."\n";
	exit;
}


use Getopt::Std;

### get arguments ###
my %args; 
getopt("fcolrp",\%args);
my $readfile=$args{f};
if (! $readfile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $orfFile=$args{c};
if (! $orfFile) {
	print "No candidate ORF annotation file"."\n";
    exit;
}
my $outputDir=$args{o};
if (! $outputDir) {
	print "No output directory"."\n";
    exit;
}
my $orfLengthCutoff;
my $orfReadCutoff;
my $predictPvalueCutoff;
if (exists ($args{l})) {
	$orfLengthCutoff=$args{l};
} else {
	$orfLengthCutoff=6;
}
if (exists ($args{r})) {
	$orfReadCutoff=$args{r};
} else {
	$orfReadCutoff=11;
}
if (exists ($args{p})) {
	$predictPvalueCutoff=$args{p};
} else {
	$predictPvalueCutoff=0.7;
}

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

open (AN, "$orfFile");
open (OUT, ">$outputDir/input.parameters.txt");
print OUT "orfID"."\t"."chrom"."\t"."strand"."\t"."codon5"."\t"."codon3"."\t"."length"."\t"."readNum"."\t"."f1"."\t"."f2"."\t"."f3"."\t"."entropy"."\t"."MAXentropy"."\t"."PME"."\t"."codonNum"."\t"."f1max"."\n";
my $dist=3;
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my $len1=0;
	my $tot=0;
	my @val;
	my @per;
	my @post;
	my @fra;
	for (my $i=0; $i<=2; $i++) {
		$per[$i]=0;
	}
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			if ($j >=$s1[5] && $j < $s1[6]) {
				push @post, $j;
			}	
		}
	}
	if ($s1[2] eq '-') {
		@post=reverse(@post);
	}
	if ($#post >= 5) {
	#splice @post, 0, 3;
	#splice @post, -3;
	for (my $m=0; $m<=$#post; $m++) {
		$len1++;
		my $k=$s1[1].":".$s1[2].":".$post[$m]; 
		if (exists ($read{$k})) {
			$tot+=$read{$k};
			$val[int(($len1-1)/$dist)]+=$read{$k};
			$per[($len1-1)%3]+=$read{$k};
		}
	}
	for (my $m=0; $m<=$#post; $m+=3) {
		my $k1=$s1[1].":".$s1[2].":".$post[$m];
		my $k2=$s1[1].":".$s1[2].":".$post[$m+1];
		my $k3=$s1[1].":".$s1[2].":".$post[$m+2];
		my $n1=0;
		my $n2=0;
		my $n3=0;
		if (exists ($read{$k1})) {
			$n1=$read{$k1};
		}
		if (exists ($read{$k2})) {
			$n2=$read{$k2};
		}
		if (exists ($read{$k3})) {
			$n3=$read{$k3};
		}
		if (($n1+$n2+$n3) > 0) {
			$fra[0]++;
			if ($n1 > $n2 && $n1 > $n3) {
				$fra[1]++;
			}
		}
	}
	if ($tot >= $orfReadCutoff && $len1 >= $orfLengthCutoff) {
		my $ent=0;
		my $ten=0;
		my $a=int(($len1+2)/$dist);
		my $b=$tot;
		my $t1=int(($a+$b-1)/$b);
		my @val2;
		for ($i=0; $i<=$#val; $i++) {
			if ($val[$i] > 0) {
				$val2[int($i/$t1)]+=$val[$i];
			}
		}
		for ($i=0; $i<=$#val2; $i++) {
			if ($val2[$i] > 0) {
				my $p=$val2[$i]/($tot);
				$ent+=$p*log(1/$p);
			}
		}
		my $t2=int($a/$t1);
		my $d1=int($b/$t2);
		my $d2=$b%$t2;
		my @va;
		for (my $i=0; $i<$t2; $i++) {
			$va[$i]=$d1;
		}
		for (my $j=0; $j<$d2; $j++) {
			$va[$j]++;
		}
		for (my $i=0; $i<=$#va; $i++) {
			if ($va[$i] > 0) {
				$p=$va[$i]/($tot);
				$ten+=$p*log(1/$p);
			}
		}
		my $per;
		if ($ten == 0) {
			$per=1;
		} else {
			$per=$ent/$ten;
		}
		if ($tot > 0) {
			my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$s1[5]."\t".$s1[6]."\t".$len1."\t".$tot;
			$out.="\t".sprintf("%.3f", $per[0]/$tot);
			$out.="\t".sprintf("%.3f", $per[1]/$tot);
			$out.="\t".sprintf("%.3f", $per[2]/$tot);
			$out.="\t".sprintf("%.3f", $ent);
			$out.="\t".sprintf("%.3f", $ten);
			$out.="\t".sprintf("%.3f", $per);
			$out.="\t".$fra[0];
			$out.="\t".sprintf("%.3f", $fra[1]/$fra[0]);
			print OUT $out."\n";
		}
	}
	}
}
close AN;
close OUT;

open (RI, ">$outputDir/ribORF.learning.R");
print RI "A <- read.table (\"$outputDir/input.parameters.txt\", sep=\"\\t\", header=T)"."\n";
print RI "B1 <- A[grepl('canonical', A[,1]),]"."\n";
print RI "B2 <- A[grepl('internal', A[,1]),]"."\n";
print RI "C1 <- data.frame(cbind(B1, gr=rep(1, nrow(B1))))"."\n";
print RI "C2 <- data.frame(cbind(B2, gr=rep(0, nrow(B2))))"."\n";
print RI "T1 <- rbind(C1[sample(nrow(C1), min(nrow(C1)/3, 1000)),], C2[sample(nrow(C2), min(nrow(C2)/3, 2000)),])"."\n";
print RI "T2 <- rbind(C1[sample(nrow(C1), min(nrow(C1)/3, 1000)),], C2[sample(nrow(C2), min(nrow(C2)/3, 2000)),])"."\n";
print RI "pred <- glm(gr~f1+PME+f1max, data=T1)"."\n";
print RI "E1 <- predict(pred, T2)"."\n";
print RI "E2 <- cbind(E1, T2\$gr)"."\n";
print RI "a <- seq(-0.5,1.5,0.005)"."\n";
print RI "fpr <- array()"."\n";
print RI "tpr <- array()"."\n";
print RI "out <- matrix(NA, nrow=length(a), ncol=7)"."\n";
print RI "for (i in 1:length(a)) {"."\n";
print RI "tp <- sum(E2[,1] > a[i] & E2[,2]>0)"."\n";
print RI "fp <- sum(E2[,1] > a[i] & E2[,2]==0)"."\n";
print RI "tn <- sum(E2[,1] <= a[i] & E2[,2]==0)"."\n";
print RI "fn <- sum(E2[,1] <= a[i] & E2[,2]>0)"."\n";
print RI "fpr[i] <- fp/(fp+tn)"."\n";
print RI "tpr[i] <- tp/(tp+fn)"."\n";
print RI "out[i,] <- c(a[i], tp, fp, tn, fn, fpr[i], tpr[i])"."\n";
print RI "}"."\n";
print RI "pdf (file=\"$outputDir/plot.ROC.curve.pdf\")"."\n";
print RI "plot(fpr, tpr, col=0, main=\"ROC Curve\", xlab=\"False positive rate\", ylab=\"True positive rate\")"."\n";
print RI "lines(fpr, tpr,col=2, lwd=3)"."\n";
print RI "dev.off()"."\n";
print RI "colnames(out) <- c(\"cutoff\", \"True.pos\", \"False.pos\", \"True.neg\", \"False.neg\", \"False.pos.rate\", \"True.pos.rate\")"."\n";
print RI "write.table (out, \"$outputDir/stat.cutoff.txt\", sep=\"\\t\", quote=F, row.names=F)"."\n";
print RI "pred.pvalue <- sprintf(\"%.4f\", predict(pred, A))"."\n";
print RI "out2 <- data.frame(cbind(A, pred.pvalue))"."\n";
print RI "write.table (out2, \"$outputDir/pred.pvalue.parameters.txt\", sep=\"\\t\", quote=F, row.names=F)"."\n";
close RI;

my $com1="Rscript $outputDir/ribORF.learning.R";
system ($com1);


open (IN, "$outputDir/pred.pvalue.parameters.txt");
my %cluster;
<IN>;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	if ($s1[$#s1] > $predictPvalueCutoff) {
		my @s2=split /\|/, $s1[0];
		my @s3=split /:/, $s2[2];
		my $k1=$s2[0].":".$s3[2];
		$cluster{$k1}.=$s3[1].":".$s2[4].":".$s1[6]."~";
	}
}
close IN;

my %sel;
foreach my $k (sort keys %cluster) {
	my @s1=split /~/, $cluster{$k};
	my $a1;
	my $a2;
	for (my $i=0; $i<=$#s1; $i++) {
		my @s3=split /:/, $s1[$i];
		if ($s3[1] eq 'ATG') {
			$a1.=$s1[$i]."|";
		} else {
			$a2.=$s1[$i]."|";
		} 
	}
	my $a;
	if ($a1 =~ /\|/) {
		$a=$a1;
	} else {
		$a=$a2;
	} 
	my @s2=split /\|/, $a;
	my $b=0;
	my $d;
	my $ind=0;
	for (my $i=0; $i<=$#s2 && $ind==0; $i++) {
		my @s3=split /:/, $s2[$i];
		if ($i < $#s2) {
			my @s4=split /:/, $s2[$i+1];
			$b=$s4[2];
		} else {
			$b=0;
		}
		if ($s3[2] > $b) {
			$ind=1;
			$d=$s3[0].":".$s3[1];
		}
	}
	$sel{$k.":".$d}=1;	
}

open (AN, "$outputDir/pred.pvalue.parameters.txt");
open (OUT, ">$outputDir/repre.valid.pred.pvalue.parameters.txt");
my $line=<AN>;
print OUT $line;
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\|/, $s1[0];
	my @s3=split /:/, $s2[2];
	my $k=$s2[0].":".$s3[2].":".$s3[1].":".$s2[4];
	if (exists ($sel{$k})) {
		print OUT $_."\n";
	}
}
close AN;
close OUT;

open (AN, "$orfFile");
open (OUT, ">$outputDir/repre.valid.ORF.genepred.txt");
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\|/, $s1[0];
	my @s3=split /:/, $s2[2];
	my $k=$s2[0].":".$s3[2].":".$s3[1].":".$s2[4];
	if (exists ($sel{$k})) {
		print OUT $_."\n";
	}
}
close AN;
close OUT;


