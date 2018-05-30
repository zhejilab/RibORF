if ($#ARGV < 2) {
	print "## remove adapter sequence"."\n";
	print "usage: perl removeAdaptor.pl -f fastqFile -a adapterSequence -o outputFile [-l readLengthCutoff]"."\n";
	print "fastqFile: raw fastq read sequences;"."\n";
	print "adapterSequence: 5' end sequence of adapter, 10nt is recommended;"."\n";
	print "outputFile: output file;"."\n";
	print "readLengthCutoff [optional]: minimal read length after trimming adapters, default is 15nt."."\n";
	exit;
}

use Getopt::Std;

### get arguments ###
my %args; 
getopt("faol",\%args);
my $file=$args{f};
my $adapter=$args{a};
my $outputFile=$args{o};
my $cutoff=$args{l};

if (! $cutoff) {
	$cutoff=15;
}

open (IN, "$file");
open (OUT, ">$outputFile");

my $num=0;
my $out;
my $m;
while (<IN>) {
	chomp;
	$num++;
	if ($num%4==1) {
		$out=$_."\n";
	} elsif ($num%4==2) {
		if ($_ =~ /$adapter/) {
			my @s=split /$adapter/, $_;
			$m=$s[0];
		} else {
			$m=$_;
		}
		$out.=$m."\n";
	} elsif ($num%4==3) {
		$out.=$_."\n";
	} elsif ($num%4==0) {
		if (length($m) >= $cutoff) {
			print OUT $out.substr($_, 0, length($m))."\n";
		}
	} 
}
close IN;
close OUT;
