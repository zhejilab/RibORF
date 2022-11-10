my $file = "./sta.read.dist.25,26,27,28,29,30,31,32,33,34,.txt";

open (IN, "$file");
open (OUT, ">offset.corretion.parameters.txt");
<IN>;
while (<IN>) {
	chomp;
	my @s=split /\t/, $_;
	my @va;
	push @va, $s[2];
	push @va, $s[3];
	push @va, $s[4];
	@va = sort {$a <=> $b} @va;
	if ($va[2] > 0.6 && ($va[2] - $va[1]) > 0.2) {
		my $ind;
		if ($s[4] > 0.6) {
			$ind=15;
			if ($s[0] == 25) {
				$ind=12;
			}
		} elsif  ($s[3] > 0.6) {
			$ind=14;
			if ($s[0] >= 30) {
				$ind=17;
			}
		} else {
			$ind=16;
			if ($s[0] < 28) {
				$ind=13;
			}
		}
		print OUT $f1[6]."\t".$f1[$#f1-1]."\t".$_."\t".$ind."\n";
	}
}
close IN;
close OUT;
