use warnings;
use strict;


my $rep1_entropy = $ARGV[0];
my $rep2_entropy = $ARGV[1];
#   These are output files from make_informationcontent_from_peaks.pl, where column 17 (last column) contains 'entropy' (aka relative information, aka clip_rpr * log2(clip_rpr / input_rpr) for each peak)

my $uid = $ARGV[2];
my $working_directory = $ARGV[3];

my %entropy_hash;
&read_entropy($uid,"01",$rep1_entropy);
&read_entropy($uid,"02",$rep2_entropy);

my $annotated_fi = $working_directory."IDR/".$uid.".01v02.IDR.out.0102merged.bed.annotated_proxdist_miRlncRNA";
my $outfi = $annotated_fi.".entropy";
open(OUTFI,">$outfi");
open(ANN,$annotated_fi) || die "no $annotated_fi\n";
for my $line (<ANN>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    
    my $peak = $tmp[0].":".$tmp[1]."-".$tmp[2].":".$tmp[5];
    my $geometric_mean = log(sqrt( (2 ** $entropy_hash{$peak}{$uid}{"01"}) * (2 ** $entropy_hash{$peak}{$uid}{"02"}) ))/log(2);
    print OUTFI "".$line."\t".sprintf("%.10f",$geometric_mean)."\n";
    
    
}
close(ANN);
close(OUTFI);



sub read_entropy {
    my $uid = shift;
    my $rep = shift;

    my $file = shift;
    open(F,$file) || die "no $file\n";
    for my $line (<F>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $stop = $tmp[2];

	my ($chr2,$pos,$str,$del) = split(/\:/,$tmp[3]);
	
	my $peak = $chr.":".$start."-".$stop.":".$str;

	my $entropy = $tmp[17];
	$entropy_hash{$peak}{$uid}{$rep} = $entropy;
    }
    close(F);


}
