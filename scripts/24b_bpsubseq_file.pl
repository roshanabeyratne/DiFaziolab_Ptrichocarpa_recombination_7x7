#! /usr/bin/perl

#  bpsubseq_file.pl: print subsequence from fasta file and input file delineating repeats

# Steve DiFazio, 5-2004
# USAGE: bpsubseq_file.pl -b=<buffer (subtract/add to start/end> <seqfile> <subseqfile>
# Where seqfile is a fasta format file containing the target sequences
# and <subseqfile> is a tab-delimited text file containing <sequence> <start> <end> <complement>
# where <sequence> is the name of the sequence (everything after the > and up to the first whitespace in the fasta header, e.g., Chr01)
# and complement = 1 or -1 (just use 1 for this purpose)


use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use strict;
use Getopt::Long;

my $buff;

GetOptions('b=i' => \$buff);

my $usage= "USAGE: bpsubseq_file.pl -b=<buffer (subtract/add to start/end> <seqfile> <subseqfile>
\tand <subseqfile> contains <sequence> <start> <end> <complement> 
\twhere complement = 1 or -1)\n";
my $seqfile = shift or die $usage;
my $subseqfile = shift or die $usage;
my (%def,$def,%printed,$truncseq,$compl,%compl,$subseq,$seq,%start,%end,$start,$end,$i,$length);
if($subseqfile) {
    open (FILE,$subseqfile) or die "can't open $subseqfile\n";
    while (<FILE>) {
	chomp();
	($seq,$start,$end,$compl,$def) = split(/\t/);
#	print "$seq\t$start\t$end\t$compl\n";
	$start -= $buff;
	if($start <=0) { $start = 1; }
	$end += $buff;
	push(@{$compl{$seq}},$compl);
	push(@{$start{$seq}},$start);
	push(@{$end{$seq}},$end);
	push(@{$def{$seq}},$def);
    }
}

# suck in the fasta sequence file
my $Seq_in = Bio::SeqIO->new (-file => $seqfile, -format =>  'fasta' );

while ( my $query = $Seq_in->next_seq() ) {
    # save name of query sequence
    my $quername = $query->display_id();
    my $length = $query->length();
    if(!($length)) { die "$quername length = $length\n"; }
    foreach $i (0 .. $#{$start{$quername}}) {
	if($start{$quername}[$i] !~ /remove/i) {
	    if($start{$quername}[$i] > $length) { $start{$quername}[$i] = $length; }
	    if($end{$quername}[$i] > $length) { $end{$quername}[$i] = $length; }
	    $truncseq = $query->trunc($start{$quername}[$i],$end{$quername}[$i]);
	    if($compl{$quername}[$i] == -1) { 
		$subseq = $truncseq->revcom();
	    }
	    else { $subseq = $truncseq; }
	    #    my $seqout = Bio::SeqIO->new(-STDOUT, 
	    #				    -format => "fasta");
	    print ">$quername","_",$start{$quername}[$i],"_",$end{$quername}[$i]," ",$def{$quername}[$i],"\n",$subseq->seq,"\n";
	    $printed{$quername} = 1;
	}
    }
}

foreach $seq (sort(keys(%def))) {
    if(!($printed{$seq})) {
	  warn "warning: $seq not found in $seqfile\n"; 
    }
}



