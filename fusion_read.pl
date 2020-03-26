# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 10000),
	's' => \(my $sort = ''),
	'q=i' => \(my $minimumMappingQuality = 0),
	'd=i' => \(my $maximumInnerDistance = 1000),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl fusion_read.pl [options] bwa.sam [...] > fusion_read.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -s       input is neither raw BWA output nor sorted by read name
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -d INT   maximum inner-distance [$maximumInnerDistance]

EOF
}
{
	my $parentPid = $$;
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryDirectory/fork.$hostname.$parentPid.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryDirectory/fork.$hostname.$parentPid.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryDirectory/fork.$hostname.$parentPid.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my (@samFileList) = @ARGV;
foreach my $samFile (@samFileList) {
	open(my $reader, $sort eq '' ? ($samFile =~ /\.gz$/ ? "gzip -dc $samFile | grep -v '^\@' |" : "grep -v '^\@' $samFile |") : ($samFile =~ /\.gz$/ ? "gzip -dc $samFile | grep -v '^\@' | LC_ALL=C sort -t '\t' -k1,1 |" : "grep -v '^\@' $samFile | LC_ALL=C sort -t '\t' -k1,1 |"));
	my ($readName, @tokenHashList) = ('');
	my @tokenHashListList = ();
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		if($tokenHash{'qname'} ne $readName) {
			push(@tokenHashListList, [@tokenHashList]) if(@tokenHashList);
			if(scalar(@tokenHashListList) >= $numberPerThread) {
				if($threads == 1) {
					printFusion(@tokenHashListList);
				} else {
					forkPrintSubroutine(\&printFusion, @tokenHashListList);
				}
				@tokenHashListList = ();
			}
			($readName, @tokenHashList) = ($tokenHash{'qname'});
		}
		push(@tokenHashList, \%tokenHash);
	}
	push(@tokenHashListList, [@tokenHashList]) if(@tokenHashList);
	if(@tokenHashListList) {
		if($threads == 1) {
			printFusion(@tokenHashListList);
		} else {
			forkPrintSubroutine(\&printFusion, @tokenHashListList);
		}
		@tokenHashListList = ();
	}
	close($reader);
	forkPrintWait();
}

sub printFusion {
	my (@tokenHashListList) = @_;
	foreach my $tokenHashList (@tokenHashListList) {
		my %readNameHash = ();
		my %numberTokenHashListHash = ();
		foreach(@$tokenHashList) {
			my %tokenHash = %$_;
			$readNameHash{$tokenHash{'qname'}} = 1;
			my $number = ($tokenHash{'flag'} & 192) / 64;
			push(@{$numberTokenHashListHash{$number}}, \%tokenHash);
		}
		my %numberSequenceHash = ();
		my %numberUnmappedHash = ();
		my %numberChromosomeStrandPositionListHash = ();
		my %numberChromosomeStrandPositionLengthListHash1 = ();
		my %numberChromosomeStrandPositionLengthListHash2 = ();
		foreach my $number (keys %numberTokenHashListHash) {
			my @tokenHashList = @{$numberTokenHashListHash{$number}};
			my %sequenceHash = ();
			foreach(@tokenHashList) {
				my %tokenHash = %$_;
				my $sequence = $tokenHash{'seq'};
				$sequence = getReverseComplementarySequence($sequence) if($tokenHash{'flag'} & 16);
				$sequenceHash{$sequence} = 1;
			}
			next if(scalar(($numberSequenceHash{$number}) = keys %sequenceHash) > 1);
			my @chromosomePositionCigarStrandList = ();
			foreach(@tokenHashList) {
				my %tokenHash = %$_;
				if($tokenHash{'flag'} & 4) {
					$numberUnmappedHash{$number} = 1;
				} elsif($tokenHash{'mapq'} >= $minimumMappingQuality) {
					push(@chromosomePositionCigarStrandList, [@tokenHash{'rname', 'pos', 'cigar'}, ($tokenHash{'flag'} & 16) ? '-' : '+']);
					while(defined($tokenHash{'XA:Z'}) && $tokenHash{'XA:Z'} =~ /([^,;]+),([-+][0-9]+),([^,;]+),([0-9]+);/g) {
						push(@chromosomePositionCigarStrandList, [$1, abs($2), $3, ($2 < 0) ? '-' : '+']);
					}
				}
			}
			foreach(@chromosomePositionCigarStrandList) {
				my ($chromosome, $position, $cigar, $strand) = @$_;
				my @positionList = getPositionList($position, $cigar);
				my ($start, $end) = ($position, grep {$_ ne ''} reverse(@positionList));
				if($cigar =~ /^([0-9]+)M/) {
					my $length = $1;
					push(@{$numberChromosomeStrandPositionListHash{$number}}, [$chromosome, $strand, $start]) if($strand eq '-');
					if($cigar =~ /([0-9]+)[IS]$/) {
						push(@{$numberChromosomeStrandPositionLengthListHash1{$number}}, [$chromosome, $strand, $start, $length]) if($strand eq '+');
						push(@{$numberChromosomeStrandPositionLengthListHash2{$number}}, [$chromosome, $strand, $start, $length]) if($strand eq '-');
					}
				}
				if($cigar =~ /([0-9]+)M$/) {
					my $length = $1;
					push(@{$numberChromosomeStrandPositionListHash{$number}}, [$chromosome, $strand, $end]) if($strand eq '+');
					if($cigar =~ /^([0-9]+)[IS]/) {
						push(@{$numberChromosomeStrandPositionLengthListHash1{$number}}, [$chromosome, $strand, $end, $length]) if($strand eq '-');
						push(@{$numberChromosomeStrandPositionLengthListHash2{$number}}, [$chromosome, $strand, $end, $length]) if($strand eq '+');
					}
				}
			}
		}
		my %fusionNumberHash = ();
		foreach my $number (keys %numberTokenHashListHash) {
			my $length = length(my $sequence = $numberSequenceHash{$number});
			foreach my $chromosomeStrandPositionLength1 (@{$numberChromosomeStrandPositionLengthListHash1{$number}}) {
				foreach my $chromosomeStrandPositionLength2 (@{$numberChromosomeStrandPositionLengthListHash2{$number}}) {
					my ($chromosome1, $strand1, $position1, $length1) = @$chromosomeStrandPositionLength1;
					my ($chromosome2, $strand2, $position2, $length2) = @$chromosomeStrandPositionLength2;
					my @mateChromosomeStrandPositionList = ();
					if(($number == 1 || $number == 2) && defined(my $chromosomeStrandPositionList = $numberChromosomeStrandPositionListHash{3 - $number})) {
						push(@mateChromosomeStrandPositionList, grep {$_->[0] eq $chromosome2 && $_->[1] ne $strand2 && (($strand2 eq '+' && $position2 - $length2 + 1 <= $_->[2] && ($_->[2] - 1) - $position2 <= $maximumInnerDistance) || ($strand2 eq '-' && $_->[2] <= $position2 + $length2 - 1 && ($position2 - 1) - $_->[2] <= $maximumInnerDistance))} @$chromosomeStrandPositionList);
					}
					my $apart = '';
					if($length1 + $length2 >= $length) {
						($length1, $length2) = ($length - $length2, $length - $length1);
					} else {
						$apart = 1;
					}
					my $intermediateSequence = substr($sequence, $length1, -$length2);
					$position1 = $position1 + $length1 - 1 if($strand1 eq '+');
					$position2 = $position2 + $length2 - 1 if($strand2 eq '-');
					$position1 = $position1 - $length1 + 1 if($strand1 eq '-');
					$position2 = $position2 - $length2 + 1 if($strand2 eq '+');
					if(($chromosome1 lt $chromosome2) || ($chromosome1 eq $chromosome2 && $position1 <= $position2)) {
						my $fusion = join("\t", $chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2, $intermediateSequence, $apart, ($number == 2 ? '-' : '+'));
						$fusionNumberHash{$fusion}->{$number} = 1;
						$fusionNumberHash{$fusion}->{3 - $number} = 1 if(@mateChromosomeStrandPositionList);
					} else {
						my $fusion = join("\t", $chromosome2, $position2, ($strand2 eq '+' ? '-' : $strand2 eq '-' ? '+' : ''), $chromosome1, $position1, ($strand1 eq '+' ? '-' : $strand1 eq '-' ? '+' : ''), getReverseComplementarySequence($intermediateSequence), $apart, ($number == 2 ? '+' : '-'));
						$fusionNumberHash{$fusion}->{$number} = 1;
						$fusionNumberHash{$fusion}->{3 - $number} = 1 if(@mateChromosomeStrandPositionList);
					}
				}
				if(($number == 1 || $number == 2) && scalar(@{$numberChromosomeStrandPositionLengthListHash2{$number}}) == 0) {
					my ($chromosome1, $strand1, $position1, $length1) = @$chromosomeStrandPositionLength1;
					my ($chromosome2, $strand2, $position2, $length2) = ('', '', '', '');
					my $apart = '';
					my $intermediateSequence = join('^', substr($sequence, $length1), getReverseComplementarySequence($numberSequenceHash{3 - $number}));
					$position1 = $position1 + $length1 - 1 if($strand1 eq '+');
					$position1 = $position1 - $length1 + 1 if($strand1 eq '-');
					{
						my $fusion = join("\t", $chromosome1, $position1, $strand1, $chromosome2, $position2, $strand2, $intermediateSequence, $apart, ($number == 2 ? '-' : '+'));
						$fusionNumberHash{$fusion}->{$number} = 1;
						$fusionNumberHash{$fusion}->{3 - $number} = 1 if($numberUnmappedHash{3 - $number});
					}
				}
			}
			if(($number == 1 || $number == 2) && scalar(@{$numberChromosomeStrandPositionLengthListHash1{$number}}) == 0) {
				foreach my $chromosomeStrandPositionLength2 (@{$numberChromosomeStrandPositionLengthListHash2{$number}}) {
					my ($chromosome1, $strand1, $position1, $length1) = ('', '', '', '');
					my ($chromosome2, $strand2, $position2, $length2) = @$chromosomeStrandPositionLength2;
					my @mateChromosomeStrandPositionList = ();
					if(($number == 1 || $number == 2) && defined(my $chromosomeStrandPositionList = $numberChromosomeStrandPositionListHash{3 - $number})) {
						push(@mateChromosomeStrandPositionList, grep {$_->[0] eq $chromosome2 && $_->[1] ne $strand2 && (($strand2 eq '+' && $position2 - $length2 + 1 <= $_->[2] && ($_->[2] - 1) - $position2 <= $maximumInnerDistance) || ($strand2 eq '-' && $_->[2] <= $position2 + $length2 - 1 && ($position2 - 1) - $_->[2] <= $maximumInnerDistance))} @$chromosomeStrandPositionList);
					}
					my $apart = '';
					my $intermediateSequence = substr($sequence, 0, -$length2);
					$position2 = $position2 + $length2 - 1 if($strand2 eq '-');
					$position2 = $position2 - $length2 + 1 if($strand2 eq '+');
					{
						my $fusion = join("\t", $chromosome2, $position2, ($strand2 eq '+' ? '-' : $strand2 eq '-' ? '+' : ''), $chromosome1, $position1, ($strand1 eq '+' ? '-' : $strand1 eq '-' ? '+' : ''), getReverseComplementarySequence($intermediateSequence), $apart, ($number == 2 ? '+' : '-'));
						$fusionNumberHash{$fusion}->{$number} = 1;
						$fusionNumberHash{$fusion}->{3 - $number} = 1 if(@mateChromosomeStrandPositionList);
					}
				}
			}
		}
		foreach my $fusion (sort keys %fusionNumberHash) {
			if($fusionNumberHash{$fusion}->{0}) {
				forkPrint(join("\t", $fusion, join(',', sort keys %readNameHash)), "\n");
			} elsif($fusionNumberHash{$fusion}->{1} && $fusionNumberHash{$fusion}->{2}) {
				forkPrint(join("\t", $fusion, join(',', sort keys %readNameHash)), "\n");
			}
		}
	}
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
