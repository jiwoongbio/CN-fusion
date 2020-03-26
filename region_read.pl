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
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl region_read.pl [options] bwa.sam [...] > region_read.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -s       input is neither raw BWA output nor sorted by read name
         -q INT   minimum mapping quality [$minimumMappingQuality]

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
open(my $writer, "| LC_ALL=C sort -t '\t' -k1,1 -k2,2n -k3,3n");
forkPrintParentWriter($writer);
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
					printRegion(@tokenHashListList);
				} else {
					forkPrintSubroutine(\&printRegion, @tokenHashListList);
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
			printRegion(@tokenHashListList);
		} else {
			forkPrintSubroutine(\&printRegion, @tokenHashListList);
		}
		@tokenHashListList = ();
	}
	forkPrintWait();
	close($reader);
}
forkPrintParentWriter();
close($writer);

sub printRegion {
	my (@tokenHashListList) = @_;
	foreach my $tokenHashList (@tokenHashListList) {
		foreach(@$tokenHashList) {
			my %tokenHash = %$_;
			if($tokenHash{'flag'} & 4) {
			} elsif($tokenHash{'mapq'} >= $minimumMappingQuality) {
				{
					my ($chromosome, $position, $cigar, $strand) = (@tokenHash{'rname', 'pos', 'cigar'}, ($tokenHash{'flag'} & 16) ? '-' : '+');
					$strand = $strand eq '+' ? '-' : '+' if($tokenHash{'flag'} & 128);
					my @positionList = grep {$_ ne ''} getPositionList($position, $cigar);
					my ($start, $end) = @positionList[0, -1];
					forkPrint(join("\t", $chromosome, $start, $end, $strand, $tokenHash{'qname'}), "\n");
				}
				while(defined($tokenHash{'XA:Z'}) && $tokenHash{'XA:Z'} =~ /([^,;]+),([-+][0-9]+),([^,;]+),([0-9]+);/g) {
					my ($chromosome, $position, $cigar, $strand) = ($1, abs($2), $3, ($2 < 0) ? '-' : '+');
					$strand = $strand eq '+' ? '-' : '+' if($tokenHash{'flag'} & 128);
					my @positionList = grep {$_ ne ''} getPositionList($position, $cigar);
					my ($start, $end) = @positionList[0, -1];
					forkPrint(join("\t", $chromosome, $start, $end, $strand, $tokenHash{'qname'}), "\n");
				}
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
