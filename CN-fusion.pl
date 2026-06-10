# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Bio::DB::Fasta;
use IPC::Open2;
use List::Util qw(sum);
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my @codonList = ();
my @regionReadFileList = ();
my @fastaFileList = ();
if(-e (my $file = "$dataPath/rna_genomic.fasta")) {
	push(@fastaFileList, $file);
}
my @proteinFastaFileList = ();
if(-e (my $file = "$dataPath/protein.fasta")) {
	push(@proteinFastaFileList, $file);
}
my @proteinCodingFileList = ();
if(-e (my $file = "$dataPath/protein_coding.txt")) {
	push(@proteinCodingFileList, $file);
}
my @proteinFeatureFileList = ();
if(-e (my $file = "$dataPath/protein_feature.txt")) {
	push(@proteinFeatureFileList, $file);
}
my @spliceJunctionFileList = ();
if(-e (my $file = "$dataPath/splice_junction.txt")) {
	push(@spliceJunctionFileList, $file);
}
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 10000),
	's' => \(my $sort = ''),
	'S=s' => \(my $stranded = ''),
	'C=s' => \@codonList,
	'c=i' => \(my $minimumCount = 3),
	'd=i' => \(my $minimumBreadthDepth = 1),
	'r=s' => \@regionReadFileList,
	'fastaFile=s' => \@fastaFileList,
	'proteinCodingFile=s' => \@proteinCodingFileList,
	'proteinFastaFile=s' => \@proteinFastaFileList,
	'proteinFeatureFile=s' => \@proteinFeatureFileList,
	'spliceJunctionFile=s' => \@spliceJunctionFileList,
	'keepSelfFusion' => \(my $keepSelfFusion = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl CN-fusion.pl [options] fusion_read.txt [...] > CN-fusion.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -s       input is neither raw fusion_read output nor sorted by read name
         -S STR   stranded, "f" or "r"
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 1 (standard)]
         -c INT   minimum count [$minimumCount]
         -d INT   minimum depth for breadth calculation [$minimumBreadthDepth]
         -r FILE  region read file
         --numberPerThread INT
                  number of fusion candidates per thread [$numberPerThread]
         --fastaFile FILE
                  transcript/genomic FASTA file
         --proteinCodingFile FILE
                  protein coding annotation file
         --proteinFastaFile FILE
                  protein FASTA file
         --proteinFeatureFile FILE
                  protein feature annotation file
         --spliceJunctionFile FILE
                  splice junction annotation file
         --keepSelfFusion
                  keep same-transcript/same-gene inframe fusion candidates otherwise excluded

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
			open($writer, "> $temporaryDirectory/fork.$hostname.$parentPid.$$") or die "Can't open '$temporaryDirectory/fork.$hostname.$parentPid.$$': $!";
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
				die "Child process $pid failed with status $?" if $? != 0;
				open(my $reader, "$temporaryDirectory/fork.$hostname.$parentPid.$pid") or die "Can't open '$temporaryDirectory/fork.$hostname.$parentPid.$pid': $!";
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
{
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. int(length($sequence) / 3) - 1);
	}
}
my %seqHash = ();
foreach my $file (@fastaFileList) {
	my $db = Bio::DB::Fasta->new($file);
	foreach my $sequenceId ($db->get_all_primary_ids) {
		$seqHash{$sequenceId} = $db->get_Seq_by_id($sequenceId);
	}
}
my %proteinSequenceHash = ();
foreach my $file (@proteinFastaFileList) {
	my $proteinId = '';
	open(my $reader, ($file =~ /\.gz$/ ? "gzip -dc $file |" : $file)) or die "Can't open '$file': $!";
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($proteinId = $1));
		$proteinSequenceHash{$proteinId} .= $line;
	}
	close($reader);
	s/U/*/g foreach(values %proteinSequenceHash);
	s/\*?$/*/ foreach(values %proteinSequenceHash);
}
my %transcriptProteinCodingListHash = ();
foreach my $file (@proteinCodingFileList) {
	open(my $reader, ($file =~ /\.gz$/ ? "gzip -dc $file |" : $file)) or die "Can't open '$file': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my ($gene, $transcriptId, $proteinId, $codingRegions, $frame) = split(/\t/, $line, -1);
		push(@{$transcriptProteinCodingListHash{$transcriptId}}, [$gene, $proteinId, $codingRegions, $frame]);
	}
	close($reader);
}
my %proteinFeatureListHash = ();
foreach my $file (@proteinFeatureFileList) {
	open(my $reader, ($file =~ /\.gz$/ ? "gzip -dc $file |" : $file)) or die "Can't open '$file': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my ($proteinId, $start, $end, $proteinFeature) = split(/\t/, $line, -1);
		push(@{$proteinFeatureListHash{$proteinId}}, [$start, $end, $proteinFeature]);
	}
	close($reader);
}
my %transcriptPositionSpliceJunctionHash = ();
foreach my $file (@spliceJunctionFileList) {
	open(my $reader, ($file =~ /\.gz$/ ? "gzip -dc $file |" : $file)) or die "Can't open '$file': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my ($transcriptId, $position, $spliceJunction1, $spliceJunction2) = split(/\t/, $line, -1);
		next if(defined($transcriptPositionSpliceJunctionHash{$transcriptId}->{$position}));
		$transcriptPositionSpliceJunctionHash{$transcriptId}->{$position} = [$spliceJunction1, $spliceJunction2];
	}
	close($reader);
}
my @fusionColumnList = ('sequenceId1', 'breakpointPosition1', 'breakpointStrand1', 'sequenceId2', 'breakpointPosition2', 'breakpointStrand2', 'insertedSequence', 'isApart');
my @fusionProteinCodingColumnList1 = (@fusionColumnList, 'consequenceType', 'spliceJunction1', 'spliceJunction2');
{
	push(@fusionProteinCodingColumnList1, 'gene1', 'proteinId1', 'codingRegions1', 'frame1', 'codingIndex1', 'codingStart1', 'codingEnd1');
	push(@fusionProteinCodingColumnList1, 'gene2', 'proteinId2', 'codingRegions2', 'frame2', 'codingIndex2', 'codingStart2', 'codingEnd2');
}
my @fusionProteinCodingColumnList2 = (@fusionColumnList, 'fusionReadCount', 'consequenceType', 'spliceJunction1', 'spliceJunction2', 'codingStart1', 'codingEnd1', 'codingStart2', 'codingEnd2');
if(@regionReadFileList) {
	push(@fusionProteinCodingColumnList2, 'gene1', 'proteinId1', 'proteinMatchLengthFromN1', 'proteinMatchLengthFromC1', 'proteinMatchFraction1', 'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1', 'fusionReadFraction1', 'meanDepth1', 'breadth1');
	push(@fusionProteinCodingColumnList2, 'gene2', 'proteinId2', 'proteinMatchLengthFromN2', 'proteinMatchLengthFromC2', 'proteinMatchFraction2', 'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2', 'fusionReadFraction2', 'meanDepth2', 'breadth2');
} else {
	push(@fusionProteinCodingColumnList2, 'gene1', 'proteinId1', 'proteinMatchLengthFromN1', 'proteinMatchLengthFromC1', 'proteinMatchFraction1', 'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1');
	push(@fusionProteinCodingColumnList2, 'gene2', 'proteinId2', 'proteinMatchLengthFromN2', 'proteinMatchLengthFromC2', 'proteinMatchFraction2', 'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2');
}
push(@fusionProteinCodingColumnList2, 'fusionProteinSequence');
my (@fusionReadFileList) = @ARGV;
if(@fusionReadFileList) {
	my $pid = open2(my $reader, my $writer, 'LC_ALL=C sort');
	forkPrintParentWriter($writer);
	{
		my $pid = open2(my $reader, my $writer, 'LC_ALL=C sort | uniq -c');
		foreach my $fusionReadFile (@fusionReadFileList) {
			open(my $reader, $sort eq '' ? ($fusionReadFile =~ /\.gz$/ ? "gzip -dc $fusionReadFile |" : $fusionReadFile) : ($fusionReadFile =~ /\.gz$/ ? "gzip -dc $fusionReadFile | LC_ALL=C sort -t '\t' -k10,10 |" : "LC_ALL=C sort -t '\t' -k10,10 $fusionReadFile |"));
			my ($readName, @tokenList) = ('');
			while(my $line = <$reader>) {
				chomp($line);
				my %tokenHash = ();
				@tokenHash{@fusionColumnList, 'strand', 'readName'} = split(/\t/, $line, -1);
				if($tokenHash{'readName'} ne $readName) {
					print $writer join("\t", @tokenList), "\n" if(@tokenList);
					($readName, @tokenList) = ($tokenHash{'readName'});
				}
				if($tokenHash{'breakpointStrand1'} eq '+') {
					if($stranded eq '' || ($stranded eq 'f' && $tokenHash{'strand'} eq '+') || ($stranded eq 'r' && $tokenHash{'strand'} eq '-')) {
						push(@tokenList, @tokenHash{@fusionColumnList});
					}
				}
				if($tokenHash{'breakpointStrand2'} eq '-') {
					@tokenHash{'sequenceId1', 'breakpointPosition1', 'breakpointStrand1', 'sequenceId2', 'breakpointPosition2', 'breakpointStrand2', 'insertedSequence', 'strand'} = (@tokenHash{'sequenceId2', 'breakpointPosition2'}, '+', @tokenHash{'sequenceId1', 'breakpointPosition1'}, ($tokenHash{'breakpointStrand1'} eq '+' ? '-' : $tokenHash{'breakpointStrand1'} eq '-' ? '+' : ''), getReverseComplementarySequence($tokenHash{'insertedSequence'}), ($tokenHash{'strand'} eq '+' ? '-' : $tokenHash{'strand'} eq '-' ? '+' : ''));
					if($stranded eq '' || ($stranded eq 'f' && $tokenHash{'strand'} eq '+') || ($stranded eq 'r' && $tokenHash{'strand'} eq '-')) {
						push(@tokenList, @tokenHash{@fusionColumnList});
					}
				}
			}
			print $writer join("\t", @tokenList), "\n" if(@tokenList);
			close($reader);
		}
		close($writer);
		my @tokenHashListList = ();
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ s/^ *([0-9]+) //) {
				my $count = $1;
				my @tokenHashList = ();
				my @tokenList = split(/\t/, $line, -1);
				while(@tokenList) {
					my %tokenHash = ();
					@tokenHash{@fusionColumnList, 'fusionReadCount'} = (splice(@tokenList, 0, scalar(@fusionColumnList)), $count);
					push(@tokenHashList, \%tokenHash);
				}
				push(@tokenHashListList, \@tokenHashList);
			}
			if(scalar(@tokenHashListList) >= $numberPerThread) {
				if($threads == 1) {
					printFusionProtein1(@tokenHashListList);
				} else {
					forkPrintSubroutine(\&printFusionProtein1, @tokenHashListList);
				}
				@tokenHashListList = ();
			}
		}
		if(@tokenHashListList) {
			if($threads == 1) {
				printFusionProtein1(@tokenHashListList);
			} else {
				forkPrintSubroutine(\&printFusionProtein1, @tokenHashListList);
			}
			@tokenHashListList = ();
		}
		forkPrintWait();
		close($reader);
		waitpid($pid, 0);
	}
	forkPrintParentWriter();
	close($writer);
	{
		open(my $writer, "| LC_ALL=C sort -t '\t' -k9,9nr");
		forkPrintParentWriter($writer);
		my ($fusionGene, $count) = ('', 0);
		my @tokenHashList = ();
		print join("\t", @fusionProteinCodingColumnList2), "\n";
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@fusionProteinCodingColumnList1, 'fusionReadCount'} = split(/\t/, $line, -1);
			if((my $currentFusionGene = join("\t", @tokenHash{@fusionProteinCodingColumnList1})) ne $fusionGene) {
				if($fusionGene ne '' && $count >= $minimumCount) {
					my %tokenHash = ();
					@tokenHash{@fusionProteinCodingColumnList1, 'fusionReadCount'} = (split(/\t/, $fusionGene, -1), $count);
					push(@tokenHashList, \%tokenHash);
				}
				if(scalar(@tokenHashList) >= $numberPerThread) {
					if($threads == 1) {
						printFusionProtein2(@tokenHashList);
					} else {
						forkPrintSubroutine(\&printFusionProtein2, @tokenHashList);
					}
					@tokenHashList = ();
				}
				($fusionGene, $count) = ($currentFusionGene, 0);
			}
			$count += $tokenHash{'fusionReadCount'};
		}
		if($fusionGene ne '' && $count >= $minimumCount) {
			my %tokenHash = ();
			@tokenHash{@fusionProteinCodingColumnList1, 'fusionReadCount'} = (split(/\t/, $fusionGene, -1), $count);
			push(@tokenHashList, \%tokenHash);
		}
		if(@tokenHashList) {
			if($threads == 1) {
				printFusionProtein2(@tokenHashList);
			} else {
				forkPrintSubroutine(\&printFusionProtein2, @tokenHashList);
			}
			@tokenHashList = ();
		}
		forkPrintWait();
		forkPrintParentWriter();
		close($writer);
	}
	close($reader);
	waitpid($pid, 0);
}

{
	my %sequenceIdProteinCodingTokenHashListHash = ();
	sub getSequenceIdProteinCodingTokenHashListHash1 {
		my ($sequenceId) = @_;
		my @proteinCodingTokenHashList = ();
		if(defined(my $proteinCodingTokenHashList = $sequenceIdProteinCodingTokenHashListHash{$sequenceId})) {
			@proteinCodingTokenHashList = @$proteinCodingTokenHashList;
		} else {
			foreach my $proteinCoding (@{$transcriptProteinCodingListHash{$sequenceId}}) {
				my $proteinCodingTokenHash = {};
				@$proteinCodingTokenHash{'gene1', 'proteinId1', 'codingRegions1', 'frame1'} = @$proteinCoding;
				my @codingPositionList = getCodingPositionList(@$proteinCodingTokenHash{'codingRegions1', 'frame1'});
				$proteinCodingTokenHash->{'codingStart1'} = $codingPositionList[0];
				$proteinCodingTokenHash->{'codingEnd1'} = $codingPositionList[-1];
				$proteinCodingTokenHash->{'codingPositionIndexHash1'}->{$codingPositionList[$_]} = $_ foreach(0 .. $#codingPositionList);
				push(@proteinCodingTokenHashList, $proteinCodingTokenHash);
			}
			@proteinCodingTokenHashList = sort {$a->{'codingStart1'} <=> $b->{'codingStart1'} || $a->{'codingEnd1'} <=> $b->{'codingEnd1'}} @proteinCodingTokenHashList;
			$sequenceIdProteinCodingTokenHashListHash{$sequenceId} = \@proteinCodingTokenHashList;
		}
		return @proteinCodingTokenHashList;
	}
}
{
	my %sequenceIdProteinCodingTokenHashListHash = ();
	sub getSequenceIdProteinCodingTokenHashListHash2 {
		my ($sequenceId) = @_;
		my @proteinCodingTokenHashList = ();
		if(defined(my $proteinCodingTokenHashList = $sequenceIdProteinCodingTokenHashListHash{$sequenceId})) {
			@proteinCodingTokenHashList = @$proteinCodingTokenHashList;
		} else {
			foreach my $proteinCoding (@{$transcriptProteinCodingListHash{$sequenceId}}) {
				my $proteinCodingTokenHash = {};
				@$proteinCodingTokenHash{'gene2', 'proteinId2', 'codingRegions2', 'frame2'} = @$proteinCoding;
				my @codingPositionList = getCodingPositionList(@$proteinCodingTokenHash{'codingRegions2', 'frame2'});
				$proteinCodingTokenHash->{'codingStart2'} = $codingPositionList[0];
				$proteinCodingTokenHash->{'codingEnd2'} = $codingPositionList[-1];
				$proteinCodingTokenHash->{'codingPositionIndexHash2'}->{$codingPositionList[$_]} = $_ foreach(0 .. $#codingPositionList);
				push(@proteinCodingTokenHashList, $proteinCodingTokenHash);
			}
			@proteinCodingTokenHashList = sort {$a->{'codingStart2'} <=> $b->{'codingStart2'} || $a->{'codingEnd2'} <=> $b->{'codingEnd2'}} @proteinCodingTokenHashList;
			$sequenceIdProteinCodingTokenHashListHash{$sequenceId} = \@proteinCodingTokenHashList;
		}
		return @proteinCodingTokenHashList;
	}
}
sub printFusionProtein1 {
	my (@tokenHashListList) = @_;
	foreach(@tokenHashListList) {
		my @tokenHashList = @$_;
		my @preferredConsequenceTokenHashList = ();
		my @fallbackConsequenceTokenHashList = ();
		foreach my $tokenHash (@tokenHashList) {
			next if($tokenHash->{'sequenceId2'} eq '');
			if((my $insertedSequenceLength = length($tokenHash->{'insertedSequence'})) > 0 && $tokenHash->{'isApart'} eq '' && $tokenHash->{'breakpointStrand1'} eq '+' && $tokenHash->{'breakpointStrand2'} eq '+') {
				if(defined(my $positionSpliceJunctionHash1 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'sequenceId1'}}) && defined(my $positionSpliceJunctionHash2 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'sequenceId2'}})) {
					foreach my $index (0 .. $insertedSequenceLength) {
						my $position1 = $tokenHash->{'breakpointPosition1'} + $index;
						my $position2 = $tokenHash->{'breakpointPosition2'} - ($insertedSequenceLength - $index);
						if(defined($positionSpliceJunctionHash1->{$position1}) && defined($positionSpliceJunctionHash2->{$position2 - 1}) && defined(my $spliceJunction1 = $positionSpliceJunctionHash1->{$position1}->[0]) && defined(my $spliceJunction2 = $positionSpliceJunctionHash2->{$position2 - 1}->[1])) {
							$tokenHash->{'breakpointPosition1'} = $position1;
							$tokenHash->{'breakpointPosition2'} = $position2;
							$tokenHash->{'spliceJunction1'} = $spliceJunction1;
							$tokenHash->{'spliceJunction2'} = $spliceJunction2;
							$tokenHash->{'insertedSequence'} = '';
							last;
						}
					}
				}
			}
			if((my $insertedSequenceLength = length($tokenHash->{'insertedSequence'})) > 0 && $tokenHash->{'isApart'} eq '' && $tokenHash->{'breakpointStrand1'} eq '+') {
				if(defined(my $positionSpliceJunctionHash1 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'sequenceId1'}})) {
					foreach my $index (0 .. $insertedSequenceLength) {
						my $position1 = $tokenHash->{'breakpointPosition1'} + $index;
						if(defined($positionSpliceJunctionHash1->{$position1}) && defined(my $spliceJunction1 = $positionSpliceJunctionHash1->{$position1}->[0])) {
							$tokenHash->{'breakpointPosition1'} = $position1;
							$tokenHash->{'spliceJunction1'} = $spliceJunction1;
							$tokenHash->{'insertedSequence'} = substr($tokenHash->{'insertedSequence'}, $index);
							last;
						}
					}
				}
			}
			if($tokenHash->{'breakpointStrand1'} eq '+' && defined($transcriptProteinCodingListHash{$tokenHash->{'sequenceId1'}})) {
				my @proteinCodingTokenHashList1 = getSequenceIdProteinCodingTokenHashListHash1($tokenHash->{'sequenceId1'});
				foreach(@proteinCodingTokenHashList1) {
					my %proteinCodingTokenHash1 = %$_;
					if(defined($proteinCodingTokenHash1{'codingIndex1'} = $proteinCodingTokenHash1{'codingPositionIndexHash1'}->{$tokenHash->{'breakpointPosition1'}})) {
						if($tokenHash->{'breakpointStrand2'} eq '+' && defined($transcriptProteinCodingListHash{$tokenHash->{'sequenceId2'}})) {
							my @inframeProteinCodingTokenHashList2 = ();
							my @proteinCodingTokenHashList2 = getSequenceIdProteinCodingTokenHashListHash2($tokenHash->{'sequenceId2'});
							foreach(@proteinCodingTokenHashList2) {
								my %proteinCodingTokenHash2 = %$_;
								if(defined($proteinCodingTokenHash2{'codingIndex2'} = $proteinCodingTokenHash2{'codingPositionIndexHash2'}->{$tokenHash->{'breakpointPosition2'}})) {
									if(($proteinCodingTokenHash1{'codingIndex1'} + 1 + length($tokenHash->{'insertedSequence'})) % 3 == $proteinCodingTokenHash2{'codingIndex2'} % 3) {
										push(@inframeProteinCodingTokenHashList2, \%proteinCodingTokenHash2);
									}
								}
							}
							if(@inframeProteinCodingTokenHashList2) {
								push(@preferredConsequenceTokenHashList, {%$tokenHash, %proteinCodingTokenHash1, %$_, 'consequenceType' => 'inframe'}) foreach(@inframeProteinCodingTokenHashList2);
							} else {
								push(@fallbackConsequenceTokenHashList, {%$tokenHash, %proteinCodingTokenHash1, 'consequenceType' => 'partial_coding'});
							}
						} else {
							push(@fallbackConsequenceTokenHashList, {%$tokenHash, %proteinCodingTokenHash1, 'consequenceType' => 'partial_coding'});
						}
					}
				}
				if($tokenHash->{'breakpointPosition1'} < $proteinCodingTokenHashList1[0]->{'codingStart1'}) {
					if($tokenHash->{'breakpointStrand2'} eq '+' && defined($transcriptProteinCodingListHash{$tokenHash->{'sequenceId2'}})) {
						my @proteinCodingTokenHashList2 = getSequenceIdProteinCodingTokenHashListHash2($tokenHash->{'sequenceId2'});
						foreach(@proteinCodingTokenHashList2) {
							my %proteinCodingTokenHash2 = %$_;
							if($tokenHash->{'breakpointPosition2'} < $proteinCodingTokenHash2{'codingStart2'}) {
								push(@preferredConsequenceTokenHashList, {%$tokenHash, %{$proteinCodingTokenHashList1[0]}, %proteinCodingTokenHash2, 'consequenceType' => '5utr_swap'});
								last;
							}
						}
					}
				}
				if($tokenHash->{'breakpointPosition1'} > $proteinCodingTokenHashList1[0]->{'codingEnd1'}) {
					if($tokenHash->{'breakpointStrand2'} eq '+' && defined($transcriptProteinCodingListHash{$tokenHash->{'sequenceId2'}})) {
						my @proteinCodingTokenHashList2 = getSequenceIdProteinCodingTokenHashListHash2($tokenHash->{'sequenceId2'});
						foreach(@proteinCodingTokenHashList2) {
							my %proteinCodingTokenHash2 = %$_;
							if($tokenHash->{'breakpointPosition2'} < $proteinCodingTokenHash2{'codingStart2'}) {
								push(@fallbackConsequenceTokenHashList, {%$tokenHash, %{$proteinCodingTokenHashList1[0]}, %proteinCodingTokenHash2, 'consequenceType' => '3utr_swap'});
								last;
							}
						}
					}
				}
			}
		}
		if(@preferredConsequenceTokenHashList) {
			if(my @identicalTranscriptTokenHashList = grep {$_->{'sequenceId1'} eq $_->{'sequenceId2'}} @preferredConsequenceTokenHashList) {
				if($keepSelfFusion) {
					@preferredConsequenceTokenHashList = @identicalTranscriptTokenHashList;
				} else {
					next;
				}
			}
			if(my @identicalProteinTokenHashList = grep {defined($_->{'gene1'}) && defined($_->{'gene2'}) && $_->{'gene1'} eq $_->{'gene2'}} @preferredConsequenceTokenHashList) {
				if($keepSelfFusion) {
					@preferredConsequenceTokenHashList = @identicalProteinTokenHashList;
				} else {
					next;
				}
			}
			foreach my $tokenHash (@preferredConsequenceTokenHashList) {
				forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList1, 'fusionReadCount'}), "\n");
			}
		} elsif(@fallbackConsequenceTokenHashList) {
			if(my @identicalTranscriptTokenHashList = grep {$_->{'sequenceId1'} eq $_->{'sequenceId2'}} @fallbackConsequenceTokenHashList) {
				@fallbackConsequenceTokenHashList = @identicalTranscriptTokenHashList;
			}
			if(my @identicalProteinTokenHashList = grep {defined($_->{'gene1'}) && defined($_->{'gene2'}) && $_->{'gene1'} eq $_->{'gene2'}} @fallbackConsequenceTokenHashList) {
				@fallbackConsequenceTokenHashList = @identicalProteinTokenHashList;
			}
			foreach my $tokenHash (@fallbackConsequenceTokenHashList) {
				forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList1, 'fusionReadCount'}), "\n");
			}
		}
	}
}

sub printFusionProtein2 {
	my (@tokenHashList) = @_;
	my %matchLengthListHash = ();
	my %proteinFeaturesListHash = ();
	my %readCountDepthBreadthHash = ();
	foreach my $tokenHash (@tokenHashList) {
		if($tokenHash->{'consequenceType'} eq 'inframe' || $tokenHash->{'consequenceType'} eq 'partial_coding') {
			my $fusionProteinSequence = '';
			if($tokenHash->{'consequenceType'} eq 'inframe') {
				my $codingSequence1 = getPositionSequence($tokenHash->{'sequenceId1'}, my @codingPositionList1 = getCodingPositionList(@$tokenHash{'codingRegions1', 'frame1'}));
				my $codingSequence2 = getPositionSequence($tokenHash->{'sequenceId2'}, my @codingPositionList2 = getCodingPositionList(@$tokenHash{'codingRegions2', 'frame2'}));
				$fusionProteinSequence = translate(my $fusionCodingSequence = join('', substr($codingSequence1, 0, $tokenHash->{'codingIndex1'} + 1), $tokenHash->{'insertedSequence'}, substr($codingSequence2, $tokenHash->{'codingIndex2'})));
				next if($fusionProteinSequence =~ /\*./);
			}
			if($tokenHash->{'consequenceType'} eq 'partial_coding') {
				my $codingSequence1 = getPositionSequence($tokenHash->{'sequenceId1'}, my @codingPositionList1 = getCodingPositionList(@$tokenHash{'codingRegions1', 'frame1'}));
				$fusionProteinSequence = translate(my $fusionCodingSequence = join('', substr($codingSequence1, 0, $tokenHash->{'codingIndex1'} + 1), $tokenHash->{'insertedSequence'}));
				next if($fusionProteinSequence =~ /\*/);
				if($tokenHash->{'breakpointStrand2'} eq '+') {
					my $position = $tokenHash->{'breakpointPosition2'} + (3 - ($tokenHash->{'codingIndex1'} + 1 + length($tokenHash->{'insertedSequence'})) % 3) - 1;
					$fusionCodingSequence .= getSequence($tokenHash->{'sequenceId2'}, $tokenHash->{'breakpointPosition2'}, $position);
					$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
					next if($fusionProteinSequence =~ /\*$/);
					while($fusionProteinSequence !~ /\*$/ && $position + 3 <= $seqHash{$tokenHash->{'sequenceId2'}}->length()) {
						$position += 3;
						$fusionCodingSequence .= getSequence($tokenHash->{'sequenceId2'}, $position - 2, $position);
						$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
					}
					$tokenHash->{'codingEnd2'} = $position;
				}
				if($tokenHash->{'breakpointStrand2'} eq '-') {
					my $position = $tokenHash->{'breakpointPosition2'} - (3 - ($tokenHash->{'codingIndex1'} + 1 + length($tokenHash->{'insertedSequence'})) % 3) + 1;
					$fusionCodingSequence .= getSequence($tokenHash->{'sequenceId2'}, $tokenHash->{'breakpointPosition2'}, $position);
					$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
					next if($fusionProteinSequence =~ /\*$/);
					while($fusionProteinSequence !~ /\*$/ && $position - 3 >= 1) {
						$position -= 3;
						$fusionCodingSequence .= getSequence($tokenHash->{'sequenceId2'}, $position + 2, $position);
						$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
					}
					$tokenHash->{'codingEnd2'} = $position;
				}
			}
			$tokenHash->{'fusionProteinSequence'} = $fusionProteinSequence;
			my $fusionProteinSequenceLength = scalar(my @fusionAAList = split(//, $fusionProteinSequence));
			if($tokenHash->{'proteinId1'} ne '') {
				die "Protein sequence not found: $tokenHash->{'proteinId1'}\n" unless(defined(my $proteinSequence1 = $proteinSequenceHash{$tokenHash->{'proteinId1'}}));
				unless(defined(my $matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId1'}}->{$fusionProteinSequence})) {
					my $proteinSequenceLength = scalar(my @aaList = split(//, $proteinSequence1));
					my $proteinMatchLengthFromN = 0;
					$proteinMatchLengthFromN += 1 while($proteinMatchLengthFromN < $proteinSequenceLength && $proteinMatchLengthFromN < $fusionProteinSequenceLength && $aaList[$proteinMatchLengthFromN] eq $fusionAAList[$proteinMatchLengthFromN]);
					my $proteinMatchLengthFromC = 0;
					$proteinMatchLengthFromC += 1 while($proteinMatchLengthFromC < $proteinSequenceLength && $proteinMatchLengthFromC < $fusionProteinSequenceLength && $aaList[-1 - $proteinMatchLengthFromC] eq $fusionAAList[-1 - $proteinMatchLengthFromC]);
					$matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId1'}}->{$fusionProteinSequence} = [$proteinMatchLengthFromN, $proteinMatchLengthFromC, ($fusionProteinSequenceLength > 0 ? ($proteinMatchLengthFromN + $proteinMatchLengthFromC) / $fusionProteinSequenceLength : 'NA')];
					@$tokenHash{'proteinMatchLengthFromN1', 'proteinMatchLengthFromC1', 'proteinMatchFraction1'} = @$matchLengthList;
				} else {
					@$tokenHash{'proteinMatchLengthFromN1', 'proteinMatchLengthFromC1', 'proteinMatchFraction1'} = @$matchLengthList;
				}
				unless(defined(my $proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId1'}}->{$tokenHash->{'proteinMatchLengthFromN1'}}->{$tokenHash->{'proteinMatchLengthFromC1'}})) {
					my @proteinFeatureListList = ([], [], []);
					if(defined(my $proteinFeatureList = $proteinFeatureListHash{$tokenHash->{'proteinId1'}})) {
						my $matchEndN = $tokenHash->{'proteinMatchLengthFromN1'};
						my $matchStartC = length($proteinSequence1) - $tokenHash->{'proteinMatchLengthFromC1'} + 1;
						push(@{$proteinFeatureListList[$_->[0]]}, $_->[1]) foreach(map {[($_->[1] <= $matchEndN || $matchStartC <= $_->[0] ? 0 : $_->[0] <= $matchEndN || $matchStartC <= $_->[1] ? 1 : 2), $_->[2]]} @$proteinFeatureList);
					}
					$proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId1'}}->{$tokenHash->{'proteinMatchLengthFromN1'}}->{$tokenHash->{'proteinMatchLengthFromC1'}} = [map {join(' | ', @$_)} @proteinFeatureListList];
					@$tokenHash{'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1'} = @$proteinFeaturesList;
				} else {
					@$tokenHash{'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1'} = @$proteinFeaturesList;
				}
			}
			if($tokenHash->{'proteinId2'} ne '') {
				die "Protein sequence not found: $tokenHash->{'proteinId2'}\n" unless(defined(my $proteinSequence2 = $proteinSequenceHash{$tokenHash->{'proteinId2'}}));
				unless(defined(my $matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId2'}}->{$fusionProteinSequence})) {
					my $proteinSequenceLength = scalar(my @aaList = split(//, $proteinSequence2));
					my $proteinMatchLengthFromN = 0;
					$proteinMatchLengthFromN += 1 while($proteinMatchLengthFromN < $proteinSequenceLength && $proteinMatchLengthFromN < $fusionProteinSequenceLength && $aaList[$proteinMatchLengthFromN] eq $fusionAAList[$proteinMatchLengthFromN]);
					my $proteinMatchLengthFromC = 0;
					$proteinMatchLengthFromC += 1 while($proteinMatchLengthFromC < $proteinSequenceLength && $proteinMatchLengthFromC < $fusionProteinSequenceLength && $aaList[-1 - $proteinMatchLengthFromC] eq $fusionAAList[-1 - $proteinMatchLengthFromC]);
					$matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId2'}}->{$fusionProteinSequence} = [$proteinMatchLengthFromN, $proteinMatchLengthFromC, ($fusionProteinSequenceLength > 0 ? ($proteinMatchLengthFromN + $proteinMatchLengthFromC) / $fusionProteinSequenceLength : 'NA')];
					@$tokenHash{'proteinMatchLengthFromN2', 'proteinMatchLengthFromC2', 'proteinMatchFraction2'} = @$matchLengthList;
				} else {
					@$tokenHash{'proteinMatchLengthFromN2', 'proteinMatchLengthFromC2', 'proteinMatchFraction2'} = @$matchLengthList;
				}
				unless(defined(my $proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId2'}}->{$tokenHash->{'proteinMatchLengthFromN2'}}->{$tokenHash->{'proteinMatchLengthFromC2'}})) {
					my @proteinFeatureListList = ([], [], []);
					if(defined(my $proteinFeatureList = $proteinFeatureListHash{$tokenHash->{'proteinId2'}})) {
						my $matchEndN = $tokenHash->{'proteinMatchLengthFromN2'};
						my $matchStartC = length($proteinSequence2) - $tokenHash->{'proteinMatchLengthFromC2'} + 1;
						push(@{$proteinFeatureListList[$_->[0]]}, $_->[1]) foreach(map {[($_->[1] <= $matchEndN || $matchStartC <= $_->[0] ? 0 : $_->[0] <= $matchEndN || $matchStartC <= $_->[1] ? 1 : 2), $_->[2]]} @$proteinFeatureList);
					}
					$proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId2'}}->{$tokenHash->{'proteinMatchLengthFromN2'}}->{$tokenHash->{'proteinMatchLengthFromC2'}} = [map {join(' | ', @$_)} @proteinFeatureListList];
					@$tokenHash{'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2'} = @$proteinFeaturesList;
				} else {
					@$tokenHash{'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2'} = @$proteinFeaturesList;
				}
			}
		}
		if(@regionReadFileList) {
			my $key = join("\t", my ($sequenceId, $position, $strand, $start, $end) = (@$tokenHash{'sequenceId1', 'breakpointPosition1', 'breakpointStrand1'}, sort {$a <=> $b} @$tokenHash{'codingStart1', 'breakpointPosition1'}));
			unless(defined($readCountDepthBreadthHash{$key})) {
				my %readNameHash = ();
				my %positionDepthHash = ();
				foreach my $regionReadFile (@regionReadFileList) {
					open(my $reader, "tabix $regionReadFile $sequenceId:$start-$end |");
					while(my $line = <$reader>) {
						chomp($line);
						my ($sequenceId, $readStart, $readEnd, $readStrand, $readName) = split(/\t/, $line, -1);
						next if($stranded eq 'f' && $readStrand ne $strand);
						next if($stranded eq 'r' && $readStrand eq $strand);
						$readNameHash{$readName} = 1 if($readStart <= $position && $position <= $readEnd);
						$positionDepthHash{$_} += 1 foreach($readStart .. $readEnd);
					}
					close($reader);
				}
				my @readNameList = keys %readNameHash;
				my $readCount = scalar(@readNameList);
				my @depthList = grep {defined} @positionDepthHash{$start .. $end};
				my $meanDepth = sum(0, @depthList) / (my $length = $end - $start + 1);
				my $breadth = scalar(grep {$_ >= $minimumBreadthDepth} @depthList) / $length;
				$readCountDepthBreadthHash{$key} = [$readCount, $meanDepth, $breadth];
			}
			my ($readCount, $meanDepth, $breadth) = @{$readCountDepthBreadthHash{$key}};
			@$tokenHash{'fusionReadFraction1', 'meanDepth1', 'breadth1'} = (($readCount > 0 ? ($tokenHash->{'fusionReadCount'} / $readCount) : 'NA'), $meanDepth, $breadth);
		}
		if(@regionReadFileList) {
			my $key = join("\t", my ($sequenceId, $position, $strand, $start, $end) = (@$tokenHash{'sequenceId2', 'breakpointPosition2', 'breakpointStrand2'}, sort {$a <=> $b} @$tokenHash{'breakpointPosition2', 'codingEnd2'}));
			unless(defined($readCountDepthBreadthHash{$key})) {
				my %readNameHash = ();
				my %positionDepthHash = ();
				foreach my $regionReadFile (@regionReadFileList) {
					open(my $reader, "tabix $regionReadFile $sequenceId:$start-$end |");
					while(my $line = <$reader>) {
						chomp($line);
						my ($sequenceId, $readStart, $readEnd, $readStrand, $readName) = split(/\t/, $line, -1);
						next if($stranded eq 'f' && $readStrand ne $strand);
						next if($stranded eq 'r' && $readStrand eq $strand);
						$readNameHash{$readName} = 1 if($readStart <= $position && $position <= $readEnd);
						$positionDepthHash{$_} += 1 foreach($readStart .. $readEnd);
					}
					close($reader);
				}
				my @readNameList = keys %readNameHash;
				my $readCount = scalar(@readNameList);
				my @depthList = grep {defined} @positionDepthHash{$start .. $end};
				my $meanDepth = sum(0, @depthList) / (my $length = $end - $start + 1);
				my $breadth = scalar(grep {$_ >= $minimumBreadthDepth} @depthList) / $length;
				$readCountDepthBreadthHash{$key} = [$readCount, $meanDepth, $breadth];
			}
			my ($readCount, $meanDepth, $breadth) = @{$readCountDepthBreadthHash{$key}};
			@$tokenHash{'fusionReadFraction2', 'meanDepth2', 'breadth2'} = (($readCount > 0 ? ($tokenHash->{'fusionReadCount'} / $readCount) : 'NA'), $meanDepth, $breadth);
		}
		forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList2}), "\n");
	}
}

sub getPositionSequence {
	my ($sequenceId, @positionList) = @_;
	my @startList = ();
	my @endList = ();
	foreach my $index (1 .. $#positionList) {
		my ($end, $start) = @positionList[$index - 1, $index];
		if($end + 1 != $start) {
			push(@startList, $start);
			push(@endList, $end);
		}
	}
	@startList = ($positionList[0], @startList);
	@endList = (@endList, $positionList[-1]);
	return join('', map {getSequence($sequenceId, $startList[$_], $endList[$_])} 0 .. $#startList);
}

sub getCodingPositionList {
	my ($codingRegions, $frame) = @_;
	my @codingPositionList = eval($codingRegions);
	die "Invalid codingRegions: $codingRegions\n$@" if($@);
	@codingPositionList = @codingPositionList[$frame .. $#codingPositionList] if($frame > 0);
	return @codingPositionList;
}

sub getSequence {
	my ($sequenceId, $start, $end) = @_;
	return $end < $start ? getReverseComplementarySequence(uc($seqHash{$sequenceId}->subseq($end, $start))) : uc($seqHash{$sequenceId}->subseq($start, $end));
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGTacgt/TGCAtgca/;
	return $sequence;
}
