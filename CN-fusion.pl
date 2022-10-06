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
system("mkdir -p $dataPath");

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 10000),
	's' => \(my $sort = ''),
	'S=s' => \(my $stranded = ''),
	'C=s' => \@codonList,
	'c=i' => \(my $minimumCount = 3),
	'd=i' => \(my $coverageDepth = 1),
	'r=s' => \(my $regionReadFile = ''),
	'fastaFile=s' => \(my $fastaFile = "$dataPath/rna_genomic.fasta"),
	'proteinFastaFile=s' => \(my $proteinFastaFile = "$dataPath/protein.fasta"),
	'proteinCodingFile=s' => \(my $proteinCodingFile = "$dataPath/protein_coding.txt"),
	'proteinFeatureFile=s' => \(my $proteinFeatureFile = "$dataPath/protein_feature.txt"),
	'spliceJunctionFile=s' => \(my $spliceJunctionFile = "$dataPath/splice_junction.txt"),
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
         -d INT   coverage depth [$coverageDepth]
         -r FILE  region read file

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
my $db = Bio::DB::Fasta->new($fastaFile);
my %transcriptProteinCodingHash = ();
{
	open(my $reader, ($proteinCodingFile =~ /\.gz$/ ? "gzip -dc $proteinCodingFile |" : $proteinCodingFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($gene, $transcriptId, $proteinId, $codingRegions, $frame) = split(/\t/, $line, -1);
		$transcriptProteinCodingHash{$transcriptId} = [$gene, $proteinId, $codingRegions, $frame];
	}
	close($reader);
}
my %proteinSequenceHash = ();
{
	my $proteinId = '';
	open(my $reader, ($proteinFastaFile =~ /\.gz$/ ? "gzip -dc $proteinFastaFile |" : $proteinFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($proteinId = $1));
		$proteinSequenceHash{$proteinId} .= $line;
	}
	close($reader);
	s/U/*/g foreach(values %proteinSequenceHash);
	s/\*?$/*/ foreach(values %proteinSequenceHash);
}
my %proteinFeatureListHash = ();
{
	open(my $reader, ($proteinFeatureFile =~ /\.gz$/ ? "gzip -dc $proteinFeatureFile |" : $proteinFeatureFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($proteinId, $start, $end, $proteinFeature) = split(/\t/, $line, -1);
		push(@{$proteinFeatureListHash{$proteinId}}, [$start, $end, $proteinFeature]);
	}
	close($reader);
}
my %transcriptPositionSpliceJunctionHash = ();
{
	open(my $reader, ($spliceJunctionFile =~ /\.gz$/ ? "gzip -dc $spliceJunctionFile |" : $spliceJunctionFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($transcriptId, $position, $spliceJunction1, $spliceJunction2) = split(/\t/, $line, -1);
		next if(defined($transcriptPositionSpliceJunctionHash{$transcriptId}->{$position}));
		$transcriptPositionSpliceJunctionHash{$transcriptId}->{$position} = [$spliceJunction1, $spliceJunction2];
	}
	close($reader);
}
my @fusionColumnList = ('chromosome1', 'position1', 'strand1', 'chromosome2', 'position2', 'strand2', 'intermediateSequence', 'apart');
my @fusionProteinCodingColumnList1 = (@fusionColumnList, 'type', 'spliceJunction1', 'spliceJunction2');
{
	push(@fusionProteinCodingColumnList1, 'gene1', 'proteinId1', 'codingRegions1', 'frame1', 'index1');
	push(@fusionProteinCodingColumnList1, 'gene2', 'proteinId2', 'codingRegions2', 'frame2', 'index2');
}
my @fusionProteinCodingColumnList2 = (@fusionColumnList, 'count', 'type', 'spliceJunction1', 'spliceJunction2', 'codingStart1', 'codingEnd2');
if($regionReadFile ne '') {
	push(@fusionProteinCodingColumnList2, 'gene1', 'proteinId1', 'matchLengthN1', 'matchLengthC1', 'matchCoverage1', 'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1', 'ratio1', 'depth1', 'coverage1');
	push(@fusionProteinCodingColumnList2, 'gene2', 'proteinId2', 'matchLengthN2', 'matchLengthC2', 'matchCoverage2', 'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2', 'ratio2', 'depth2', 'coverage2');
} else {
	push(@fusionProteinCodingColumnList2, 'gene1', 'proteinId1', 'matchLengthN1', 'matchLengthC1', 'matchCoverage1', 'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1');
	push(@fusionProteinCodingColumnList2, 'gene2', 'proteinId2', 'matchLengthN2', 'matchLengthC2', 'matchCoverage2', 'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2');
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
				if($tokenHash{'strand1'} eq '+') {
					if($stranded eq '' || ($stranded eq 'f' && $tokenHash{'strand'} eq '+') || ($stranded eq 'r' && $tokenHash{'strand'} eq '-')) {
						push(@tokenList, @tokenHash{@fusionColumnList});
					}
				}
				if($tokenHash{'strand2'} eq '-') {
					@tokenHash{'chromosome1', 'position1', 'strand1', 'chromosome2', 'position2', 'strand2', 'intermediateSequence', 'strand'} = (@tokenHash{'chromosome2', 'position2'}, '+', @tokenHash{'chromosome1', 'position1'}, ($tokenHash{'strand1'} eq '+' ? '-' : $tokenHash{'strand1'} eq '-' ? '+' : ''), getReverseComplementarySequence($tokenHash{'intermediateSequence'}), ($tokenHash{'strand'} eq '+' ? '-' : $tokenHash{'strand'} eq '-' ? '+' : ''));
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
					@tokenHash{@fusionColumnList, 'count'} = (splice(@tokenList, 0, scalar(@fusionColumnList)), $count);
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
			@tokenHash{@fusionProteinCodingColumnList1, 'count'} = split(/\t/, $line, -1);
			if((my $currentFusionGene = join("\t", @tokenHash{@fusionProteinCodingColumnList1})) ne $fusionGene) {
				if($fusionGene ne '' && $count >= $minimumCount) {
					my %tokenHash = ();
					@tokenHash{@fusionProteinCodingColumnList1, 'count'} = (split(/\t/, $fusionGene, -1), $count);
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
			$count += $tokenHash{'count'};
		}
		if($fusionGene ne '' && $count >= $minimumCount) {
			my %tokenHash = ();
			@tokenHash{@fusionProteinCodingColumnList1, 'count'} = (split(/\t/, $fusionGene, -1), $count);
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

sub printFusionProtein1 {
	my (@tokenHashListList) = @_;
	my %transcriptCodingPositionIndexHash1 = ();
	my %transcriptCodingPositionIndexHash2 = ();
	foreach(@tokenHashListList) {
		my @tokenHashList = @$_;
		my @codingTokenHashList = ();
		my @noncodingTokenHashList = ();
		foreach my $tokenHash (@tokenHashList) {
			if((my $intermediateSequenceLength = length($tokenHash->{'intermediateSequence'})) > 0 && $tokenHash->{'apart'} eq '' && $tokenHash->{'strand1'} eq '+' && $tokenHash->{'strand2'} eq '+') {
				if(defined(my $positionSpliceJunctionHash1 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'chromosome1'}}) && defined(my $positionSpliceJunctionHash2 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'chromosome2'}})) {
					foreach my $index (0 .. $intermediateSequenceLength) {
						my $position1 = $tokenHash->{'position1'} + $index;
						my $position2 = $tokenHash->{'position2'} - ($intermediateSequenceLength - $index);
						if(defined(my $spliceJunction1 = $positionSpliceJunctionHash1->{$position1}->[0]) && defined(my $spliceJunction2 = $positionSpliceJunctionHash2->{$position2 - 1}->[1])) {
							$tokenHash->{'position1'} = $position1;
							$tokenHash->{'position2'} = $position2;
							$tokenHash->{'spliceJunction1'} = $spliceJunction1;
							$tokenHash->{'spliceJunction2'} = $spliceJunction2;
							$tokenHash->{'intermediateSequence'} = '';
							last;
						}
					}
				}
			}
			if((my $intermediateSequenceLength = length($tokenHash->{'intermediateSequence'})) > 0 && $tokenHash->{'apart'} eq '' && $tokenHash->{'strand1'} eq '+') {
				if(defined(my $positionSpliceJunctionHash1 = $transcriptPositionSpliceJunctionHash{$tokenHash->{'chromosome1'}})) {
					foreach my $index (0 .. $intermediateSequenceLength) {
						my $position1 = $tokenHash->{'position1'} + $index;
						if(defined(my $spliceJunction1 = $positionSpliceJunctionHash1->{$position1}->[0])) {
							$tokenHash->{'position1'} = $position1;
							$tokenHash->{'spliceJunction1'} = $spliceJunction1;
							$tokenHash->{'intermediateSequence'} = substr($tokenHash->{'intermediateSequence'}, $index);
							last;
						}
					}
				}
			}
			if($tokenHash->{'strand1'} eq '+' && defined(my $proteinCoding1 = $transcriptProteinCodingHash{$tokenHash->{'chromosome1'}})) {
				@$tokenHash{'gene1', 'proteinId1', 'codingRegions1', 'frame1'} = @$proteinCoding1;
				my $codingPositionIndexHash1;
				unless(defined($codingPositionIndexHash1 = $transcriptCodingPositionIndexHash1{$tokenHash->{'chromosome1'}})) {
					my @codingPositionList1 = getCodingPositionList(@$tokenHash{'codingRegions1', 'frame1'});
					my %codingPositionIndexHash1 = map {$codingPositionList1[$_] => $_} 0 .. $#codingPositionList1;
					$codingPositionIndexHash1 = $transcriptCodingPositionIndexHash1{$tokenHash->{'chromosome1'}} = \%codingPositionIndexHash1;
				}
				if(defined($tokenHash->{'index1'} = $codingPositionIndexHash1->{$tokenHash->{'position1'}})) {
					if($tokenHash->{'strand2'} eq '+' && defined(my $proteinCoding2 = $transcriptProteinCodingHash{$tokenHash->{'chromosome2'}})) {
						@$tokenHash{'gene2', 'proteinId2', 'codingRegions2', 'frame2'} = @$proteinCoding2;
						my $codingPositionIndexHash2;
						unless(defined($codingPositionIndexHash2 = $transcriptCodingPositionIndexHash2{$tokenHash->{'chromosome2'}})) {
							my @codingPositionList2 = getCodingPositionList(@$tokenHash{'codingRegions2', 'frame2'});
							splice(@codingPositionList2, -3) if(translate(getSequence($tokenHash->{'chromosome2'}, @codingPositionList2[-3 .. -1])) eq '*');
							my %codingPositionIndexHash2 = map {$codingPositionList2[$_] => $_} 0 .. $#codingPositionList2;
							$codingPositionIndexHash2 = $transcriptCodingPositionIndexHash2{$tokenHash->{'chromosome2'}} = \%codingPositionIndexHash2;
						}
						if(defined($tokenHash->{'index2'} = $codingPositionIndexHash2->{$tokenHash->{'position2'}})) {
							if(($tokenHash->{'index1'} + 1 + length($tokenHash->{'intermediateSequence'})) % 3 == $tokenHash->{'index2'} % 3) {
								push(@codingTokenHashList, $tokenHash);
							} else {
								push(@noncodingTokenHashList, $tokenHash);
							}
						} else {
							push(@noncodingTokenHashList, $tokenHash);
						}
					} else {
						push(@noncodingTokenHashList, $tokenHash);
					}
				}
			}
		}
		if(@codingTokenHashList) {
			if(my @identicalTranscriptCodingTokenHashList = grep {$_->{'chromosome1'} eq $_->{'chromosome2'}} @codingTokenHashList) {
				@codingTokenHashList = @identicalTranscriptCodingTokenHashList;
				next;
			}
			if(my @identicalProteinCodingTokenHashList = grep {defined($_->{'gene1'}) && defined($_->{'gene2'}) && $_->{'gene1'} eq $_->{'gene2'}} @codingTokenHashList) {
				@codingTokenHashList = @identicalProteinCodingTokenHashList;
				next;
			}
			foreach my $tokenHash (@codingTokenHashList) {
				$tokenHash->{'type'} = 'coding';
				forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList1, 'count'}), "\n");
			}
		} elsif(@noncodingTokenHashList) {
			if(my @identicalTranscriptCodingNoncodingTokenHashList = grep {$_->{'chromosome1'} eq $_->{'chromosome2'}} @noncodingTokenHashList) {
				@noncodingTokenHashList = @identicalTranscriptCodingNoncodingTokenHashList;
			}
			if(my @identicalProteinCodingNoncodingTokenHashList = grep {defined($_->{'gene1'}) && defined($_->{'gene2'}) && $_->{'gene1'} eq $_->{'gene2'}} @noncodingTokenHashList) {
				@noncodingTokenHashList = @identicalProteinCodingNoncodingTokenHashList;
			}
			foreach my $tokenHash (grep {$_->{'chromosome2'} ne ''} @noncodingTokenHashList) {
				$tokenHash->{'type'} = 'noncoding';
				forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList1, 'count'}), "\n");
			}
		}
	}
}

sub printFusionProtein2 {
	my (@tokenHashList) = @_;
	my %matchLengthListHash = ();
	my %proteinFeaturesListHash = ();
	my %readCountDepthCoverageHash = ();
	foreach my $tokenHash (@tokenHashList) {
		my $fusionProteinSequence = '';
		if($tokenHash->{'type'} eq 'coding') {
			my $codingSequence1 = getSequence($tokenHash->{'chromosome1'}, my @codingPositionList1 = getCodingPositionList(@$tokenHash{'codingRegions1', 'frame1'}));
			$tokenHash->{'codingStart1'} = $codingPositionList1[0];
			my $codingSequence2 = getSequence($tokenHash->{'chromosome2'}, my @codingPositionList2 = getCodingPositionList(@$tokenHash{'codingRegions2', 'frame2'}));
			$tokenHash->{'codingEnd2'} = $codingPositionList2[-1];
			$fusionProteinSequence = translate(my $fusionCodingSequence = join('', substr($codingSequence1, 0, $tokenHash->{'index1'} + 1), $tokenHash->{'intermediateSequence'}, substr($codingSequence2, $tokenHash->{'index2'})));
			next if($fusionProteinSequence =~ /\*./);
		}
		if($tokenHash->{'type'} eq 'noncoding') {
			my $codingSequence1 = getSequence($tokenHash->{'chromosome1'}, my @codingPositionList1 = getCodingPositionList(@$tokenHash{'codingRegions1', 'frame1'}));
			$tokenHash->{'codingStart1'} = $codingPositionList1[0];
			$fusionProteinSequence = translate(my $fusionCodingSequence = join('', substr($codingSequence1, 0, $tokenHash->{'index1'} + 1), $tokenHash->{'intermediateSequence'}));
			next if($fusionProteinSequence =~ /\*/);
			if($tokenHash->{'strand2'} eq '+') {
				my $position = $tokenHash->{'position2'} + (3 - ($tokenHash->{'index1'} + 1 + length($tokenHash->{'intermediateSequence'})) % 3);
				$fusionCodingSequence .= uc($db->seq($tokenHash->{'chromosome2'}, $tokenHash->{'position2'} + 1, $position));
				$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
				next if($fusionProteinSequence =~ /\*$/);
				while($fusionProteinSequence !~ /\*$/ && $position + 3 <= $db->length($tokenHash->{'chromosome2'})) {
					$position += 3;
					$fusionCodingSequence .= uc($db->seq($tokenHash->{'chromosome2'}, $position - 2, $position));
					$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
				}
				$tokenHash->{'codingEnd2'} = $position;
			}
			if($tokenHash->{'strand2'} eq '-') {
				my $position = $tokenHash->{'position2'} - (3 - ($tokenHash->{'index1'} + 1 + length($tokenHash->{'intermediateSequence'})) % 3);
				$fusionCodingSequence .= uc($db->seq($tokenHash->{'chromosome2'}, $tokenHash->{'position2'} - 1, $position));
				$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
				next if($fusionProteinSequence =~ /\*$/);
				while($fusionProteinSequence !~ /\*$/ && $position - 3 >= 1) {
					$position -= 3;
					$fusionCodingSequence .= uc($db->seq($tokenHash->{'chromosome2'}, $position + 2, $position));
					$fusionProteinSequence .= translate(substr($fusionCodingSequence, -3));
				}
				$tokenHash->{'codingEnd2'} = $position;
			}
		}
		$tokenHash->{'fusionProteinSequence'} = $fusionProteinSequence;
		my $fusionProteinSequenceLength = scalar(my @fusionAAList = split(//, $fusionProteinSequence));
		if($tokenHash->{'proteinId1'} ne '') {
			unless(defined(my $matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId1'}}->{$fusionProteinSequence})) {
				my $proteinSequenceLength = scalar(my @aaList = split(//, $proteinSequenceHash{$tokenHash->{'proteinId1'}}));
				my $matchLengthN = $tokenHash->{'index1'} ne '' ? int(($tokenHash->{'index1'} + 1) / 3) : 0;
				$matchLengthN += 1 while($matchLengthN < $proteinSequenceLength && $matchLengthN < $fusionProteinSequenceLength && $aaList[$matchLengthN] eq $fusionAAList[$matchLengthN]);
				my $matchLengthC = 0;
				$matchLengthC += 1 while($matchLengthC < $proteinSequenceLength && $matchLengthC < $fusionProteinSequenceLength && $aaList[-1 - $matchLengthC] eq $fusionAAList[-1 - $matchLengthC]);
				$matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId1'}}->{$fusionProteinSequence} = [$matchLengthN, $matchLengthC, ($matchLengthN + $matchLengthC) / $fusionProteinSequenceLength];
				@$tokenHash{'matchLengthN1', 'matchLengthC1', 'matchCoverage1'} = @$matchLengthList;
			} else {
				@$tokenHash{'matchLengthN1', 'matchLengthC1', 'matchCoverage1'} = @$matchLengthList;
			}
			unless(defined(my $proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId1'}}->{$tokenHash->{'matchLengthN1'}}->{$tokenHash->{'matchLengthC1'}})) {
				my @proteinFeatureListList = ([], [], []);
				if(defined(my $proteinFeatureList = $proteinFeatureListHash{$tokenHash->{'proteinId1'}})) {
					push(@{$proteinFeatureListList[$_->[0]]}, $_->[1]) foreach(map {[($_->[1] <= $tokenHash->{'matchLengthN1'} || $tokenHash->{'matchLengthC1'} <= $_->[0] ? 0 : $_->[0] <= $tokenHash->{'matchLengthN1'} || $tokenHash->{'matchLengthC1'} <= $_->[1] ? 1 : 2), $_->[2]]} @$proteinFeatureList);
				}
				$proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId1'}}->{$tokenHash->{'matchLengthN1'}}->{$tokenHash->{'matchLengthC1'}} = [map {defined($_) ? join(' | ', @$_) : ''} @proteinFeatureListList];
				@$tokenHash{'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1'} = @$proteinFeaturesList;
			} else {
				@$tokenHash{'includedProteinFeatures1', 'partiallyIncludedProteinFeatures1', 'excludedProteinFeatures1'} = @$proteinFeaturesList;
			}
		}
		if($tokenHash->{'proteinId2'} ne '') {
			unless(defined(my $matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId2'}}->{$fusionProteinSequence})) {
				my $proteinSequenceLength = scalar(my @aaList = split(//, $proteinSequenceHash{$tokenHash->{'proteinId2'}}));
				my $matchLengthN = 0;
				$matchLengthN += 1 while($matchLengthN < $proteinSequenceLength && $matchLengthN < $fusionProteinSequenceLength && $aaList[$matchLengthN] eq $fusionAAList[$matchLengthN]);
				my $matchLengthC = $tokenHash->{'index2'} ne '' ? $proteinSequenceLength - int(($tokenHash->{'index2'} + 5) / 3) + 1 : 0;
				$matchLengthC += 1 while($matchLengthC < $proteinSequenceLength && $matchLengthC < $fusionProteinSequenceLength && $aaList[-1 - $matchLengthC] eq $fusionAAList[-1 - $matchLengthC]);
				$matchLengthList = $matchLengthListHash{$tokenHash->{'proteinId2'}}->{$fusionProteinSequence} = [$matchLengthN, $matchLengthC, ($matchLengthN + $matchLengthC) / $fusionProteinSequenceLength];
				@$tokenHash{'matchLengthN2', 'matchLengthC2', 'matchCoverage2'} = @$matchLengthList;
			} else {
				@$tokenHash{'matchLengthN2', 'matchLengthC2', 'matchCoverage2'} = @$matchLengthList;
			}
			unless(defined(my $proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId2'}}->{$tokenHash->{'matchLengthN2'}}->{$tokenHash->{'matchLengthC2'}})) {
				my @proteinFeatureListList = ([], [], []);
				if(defined(my $proteinFeatureList = $proteinFeatureListHash{$tokenHash->{'proteinId2'}})) {
					push(@{$proteinFeatureListList[$_->[0]]}, $_->[1]) foreach(map {[($_->[1] <= $tokenHash->{'matchLengthN2'} || $tokenHash->{'matchLengthC2'} <= $_->[0] ? 0 : $_->[0] <= $tokenHash->{'matchLengthN2'} || $tokenHash->{'matchLengthC2'} <= $_->[1] ? 1 : 2), $_->[2]]} @$proteinFeatureList);
				}
				$proteinFeaturesList = $proteinFeaturesListHash{$tokenHash->{'proteinId2'}}->{$tokenHash->{'matchLengthN2'}}->{$tokenHash->{'matchLengthC2'}} = [map {defined($_) ? join(' | ', @$_) : ''} @proteinFeatureListList];
				@$tokenHash{'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2'} = @$proteinFeaturesList;
			} else {
				@$tokenHash{'includedProteinFeatures2', 'partiallyIncludedProteinFeatures2', 'excludedProteinFeatures2'} = @$proteinFeaturesList;
			}
		}
		if($regionReadFile ne '') {
			my $key = join("\t", my ($chromosome, $position, $strand, $start, $end) = (@$tokenHash{'chromosome1', 'position1', 'strand1'}, sort {$a <=> $b} @$tokenHash{'codingStart1', 'position1'}));
			unless(defined($readCountDepthCoverageHash{$key})) {
				my %readNameHash = ();
				my %positionDepthHash = ();
				open(my $reader, "tabix $regionReadFile $chromosome:$start-$end |");
				while(my $line = <$reader>) {
					chomp($line);
					my ($chromosome, $readStart, $readEnd, $readStrand, $readName) = split(/\t/, $line, -1);
					next if($stranded eq 'f' && $readStrand ne $strand);
					next if($stranded eq 'r' && $readStrand eq $strand);
					$readNameHash{$readName} = 1 if($readStart <= $position && $position <= $readEnd);
					$positionDepthHash{$_} += 1 foreach($readStart .. $readEnd);
				}
				close($reader);
				my @readNameList = keys %readNameHash;
				my $readCount = scalar(@readNameList);
				my @depthList = grep {defined} @positionDepthHash{$start .. $end};
				my $depth = sum(0, @depthList) / (my $length = $end - $start + 1);
				my $coverage = scalar(grep {$_ >= $coverageDepth} @depthList) / $length;
				$readCountDepthCoverageHash{$key} = [$readCount, $depth, $coverage];
			}
			my ($readCount, $depth, $coverage) = @{$readCountDepthCoverageHash{$key}};
			@$tokenHash{'ratio1', 'depth1', 'coverage1'} = ($tokenHash->{'count'} / $readCount, $depth, $coverage);
		}
		if($regionReadFile ne '') {
			my $key = join("\t", my ($chromosome, $position, $strand, $start, $end) = (@$tokenHash{'chromosome2', 'position2', 'strand2'}, sort {$a <=> $b} @$tokenHash{'position2', 'codingEnd2'}));
			unless(defined($readCountDepthCoverageHash{$key})) {
				my %readNameHash = ();
				my %positionDepthHash = ();
				open(my $reader, "tabix $regionReadFile $chromosome:$start-$end |");
				while(my $line = <$reader>) {
					chomp($line);
					my ($chromosome, $readStart, $readEnd, $readStrand, $readName) = split(/\t/, $line, -1);
					next if($stranded eq 'f' && $readStrand ne $strand);
					next if($stranded eq 'r' && $readStrand eq $strand);
					$readNameHash{$readName} = 1 if($readStart <= $position && $position <= $readEnd);
					$positionDepthHash{$_} += 1 foreach($readStart .. $readEnd);
				}
				close($reader);
				my @readNameList = keys %readNameHash;
				my $readCount = scalar(@readNameList);
				my @depthList = grep {defined} @positionDepthHash{$start .. $end};
				my $depth = sum(0, @depthList) / (my $length = $end - $start + 1);
				my $coverage = scalar(grep {$_ >= $coverageDepth} @depthList) / $length;
				$readCountDepthCoverageHash{$key} = [$readCount, $depth, $coverage];
			}
			my ($readCount, $depth, $coverage) = @{$readCountDepthCoverageHash{$key}};
			@$tokenHash{'ratio2', 'depth2', 'coverage2'} = ($tokenHash->{'count'} / $readCount, $depth, $coverage);
		}
		forkPrint(join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@fusionProteinCodingColumnList2}), "\n");
	}
}

sub getSequence {
	my ($chromosome, @positionList) = @_;
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
	return join('', map {uc($db->seq($chromosome, $startList[$_], $endList[$_]))} 0 .. $#startList);
}

sub getCodingPositionList {
	my ($codingRegions, $frame) = @_;
	my @codingPositionList = eval($codingRegions);
	@codingPositionList = @codingPositionList[$frame .. $#codingPositionList] if($frame > 0);
	return @codingPositionList;
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
