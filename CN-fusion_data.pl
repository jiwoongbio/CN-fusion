# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use IPC::Open2;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

my @gpffFileList = ();
GetOptions(
	'h' => \(my $help = ''),
	'p=s' => \@gpffFileList,
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl CN-fusion_data.pl ftp://ftp.ncbi.nlm.nih.gov/genomes/...

Options: -h       display this help message
         -p FILE  GenPept format file for additional protein features

EOF
}
my (@genomeAssemblyFTPList) = @ARGV;
my @genomeAssemblyList = ();
foreach my $genomeAssemblyFTP (@genomeAssemblyFTPList) {
	system("lftp -c 'mirror -p -L $genomeAssemblyFTP $dataPath/'");
	(my $genomeAssembly = $genomeAssemblyFTP) =~ s/^.*\///;
	push(@genomeAssemblyList, $genomeAssembly);
}
{
	my $fastaFile = "$dataPath/rna_genomic.fasta";
	open(my $writer, "> $fastaFile");
	foreach my $genomeAssembly (@genomeAssemblyList) {
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_rna.fna.gz |");
		print $writer $_ while(<$reader>);
		close($reader);
	}
	foreach my $genomeAssembly (@genomeAssemblyList) {
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_genomic.fna.gz |");
		print $writer $_ while(<$reader>);
		close($reader);
	}
	close($writer);
	system("bwa index $fastaFile");
	my $db = Bio::DB::Fasta->new($fastaFile);
}
{
	my $proteinFastaFile = "$dataPath/protein.fasta";
	open(my $writer, "> $proteinFastaFile");
	foreach my $genomeAssembly (@genomeAssemblyList) {
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_protein.faa.gz |");
		print $writer $_ while(<$reader>);
		close($reader);
	}
	close($writer);
}
{
	my $proteinCodingFile = "$dataPath/protein_coding.txt";
	my $terminationCodons = 'TAG,TAA,TGA';
	my %terminationCodonHash = map {$_ => 1} split(/,/, $terminationCodons);
	open(my $writer, "| sort -u > $proteinCodingFile");
	foreach my $genomeAssembly (@genomeAssemblyList) {
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_rna.gbff.gz |");
		my $string = '';
		while(my $line = <$reader>) {
			$string .= $line;
			if($line eq "//\n") {
				my $seqInput_object = Bio::SeqIO->new(-string => $string, -format => 'genbank');
				while(my $seq_object = $seqInput_object->next_seq()) {
					my $transcriptId = join('.', $seq_object->display_id, $seq_object->version);
					my $transcriptSequence = $seq_object->seq;
					my $gene = join(',', map {$_->get_tag_values('gene')} grep {$_->primary_tag eq 'gene'} $seq_object->get_SeqFeatures);
					foreach my $feat_object (grep {$_->primary_tag eq 'CDS'} $seq_object->get_SeqFeatures) {
						my ($proteinId) = $feat_object->get_tag_values('protein_id');
						my ($codonStart) = $feat_object->get_tag_values('codon_start');
						my $codingRegions = join(',', map {join('..', $_->start, $_->end)} $feat_object->location->each_Location);
						if($codingRegions =~ s/([0-9]+)$//) {
							if(!$terminationCodonHash{substr($transcriptSequence, $1 - 3, 3)} && $terminationCodonHash{substr($transcriptSequence, $1, 3)}) {
								$codingRegions .= $1 + 3;
							} else {
								$codingRegions .= $1;
							}
						}
						my $frame = $codonStart - 1;
						print $writer join("\t", $gene, $transcriptId, $proteinId, $codingRegions, $frame), "\n";
					}
				}
				$string = '';
			}
		}
		close($reader);
	}
	close($writer);
}
{
	my $proteinFeatureFile = "$dataPath/protein_feature.txt";
	open(my $writer, "| sort -u > $proteinFeatureFile");
	foreach my $genomeAssembly (@genomeAssemblyList) {
		my $string = '';
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_protein.gpff.gz |");
		while(my $line = <$reader>) {
			$string .= $line;
			if($line eq "//\n") {
				my $seqInput_object = Bio::SeqIO->new(-string => $string, -format => 'genbank');
				while(my $seq_object = $seqInput_object->next_seq()) {
					my $proteinId = join('.', $seq_object->display_id, $seq_object->version);
					foreach my $feat_object (grep {$_->primary_tag eq 'Region' && $_->has_tag('region_name')} $seq_object->get_SeqFeatures) {
						my ($name) = $feat_object->get_tag_values('region_name');
						my ($start, $end) = ($feat_object->location->start, $feat_object->location->end);
						print $writer join("\t", $proteinId, $start, $end, $name), "\n";
					}
				}
				$string = '';
			}
		}
		close($reader);
	}
	foreach my $gpffFile (@gpffFileList) {
		my $string = '';
		open(my $reader, ($gpffFile =~ /\.gz$/ ? "gzip -dc $gpffFile |" : $gpffFile));
		while(my $line = <$reader>) {
			$string .= $line;
			if($line eq "//\n") {
				my $seqInput_object = Bio::SeqIO->new(-string => $string, -format => 'genbank');
				while(my $seq_object = $seqInput_object->next_seq()) {
					my $proteinId = join('.', $seq_object->display_id, $seq_object->version);
					foreach my $feat_object (grep {$_->primary_tag eq 'Region' && $_->has_tag('region_name')} $seq_object->get_SeqFeatures) {
						my ($name) = $feat_object->get_tag_values('region_name');
						my ($start, $end) = ($feat_object->location->start, $feat_object->location->end);
						print $writer join("\t", $proteinId, $start, $end, $name), "\n";
					}
				}
				$string = '';
			}
		}
		close($reader);
	}
	close($writer);
}
{
	my $pid = open2(my $reader, my $writer, "LC_ALL=C sort -t '\t' -k1,1 -k2,2 -k3,3 -k4,4n -k5,5n | uniq");
	foreach my $genomeAssembly (@genomeAssemblyList) {
		open(my $reader, "gzip -dc $dataPath/$genomeAssembly/$genomeAssembly\_genomic.gtf.gz |");
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			my %tokenHash = ();
			@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
			$tokenHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
			print $writer join("\t", @tokenHash{'transcript_id', 'chromosome', 'strand', 'start', 'end'}), "\n" if($tokenHash{'feature'} eq 'exon');
		}
		close($reader);
	}
	close($writer);
	{
		open(my $writer, "> $dataPath/splice_junction.txt");
		my @tokenListList = ();
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
			if(@tokenListList && grep {$tokenList[$_] ne $tokenListList[0]->[$_]} (0, 1, 2)) {
				my ($transcriptId, $chromosome, $strand) = @{$tokenListList[0]}[0, 1, 2];
				$transcriptId =~ s/_[0-9]+$//;
				if($strand eq '+') {
					my $position = 0;
					foreach my $index (0 .. $#tokenListList) {
						if($index > 0) {
							my $spliceJunction1 = "$chromosome:$tokenListList[$index - 1]->[4]";
							my $spliceJunction2 = "$chromosome:$tokenListList[$index]->[3]";
							print $writer join("\t", $transcriptId, $position, $spliceJunction1, $spliceJunction2), "\n";
						}
						$position += $tokenListList[$index]->[4] - $tokenListList[$index]->[3] + 1;
					}
				}
				if($strand eq '-') {
					my $position = 0;
					foreach my $index (reverse(0 .. $#tokenListList)) {
						if($index < $#tokenListList) {
							my $spliceJunction1 = "$chromosome:$tokenListList[$index + 1]->[3]";
							my $spliceJunction2 = "$chromosome:$tokenListList[$index]->[4]";
							print $writer join("\t", $transcriptId, $position, $spliceJunction1, $spliceJunction2), "\n";
						}
						$position += $tokenListList[$index]->[4] - $tokenListList[$index]->[3] + 1;
					}
				}
				@tokenListList = ();
			}
			push(@tokenListList, \@tokenList);
		}
		if(@tokenListList) {
			my ($transcriptId, $chromosome, $strand) = @{$tokenListList[0]}[0, 1, 2];
			$transcriptId =~ s/_[0-9]+$//;
			if($strand eq '+') {
				my $position = 0;
				foreach my $index (0 .. $#tokenListList) {
					if($index > 0) {
						my $spliceJunction1 = "$chromosome:$tokenListList[$index - 1]->[4]";
						my $spliceJunction2 = "$chromosome:$tokenListList[$index]->[3]";
						print $writer join("\t", $transcriptId, $position, $spliceJunction1, $spliceJunction2), "\n";
					}
					$position += $tokenListList[$index]->[4] - $tokenListList[$index]->[3] + 1;
				}
			}
			if($strand eq '-') {
				my $position = 0;
				foreach my $index (reverse(0 .. $#tokenListList)) {
					if($index < $#tokenListList) {
						my $spliceJunction1 = "$chromosome:$tokenListList[$index + 1]->[3]";
						my $spliceJunction2 = "$chromosome:$tokenListList[$index]->[4]";
						print $writer join("\t", $transcriptId, $position, $spliceJunction1, $spliceJunction2), "\n";
					}
					$position += $tokenListList[$index]->[4] - $tokenListList[$index]->[3] + 1;
				}
			}
		}
		close($writer);
	}
	close($reader);
	waitpid($pid, 0);
}
