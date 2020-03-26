# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
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
	system("diamond makedb --in $proteinFastaFile --db $proteinFastaFile.dmnd");
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
