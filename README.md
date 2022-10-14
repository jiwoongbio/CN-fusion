## Requirements

1. Perl - https://www.perl.org

2. BioPerl - https://bioperl.org
   * Bio::DB::Fasta
   * Bio::SeqIO

3. Burrows-Wheeler Aligner (BWA) - https://bio-bwa.sourceforge.net

4. SAMtools - https://www.htslib.org
   * bgzip
   * tabix

5. Basic linux commands: bash, rm, gzip, sort, echo, find

6. lftp - https://lftp.yar.ru


## Usage example (human genome GRCh38)

1. Prepare reference data

```
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot refseq/H_sapiens/mRNA_Prot'
perl CN-fusion_data.pl `for file in refseq/H_sapiens/mRNA_Prot/*.gpff.gz; do echo "-p $file"; done` ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13
```

2. Align sequencing reads

```
bwa mem -t $THREADS -T 19 -h 200 -Y data/rna_genomic.fasta $SAMPLE.1.fastq.gz $SAMPLE.2.fastq.gz | gzip > $SAMPLE.rna_genomic.sam.gz
```

3. Extract fusion reads

```
perl fusion_read.pl -p $THREADS $SAMPLE.rna_genomic.sam.gz | gzip > $SAMPLE.rna_genomic.fusion_read.txt.gz
```

4. Extract read mapping regions and build index

```
perl region_read.pl -p $THREADS $SAMPLErna_genomic.sam.gz | bgzip > $SAMPLE.rna_genomic.region_read.txt.gz
tabix -s 1 -b 2 -e 3 $SAMPLE.rna_genomic.region_read.txt.gz
```

5. Identify coding/noncoding fusions

```
time perl CN-fusion.pl -p $THREADS -r $SAMPLE.rna_genomic.region_read.txt.gz $SAMPLE.rna_genomic.fusion_read.txt.gz | gzip > $SAMPLE.CN-fusion.txt.gz
```
