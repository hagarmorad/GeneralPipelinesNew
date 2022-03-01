# GeneralPipelinesNew
## Usage:
 `python3 main.py [-h] | [-r REFSEQ_PATH] [-i FASTQ_PATH] [options] `
#### -i /input/path : provide input path  (required)
`/input/path` - the path to fastq files.

#### -r|--refseq <refseq/path/> : provide refseq path (required)
Don't worry about indexing the fasta file, it happens automatically.

#### --flu (optional)
run the pipeline for influenza viruses. \
The reference should contain the virus segments seperated by fasta headers. \
 `python3 main.py  -r path/to/influenza_ref.fa -i path/to/fastq/location --flu `
 
#### --polio (optional)
run the pipeline for polio virus. \
The reference should contain the 3 Sabins seperated by fasta headers. \
 `python3 main.py  -r path/to/polio_ref.fa -i path/to/fastq/location --polio `
 
#### --de_novo (optional)
run de-novo analysis in addition to the regular pipeline . \
The output will be seperated to contig based and fastq based analysis. \
 `python3 main.py  -r path/to/ref.fa -i path/to/fastq/location --de_novo `

## Main Steps:

# General pipeline:

* Index the reference sequence. \
`bwa index /refseq/path`

* For each fastq sample - map to reference and keep as bam file. \
For each fastq file run: \
    `bwa mem /refseq/path sample_R1 sample_R2 | samtools view -b - > BAM/sample_name.bam` \
R1 and R2: forward and reverse paired ends.
This command execute in parallel for 5 fastq files.

* Keep only mapped reads of each sample. \
For each bam file run: \
    `samtools view -b -F 260 file_name.bam > BAM/file_name.mapped.bam` 

* Sort and index each sample in BAM. \
For each mapped bam file run:
    `samtools sort BAM/file_name.mapped.bam BAM/file_name.mapped.sorted.bam`
    
* Create consensus files for each sample. \
For each mapped and sorted bam file run:
    `samtools mpileup -A  BAM/file_name.mapped.sorted.bam | ivar consensus -m 5 -p CNS5/file_name`

* Align all consensus sequences and references sequence: \
Gather all fasta sequences (consensus + reference): \
    `cat CNS5/*.fasta /refseq/path > alignment/all_not_aligned.fasta` \
Align with augur (mafft based) and save output in alignment/ directory: \
  ` augur align \
  --sequences alignment/all_not_aligned.fasta \
  --reference-sequence /refseq/path \
  --output alignment/all_aligned.fasta`

* report.csv: \
    some coverage statistics and number of mapped reads: `samtools coverage -H BAM/sample_name.mapped.sorted.bam`  \
    total number of reads in sample: `samtools view -c BAM/sample_name.bam` \
    depth of each position in sample: `samtools depth BAM/sample_name.mapped.sorted.bam`  \
Check out pipeline.sh code for specifics.

- parallel commands are executing using python multiprocessing package.

# Influenza addional steps:

* Splitting the BAM files by the references they were mapped to to multiple files. further analysis is made using the splitted files.
for each bam file run:
`bamtools split -in bam_file -reference`

* concatenate the consensus files by sample.

* concatenate the consensus files by reference segment.

* multiple sequnce alignment:
concatenate the refence segment into one sequence with no headers.
concatenate the by-sample consensus reference to one multi-fasta.
run general pipeline augur align function.



