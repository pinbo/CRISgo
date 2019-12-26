# CRISgo

The Golang implementation of CRISpyJZ (https://github.com/pinbo/CRISpyJZ).

CRISgo is a tool to process NGS data of CRISPR editings. I Just rewrote CRISpyJZ with golang to speed it up, because I found it kind of slow when you have more than 100 fastq files to process.

Just like CRISpyJZ, CRISgo searches two flanking sequences of gRNA in fastq files (plan text files), and check whether the gRNA was edited and compare the length with the wild type sequence to see whether there is an indel.

To understand how it works, check out the **README** of [CRISpyJZ](https://github.com/pinbo/CRISpyJZ) and [CRIS.py](https://github.com/patrickc01/CRIS.py), and read the original CRIS.py paper, [CRIS.py: A Versatile and High-throughput Analysis Program for CRISPR-based Genome Editing](https://www.nature.com/articles/s41598-019-40896-w).

## Differences from CRISpyJZ

1. It is written in golang, and the compiled version is about 15 to 20 times faster than CRISpyJZ.

2. The strand option was removed. So you need give all 3 sequences in the arguments from the same strand of the reference file.

## Differences among CRISgov2, CRISgov3, and CRISgov4

- **CRISgov2** and **CRISgov3** process ".fastq" files, which are uncompressed. You need to uncompress the ".fastq.gz" files first. In window system, you can use 7-Zip to unzip all the .gz files.

- **CRISgov4** processes directly on "fastq.gz" files, so no need to uncompress first.

- **CRISgov2** searches the input sequences in all the "fastq" files first, then searches using the reverse complements of the input sequences to search on those fastq files that have 0 matches in the first round search. So it is "orientation" insenstive.

- **CRISgov3** and **CRISgov4** suppose your input sequences on the 'xxxx_R1_001.fastq.gz' files, so first search on these files, then use the reverse complements of your input sequences to search 'xxxx_R2_001.fastq.gz' files.

## Usage

```sh
# using compiled binary file
./CRISgo (or CRISgo.exe on Windows) reference-file output_file_prefix gene_name left_flanking_sequence right_flanking_sequence gRNA_sequence

# You can also use 'go run' if you have go installed on your computer
go run CRISgov3.go reference-file output_file_prefix gene_name left_flanking_sequence right_flanking_sequence gRNA_sequence
```

Command line parameters:

1. **reference-file**: a fasta file with all reference sequences you will need.
2. **output_file_prefix**: two files, *output_file_prefix.csv* and *top10_reads_output_file_prefix.txt*, will be created.
3. **gene_name**: the name of the target gene in the reference fasta file.
4. **left_flanking_sequence**: left flanking sequence, a unique ~20bp sequence in your reference. Must be at **5'** of the gRNA sequence. Case insensitive.
5. **right_flanking_sequence**: right flanking sequence, a unique ~20bp sequence in your reference. Must be at **3'** of the gRNA sequence. Case insensitive.
6. **gRNA_sequence** : your gRNA sequence (or its reverse complement sequence it is on the "-" strand). Case insensitive.

## Outputs

Two files will be created:

1. A csv file with the editing efficiency and top 10 indels in each fastq file.

2. A txt file with names "top10_reads_xxxx.txt" that lists the top 10 most abundant reads in each fastq file.

See the "**Examples**" folder for some example outputs. The "**Examples**" folder has paired end fastq files from [CRIS.py](https://github.com/patrickc01/CRIS.py). You can download the proper version of CRISgo and test on your computer. 

## Notes

1. Please make sure the 3 sequences, left_flanking_sequence, right_flanking_sequence, and gRNA_sequence, are on the same strand.

2. `CRISgov3.go` supposes all the 3 sequences are in the R1 fastq files if your NGS are paired end. It will search the 3 sequences in all the R1 fastq files, then use the reverse complment of the 3 sequences to search the R2 fastq files. You can use command like `grep your-sgRNA-sequence *_R1_*.fastq | less` and `grep your-sgRNA-sequence *_R2_*.fastq| less` to make sure your flanking sequences and sgRNA sequences are on the same strand as the R1 fastq files.

3. If you do not understand, you can just use `CRISgov2.go`, which is R1/R2 insensitive.

## Get the binary file

You can go to ["**Releases**"](https://github.com/pinbo/CRISgo/releases) to download the version you need (v2 is flexible than v3, but a little slow).

You can also install golang on your computer, then run the command below to compile it for your own computer.

`go build CRISgov3.go`
