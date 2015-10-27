#fastq\_qual\_trimmer
##What it does?
It trims and filters Sanger encoded FASTQ files by quality, length and
homopolymer count.
##Installation
To install on Unix systems with git, type:
```
git clone https://github.com/ivars-silamikelis/fastq_qual_trimmer
cd fastq_qual_trimmer
make
```
binary file will appear that can be copied to `/usr/local/bin/` 
or added to your `$PATH`

##Usage
If you wish to trim by quality 25 both ends of reads in FASTQ file 
named `test.fq`, type
```
./fastq_qual_trimmer -i test.fq -q 25 
```

##Available options
Option | Argument | Description
  ---  |    ---   |    ---
**-h** | None     | Prints help message
**-i** | Filename | Specifies input file
**-q** | Integer  | Smallest allowed quality for nucleotides at read ends (default: 20)
**-m** | Integer  | Smallest allowed mean read quality (default 0)
**-l** | Integer  | Minimal allowed read length (default 1)
**-H** | Integer  | Reads with homopolymers with specified length or higher will be filtered. Use 0 to switch off (default: 0)
**-w** | None     | Use sliding window trimming. Trim nucleotides from both ends by calculating mean quality of nucleotides in the window (default: off)
**-s** | Integer  | Specify size of sliding window (default: 5)
**-b** | None     | Use sliding window and after that - default trimming approaches (default: off)


