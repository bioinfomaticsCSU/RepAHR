# RepAHR
RepAHR is used to identify repeats(repetitive sequences) in genome using Next-Generation Sequencing reads.
## Installation
The OSs must be Linux.
### Dependencies
Before the installation of RepAHR, the following tools or environment are required:
* [jre.1.8.0](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
* [python 2.7.13  or higher](https://www.python.org/downloads/release/python-2713/)
* [SPAdes.3.12.0](http://cab.spbu.ru/software/spades/)
* [Jellyfish 2.0.0](https://github.com/gmarcais/Jellyfish/tree/master/include/jellyfish)

To install RepAHR, follow the instruction below:
```
git clone https://github.com/bioinfomaticsCSU/RepAHR
```
## Edit configure files:
Before running RepAHR, need to edit 2 configure files: **parameter_config_file** and **reads_config_file**.
Some parameters need to be set in the parameter_config_file:
    
    READS_DEPTH=auto
    AVE_READS_LEN=101
    THREAD=40
    MAX_MEMORY=default
    OUTPUT_PATH=/home/zhangxk/temp
    FILTER_RATE=1
    JELLYFISH_PATH=GLOBAL
    SPADES_PATH=GLOBAL
    JAVA_HOME=GLOBAL
    MIN_CONTIG_LENGTH=100

    Explanation of parameter_config_file:

    * 'READS_DEPTH': the average depth of input NGS reads. If the user know the value of the average
    depth, fill in the rounded number after the equal sign. If the user does not know the specific
    value, fill in the 'auto' after the the equal sign, RepAHR will estimate the average depth of the
    input NGS reads, but there may be some errors.
    * 'AVE_READS_LEN': the average length of all the input NGS reads. User need to fill in the rounded
    number the equal sign.
    * 'THREAD': the running threads when running RepAHR.
    * 'MAX_MEMORY': the maximum memory the user wants to use, in GB. Filling in 'default' means RepAHR
    is allowed to use all available physical memory of the machine.
    * 'FILTER_RATE': a parameter used to obtain high-frequency reads, usually no need to modify.
    * 'OUTPUT_PATH': the path of RepAHR to give the final result and some intermediate files.
    * 'JELLYFISH_PATH': if Jellyfish is in the system environment, fill in 'GLOBAL' after the equal
    sign. If not, fill in the absolute path of Jellyfish after the equal sign.
    * 'SPADES_PATH': if SPADES is in the system environment, fill in 'GLOBAL' after the equal sign. If
    not, fill in the absolute path of SPADES after the equal sign.
    * 'JAVA_HOME': if Java is in the system environment, fill in 'GLOBAL' after the equal sign. If not,
    fill in the absolute path of JAVA_HOME after the equal sign.
    * 'MIN_CONTIG_LENGTH': the repeats which is shorter than this value will be filtered out from the
    fianl result.

The location of ecah reads file need to be given in the reads_config_file:

    1 /home/zhangxk/data/file1_1.fastq
    1 /home/zhangxk/data/file2_2.fastq
    2 /home/zhangxk/data/file2_1.fastq
    2 /home/zhangxk/data/file2_2.fastq

    Explanation of reads_config_file:

    Each line include the infromation of one reads file, the beginning of each line is the number of the
    reads file, followed by a space, and then is the absolute path of the reads file. The number of the
    two files of one paired-end reads libraray should be same in the reads_config_file. Note that 
    there cannot be both single-end reads and paired-end reads, can only be on of them.

## Running
```
python RepAHR.py -c parameter_config_file -r reads_config_file
```
'RepAHR.py' 'parameter_config_file' 'reads_config_file' are all the absolute or relative path of the 3 files.
## Result
The final result will be named final_repeat_lib.fa in the folder given in the parameter_config_file. And the raw contig file is name contig.fasta in the subfolder './spades_repeat_lib'.
