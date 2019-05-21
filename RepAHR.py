#!/usr/bin/env python
import sys
import os
import time
import argparse
from subprocess import *
READS_DEPTH="auto"
AVE_READS_LEN=0
THREAD="1"
MAX_MEMORY="250"
OUTPUT_PATH="./temp/"
THRESHOLD="84"
#FILTER_RATE="1"
JELLYFISH_PATH=""
SPADES_PATH=""
JAVA_HOME=""
MIN_CONTIG_LENGTH="100"

is_paired=False
is_single=False
single_file_list=[]
paired_file_list=[]

def introduction():
    print 'Usage: python {0} -c configureFile -r reads_listFile(fastq or fasta file)\n'.format(sys.argv[0]),
    #print 'An example of running the whole pipeline: python main.py -c config.txt -r reads_listfile.txt'

def parse_args():
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
        description='Run the pipeline of zhang')
    # Add arguments
    parser.add_argument(
        '-r', '--reads', type=str, help='Reads listfile name', required=True)
    parser.add_argument(
        '-c', '--config', type=str, help='Configuration file name', required=True)

    args = parser.parse_args()
    configFile=args.config
    listFile= args.reads
    return configFile, listFile

def main_process(configFile, listFile):
    print '*******************************************************************************************'
    print '2018 Jianxin Wang(jxwang@mail.csu.edu.cn), Xiankai Zhang(xkzhang@csu.edu.cn)'
    print 'School of Computer Science and Engineering'
    print 'Central South University'
    print 'ChangSha'
    print 'CHINA, 410083'
    print '*******************************************************************************************'
    print ''
    print 'Step-1: Enter data and parameters'
    global READS_DEPTH, THREAD, OUTPUT_PATH, THRESHOLD, JELLYFISH_PATH, SPADES_PATH, JAVA_HOME, MAX_MEMORY
    #global FILTER_RATE
    global single_file_list, paired_file_list
    global is_single,is_paired
    assert os.path.exists(configFile), "Configuration file does not exist."
    assert os.path.exists(listFile), "Reads listFile does not exist."
    read_configfile(configFile)
    read_inputfile_list(listFile)
    check_library_exist()
    print ''
    print 'Step-2: Determine reads depth'
    depth = estimate_reads_coverage()
    THRESHOLD = depth * 2.3
    THRESHOLD = str(int(THRESHOLD))
    print 'Use kmer threshold:' + THRESHOLD
    print ''
    print 'Step-3: Generate unique kmer use jellyfish'
    generate_kmer_set()
    print ''
    print 'Step-4: Generate frequent kmer set'
    generate_frequent_kmer_set()
    print ''
    print 'Step-5: Filter reads with frequent kmer'
    filter_reads_with_kmer()
    print ''
    print 'Step-6: Generate repeat contigs'
    generate_repeat_spades()
    print 'Step-7: Filter repeat libraries'
    filter_repeat_libraries()
    print 'Pipeline finished.'
###### read in parameters from configuration file. #######
def read_configfile(configFile):
    configf=open(configFile,"r")
    lines=configf.readlines()
    for line in lines:
        parts=line.strip('\n').split("=")
        if len(parts) < 2 or len(parts[0]) == 0 or len(parts[1]) == 0:
            error = 'Configuration file has errors. Maybe some parameters are empty or there are some blank lines.'
            sys.exit(error)
        if parts[0] == 'READ_DEPTH':
            if parts[1] != 'auto':
                global READS_DEPTH
                try:
                    READS_DEPTH=float(parts[1])
                except ValueError:
                    error = "Value READS_DEPTH should be an number,such like 40.0"
                    sys.exit(error)
        elif parts[0] == 'AVE_READS_LEN':
            global AVE_READS_LEN
            AVE_READS_LEN=float(parts[1])
        elif parts[0] == 'THREAD':
            global THREAD
            THREAD=parts[1]
        elif parts[0] == 'OUTPUT_PATH':
            global OUTPUT_PATH
            OUTPUT_PATH=parts[1]
            if OUTPUT_PATH[-1] != "/":
                OUTPUT_PATH = OUTPUT_PATH + "/"
        elif parts[0] == 'THRESHOLD':
            global THRESHOLD
            THRESHOLD=parts[1]
        elif parts[0] == 'MIN_CONTIG_LENGTH':
            global MIN_CONTIG_LENGTH
            MIN_CONTIG_LENGTH=parts[1]
        # elif parts[0] == 'FILTER_RATE':
        #     global FILTER_RATE
        #     FILTER_RATE = parts[1]
        elif parts[0] == 'JELLYFISH_PATH':
            if parts[1] != "GLOBAL":
                global JELLYFISH_PATH
                JELLYFISH_PATH = parts[1]
        elif parts[0] == 'SPADES_PATH':
            if parts[1] != "GLOBAL":
                global SPADES_PATH
                SPADES_PATH = parts[1]
        elif parts[0] == 'JAVA_HOME':
            if parts[1] != "GLOBAL":
                global JAVA_HOME
                JAVA_HOME = parts[1]
        elif parts[0] == 'MAX_MEMORY':
            if parts[1] != "default":
                global MAX_MEMORY
                MAX_MEMORY = parts[1]
    configf.close()



###### read in file names from list file. #######
def read_inputfile_list(listFile):
    global is_paired
    global is_single
    file_reader=open(listFile,'r')
    lines=file_reader.readlines()
    n = 0
    while n < len(lines):
        parts=lines[n].strip('\n').split()
        if len(parts) > 0 and len(parts[0]) > 0 and parts[0] == "#":
            n += 1
            continue
        elif len(parts) <= 0:
            n += 1
            continue
        if len(parts)!=2:
            error = "Input list file format is wrong." \
                    "Each line should be like these: \n" \
                    "1 paired-end-lib1-left-filePath\n" \
                    "1 paired-end-lib1-left-filePath\n" \
                    "-1 single-end-lib1-filePath\n" \
                    "Left reads and right reads from one Paired-end library should have same lib-number,and set single-end file lib-number to -1."
            sys.exit(error)

        if int(parts[0])==-1:
            if is_paired:
                error = "Input file can not contain paired-end reads and single-end reads at same time."
                sys.exit(error)
            single_file_list.append(parts[1])
            is_single=True
            if len(single_file_list) > 9:
                error = "The number of single library is up to 9."
                sys.exit(error)
        else:
            if is_single:
                error = "Input file can not contain paired-end reads and single-end reads at same time."
                sys.exit(error)
            try:
                left_id=int(parts[0])
            except ValueError:
                error = "Library number should be an integer."
                sys.exit(error)
            paired_left=parts[1]
            n+=1
            if n >=len(lines):
                error = "Paired-end library should have two files."
                sys.exit(error)
            parts2=lines[n].strip('\n').split()
            try:
                right_id = int(parts2[0])
            except ValueError:
                error = "Library number should be an integer."
                sys.exit(error)
            paired_right=parts2[1]
            if left_id != right_id:
                error = "Left read and right read from one library should be same."
                sys.exit(error)
            paired_file_list.append([paired_left, paired_right])
            if len(paired_file_list) > 9:
                error = "The number of paired library is up to 9."
                sys.exit(error)
            is_paired=True
        n+=1
    file_reader.close()

##check all the file exist
def check_library_exist():
    global is_paired, is_single
    global single_file_list,paired_file_list
    if is_single:
        for signle_file in single_file_list:
            if not os.path.exists(signle_file):
                error = "{0} does not exist.Please check your reads list file.".format(signle_file)
                sys.exit(error)
    elif is_paired:
        for paired_file in paired_file_list:
            if not os.path.exists(paired_file[0]):
                error = "{0} does not exist.Please check your reads list file.".format(paired_file[0])
                sys.exit(error)
            if not os.path.exists(paired_file[1]):
                error = "{0} does not exist.Please check your reads list file.".format(paired_file[1])
                sys.exit(error)

def estimate_reads_coverage():
    global READS_DEPTH, JELLYFISH_PATH, THREAD, OUTPUT_PATH
    global is_single, is_paired, single_file_list, paired_file_list
    if READS_DEPTH != 'auto':
        try:
            READS_DEPTH = float(READS_DEPTH)
        except ValueError:
            error = "READS_DEPTH should be a positive number."
            sys.exit(error)
        if READS_DEPTH <= 0:
            error = "READS_DEPTH should be a positive number."
            sys.exit(error)
        print "reads coverage was set by user as " + str(READS_DEPTH)
        return READS_DEPTH
    else:
        check_output_path()
        cmd = "cat "
        if is_single:
            for single_file in single_file_list:
                cmd = cmd + single_file + " "
        elif is_paired:
            for paired_file in paired_file_list:
                cmd = cmd + paired_file[0] + " " + paired_file[1] + " "

        if JELLYFISH_PATH == "":
            cmd = cmd + "| jellyfish count -s 100000000 -C -m 15 -t " + THREAD + " -o " + OUTPUT_PATH + "jf_zhang_15.db /dev/stdin"
        else:
            cmd = cmd + "| " + JELLYFISH_PATH +" count -s 100000000 -C -m 15 -t " + THREAD + " -o " + OUTPUT_PATH + "jf_zhang_15.db /dev/stdin"
        print "Running command: \n" + cmd
        Popen(cmd, shell=True, stdout=PIPE).communicate()

        if JELLYFISH_PATH == "":
            cmd2 = "jellyfish histo " +  OUTPUT_PATH + "jf_zhang_15.db > " + OUTPUT_PATH +  "jf_zhang_15.histo"
        else:
            cmd2 = JELLYFISH_PATH + " histo " +  OUTPUT_PATH + "jf_zhang_15.db > " + OUTPUT_PATH +  "jf_zhang_15.histo"
        print "Running command: \n" + cmd2
        Popen(cmd2, shell=True, stdout=PIPE).communicate()
        if not os.path.exists(OUTPUT_PATH + "jf_zhang_15.histo"):
            error = "Jellyfish did not properly finish."
            sys.exit(error)

        depth = cal_coverage_from_histo()
        if os.path.exists(OUTPUT_PATH + "jf_zhang_15.db"):
            os.remove(OUTPUT_PATH + "jf_zhang_15.db")
        return depth


def cal_coverage_from_histo():
    global AVE_READS_LEN, OUTPUT_PATH,JAVA_HOME
    if JAVA_HOME == "":
        cmd = "java"
    else:
        cmd = JAVA_HOME
    cmd = cmd + " -classpath ./bin/ FindThresholdFromHisto " + OUTPUT_PATH + \
              "jf_zhang_15.histo " + str(AVE_READS_LEN) + " " + OUTPUT_PATH + "AutoDetectionThreshold.config"
    print "Running command: \n" + cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    if not os.path.exists(OUTPUT_PATH + "AutoDetectionThreshold.config"):
        error = "Can not detect your coverage from reads.Please set an estimated coverage."
        sys.exit(error)
    threshold_file = open(OUTPUT_PATH + "AutoDetectionThreshold.config", "r")
    threshold = (int)(threshold_file.readlines()[0])
    print "Reach a preliminary threshold : " + str(threshold)
    # #last version
    # histo_file = open(OUTPUT_PATH + "jf_zhang_15.histo", "r")
    # i = 0
    # peak = -1
    # ascending = False
    # for eachLine in histo_file:
    #     parts = eachLine.strip('\n').split()
    #     if i == 5:
    #         cur = int(parts[1])
    #     if i > 5:
    #         last = cur
    #         cur = int(parts[1])
    #         if ascending:
    #             if cur < last:
    #                 break
    #             else:
    #                 peak = int(parts[0])
    #         elif cur > last:
    #             ascending = True
    #     if i >= 500:
    #         break
    #     i += 1
    # if peak == -1 or ascending == False:
    #     error = "Can not detect your coverage from reads.Please set an estimated coverage."
    #     sys.exit(error)
    # histo_file.close()
    # depth = (peak * AVE_READS_LEN)/(AVE_READS_LEN - 15 + 1)
    # depth = round(depth)
    return threshold

def generate_kmer_set():
    global JELLYFISH_PATH, THREAD, OUTPUT_PATH
    global is_single, is_paired, single_file_list, paired_file_list
    check_output_path()
    cmd = "cat "
    if is_single:
        for single_file in single_file_list:
            cmd = cmd + single_file + " "
    elif is_paired:
        for paired_file in paired_file_list:
            cmd = cmd + paired_file[0] + " " + paired_file[1] + " "
    if JELLYFISH_PATH == "":
        cmd = cmd + "| jellyfish count -s 100000000 -C -m 31 -t " + THREAD + " -o " + OUTPUT_PATH + "jf_zhang.db /dev/stdin"
    else:
        cmd = cmd + "| " + JELLYFISH_PATH + " count -s 100000000 -C -m 31 -t " + THREAD + " -o " + OUTPUT_PATH + "jf_zhang_15.db /dev/stdin"
    print "Running command: \n" + cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()

    if JELLYFISH_PATH == "":
        cmd3 = "jellyfish histo " + OUTPUT_PATH + "jf_zhang.db > " + OUTPUT_PATH + "jf_zhang_kmers.histo"
	cmd2 = "jellyfish dump " + OUTPUT_PATH + "jf_zhang.db > " + OUTPUT_PATH + "jf_zhang_kmers.fasta"
    else:
	cmd3 = JELLYFISH_PATH + " histo " + OUTPUT_PATH + "jf_zhang.db > " + OUTPUT_PATH + "jf_zhang_kmers.histo"
        cmd2 = JELLYFISH_PATH + " dump " + OUTPUT_PATH + "jf_zhang.db > " + OUTPUT_PATH + "jf_zhang_kmers.fasta"
    print "Running command: \n" + cmd2
    Popen(cmd3, shell=True, stdout=PIPE).communicate()	
    Popen(cmd2, shell=True, stdout=PIPE).communicate()
    if not os.path.exists(OUTPUT_PATH + "jf_zhang_kmers.fasta"):
        error = "Jellyfish did not properly finish."
        sys.exit(error)
    #if os.path.exists(OUTPUT_PATH + "jf_zhang.db"):
    #    os.remove(OUTPUT_PATH + "jf_zhang.db")

def generate_frequent_kmer_set():
    global THRESHOLD, OUTPUT_PATH, JAVA_HOME
    check_output_path()
    if JAVA_HOME == "":
        cmd = "java"
    else:
        cmd = JAVA_HOME
    cmd = cmd + " -classpath ./bin/ GetFreKmerWithThreshold " + OUTPUT_PATH + \
              "jf_zhang_kmers.fasta " + THRESHOLD + " " + OUTPUT_PATH + "frequent_kmers_" + THRESHOLD + ".fasta"
    print "Running command: \n" + cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    if not os.path.exists(OUTPUT_PATH + "frequent_kmers_" + THRESHOLD + ".fasta"):
        error = "Generate frequent kmer set failed."
        sys.exit(error)

def filter_reads_with_kmer():
    global OUTPUT_PATH, JAVA_HOME, is_single, is_paired, THRESHOLD
    global paired_file_list,single_file_list
    process_list = []
    if JAVA_HOME == "":
        cmd = "java"
    else:
        cmd = JAVA_HOME
    if is_paired:
        n = 1
        for paired_file in paired_file_list:
            cmd_cur = cmd + " -classpath ./bin/ FilterSequence P fastq " + paired_file[0] + " " + paired_file[1] + " " + OUTPUT_PATH + "filtered_" \
                  + THRESHOLD + "_lib" + str(n) + "_1.fastq " + OUTPUT_PATH + "filtered_" + THRESHOLD + "_lib" + str(n) + "_2.fastq " \
                  + OUTPUT_PATH + "frequent_kmers_" + THRESHOLD + ".fasta 0.5 both"
            print "Running command: \n" + cmd_cur
            Popen(cmd_cur, shell=True, stdout=PIPE).communicate()
            # process_list.append()
            n += 1
        print('1')
        # for p in process_list:
        #     print(p)
        #     while p.poll() is None:
        #         time.sleep(1)
        # print('2')
        print 'All processes done.'
    elif is_single:
        n = 1
        for single_file in single_file_list:
            cmd_cur = cmd + " -classpath ./bin/ FilterSequence S fastq " + single_file + " " + OUTPUT_PATH + "filtered_" \
                  + THRESHOLD + "_lib" + str(n) + ".fastq " + OUTPUT_PATH + "frequent_kmers_" + THRESHOLD + ".fasta 0.5 both"
            print "Running command: \n" + cmd_cur
            Popen(cmd_cur, shell=True, stdout=PIPE).communicate()
            # process_list.append()
            n += 1
        # for p in process_list:
        #     while p.poll() is None:
        #         time.sleep(1)
        print 'All processes done.'

def generate_repeat_spades():
    global OUTPUT_PATH, SPADES_PATH, THRESHOLD, MAX_MEMORY
    global is_single, is_paired, single_file_list, paired_file_list
    if SPADES_PATH == "":
        cmd = "spades.py"
    else:
        cmd = SPADES_PATH
    if is_paired:
        n = 1
        while n <= len(paired_file_list):
            cmd = cmd + " --pe" + str(n) + "-1 " + OUTPUT_PATH + "filtered_"  + THRESHOLD + "_lib" + str(n) + "_1.fastq " \
                  + " --pe" + str(n) + "-2 "+ OUTPUT_PATH + "filtered_"  + THRESHOLD + "_lib" + str(n) + "_2.fastq "
            n += 1
        cmd = cmd + "-m " + MAX_MEMORY + " -t " + THREAD + " --cov-cutoff auto -o " + OUTPUT_PATH + "spades_repeat_lib"
        print "Running command: \n" + cmd
        Popen(cmd, shell=True, stdout=PIPE).communicate()
    elif is_single:
        n = 1
        while n <= len(single_file_list):
            cmd = cmd + " --s" + str(n) + "-1 " + + OUTPUT_PATH + "filtered_" + THRESHOLD + "_lib" + str(n) + ".fastq " \
                  + OUTPUT_PATH + "filtered_"  + THRESHOLD + "_lib" + str(n) + ".fastq "
            n += 1
        cmd = cmd + "-m " + MAX_MEMORY + " -t " + THREAD + " --cov-cutoff auto -o " + OUTPUT_PATH + "spades_repeat_lib"
        print "Running command: \n" + cmd
        Popen(cmd, shell=True, stdout=PIPE).communicate()

def filter_repeat_libraries():
    global OUTPUT_PATH, JAVA_HOME, THRESHOLD, MIN_CONTIG_LENGTH
    check_output_path()
    if JAVA_HOME == "":
        cmd = "java"
    else:
        cmd = JAVA_HOME
    COV_THRESHOLD=(int)(THRESHOLD)/2
    cmd = cmd + " -classpath ./bin/ FilterContigs " + OUTPUT_PATH + \
          "spades_repeat_lib/contigs.fasta " + OUTPUT_PATH + "final_repeat_lib.fa " + MIN_CONTIG_LENGTH + " " + str(COV_THRESHOLD)
    print "Running command: \n" + cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    if not os.path.exists(OUTPUT_PATH + "final_repeat_lib.fa"):
        error = "Filter repeat libraries failed."
        sys.exit(error)

def check_output_path():
    global OUTPUT_PATH
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        introduction()
        raise SystemExit

    configFile, listFile = parse_args()
    main_process(configFile, listFile)
