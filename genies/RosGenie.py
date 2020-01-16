#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


def SUM(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def deExt(string):
    ls = string.split(".")
    ls2 = ls[0:len(ls) - 1]
    outstring = "".join(ls2)
    return outstring


def unique(ls, ls2):
    unqlist = []
    for i in ls:
        if i not in unqlist and i in ls2:
            unqlist.append(i)
    return len(unqlist)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def remove2(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    # outString = "".join(emptyList)
    return emptyList


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" # ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" # ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


def main():
    parser = argparse.ArgumentParser(
        prog="ROSGenie.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************
    
        Developed by Arkadiy Garber;
        University of Southern California, Earth Sciences
        University of Montana, Biological Sciences
        Please send comments and inquiries to arkadiyg@usc.edu            
                                                                                                                   
         /////////////,                  .//,                  .***//////////.     
       ###################.         ##############          ###################       
       #####################     ,###################     ,####################         
                       #####                  *#######    #####,                      
         /            ######    /,               ######   ######*                            
       #####        #######   #####               /#####   ##########                     
       #####      #######(    #####                #####.    ############.            
       #####    #######(      #####                #####.        ############             
       #####  ,######/        #####/               #####             /########/                  
       #####  #####/           #####(            .#####(                 *#####             
       #####  #####            .######(        /######(                   #####       
       #####  #############.     ####################      ####################  GENIE
       #####   /############       ################       ###################(              
        ###       (########            /######(            ################                   
    
        STRESS FOR LESS
    
        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/ 
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of ORFs in FASTA format")
    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)")

    parser.add_argument('-out', type=str, help="name output directory (default=rosgenie_out)",
                        default="rosgenie_out")

    parser.add_argument('--d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                              "This number should be an integer and should reflect the maximum number of "
                                              "genes in between putative iron-related genes identified by the HMM database "
                                              "(default=10)", default=10)


    parser.add_argument('--cpu', type=int, help="number of threads to allow for hmmsearch (default = 1)", default=1)

    parser.add_argument('-bams', type=str, help="a tab-delimited file with two columns: first column has the genome or "
                                                "metagenome file names; second column has the corresponding BAM file "
                                                "(provide full path to the BAM file). Use this option if you have genomes "
                                                "that each have different BAM files associated with them. If you have a set "
                                                "of bins from a single metagenome sample and, thus, have only one BAM file, "
                                                " then use the \'-bam\' option. BAM files are only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts.", default="NA")

    parser.add_argument('-bam', type=str, help="BAM file. This option is only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts. If you have more than one BAM file"
                                               "corresponding to different genomes that you are providing, then use the \'-bams\' "
                                               "option to provide a tab-delimited file that denotes which BAM file (or files) belongs "
                                               "with which genome", default="NA")

    parser.add_argument('--makeplots', type=str,
                        help="Would you like LithoGenie to make some figures from your data? y = yes, n = no (default = n). "
                             "If so, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Warning: this part of the program is currently under beta-testing, and if there are any problems running Rscript, "
                             "or installing any of the required packages, you may get a bunch of error messages at the end. "
                             "This will not crash the program, however, and you can still expect to "
                             "see the main output (CSV files) from Genie.", default="n")




    # CHECKING FOR CONDA INSTALL
    os.system("echo ${ros_hmms}/hmm-meta.txt > HMMlib.txt")
    file = open("HMMlib.txt")
    for i in file:
        location = i.rstrip()

    os.system("rm HMMlib.txt")
    try:
        bitscores = open(location)
        conda = 1
    except FileNotFoundError:
        conda = 0

    if conda == 0:
        parser.add_argument('-hmm_dir', type=str,
                            help='directory of HMMs. Provide the directory path (e.g. hmms/iron or hmms/lux) if you do not have conda installed and were not able to run the setup.sh script',
                            default="NA")

        # parser.add_argument('--R', type=str,
        #                     help="location of R scripts directory (note: this optional argument requires Rscript to be "
        #                          "installed on your system). The R scripts directory is in the same directory as the "
        #                          "Genie.py code", default="NA")

    args = parser.parse_args()

    # CHECKING ARGUMENTS AND PATHS
    cwd = os.getcwd()
    print("checking arguments")
    if conda == 0:
        if args.hmm_dir != "NA":
            print(".")
        else:
            print("You have not provided the location of the HMM library via the -hmm_dir argument. Please do so, and try "
                  "again. The HMM library is found within the same directory as the FeGenie executable.")
            print("Exiting")
            raise SystemExit

    if args.bin_dir != "NA":
        binDir = args.bin_dir + "/"
        binDirLS = os.listdir(args.bin_dir)
        print(".")
    else:
        print("Looks like you did not provide a directory of genomes/bins or assemblies.")
        print("Exiting")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print('Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
              ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit

    if conda == 1:
        os.system("echo ${ros_hmms} > HMMlib.txt")
        file = open("HMMlib.txt")
        for i in file:
            HMMdir = (i.rstrip())
            HMMdirLS = os.listdir(HMMdir)
        os.system("rm HMMlib.txt")
    else:
        HMMdir = args.hmm_dir
        HMMdirLS = os.listdir(args.hmm_dir)


    # STARTING MAIN ALGORITHM
    bits = open(HMMdir + "/hmm-meta.txt", "r")
    bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    print("\nreading in HMM bitscore cut-offs...")
    for i in bits:
        ls = i.rstrip().split("\t")
        bitDict[ls[0]]["bit"] = ls[1]
    print("...")


    os.system("mkdir " + args.out)
    count = 0
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    out = open("%s/rosgenie.csv" % (args.out), "w")
    out.write(
        "bin" + "," + "gene" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "sequence" + "\n")
    for i in binDirLS:
        HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.bin_ext:

            fastaFile = open(args.bin_dir + "/" + i, "r")
            fastaFile = fasta(fastaFile)
            os.system("mkdir " + args.bin_dir + "/" + i + "-HMM")
            count = 0
            for hmm in HMMdirLS:
                if hmm != "hmm-meta.txt":
                    count += 1
                    perc = (count / len(HMMdirLS)) * 100
                    sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc + 8))
                    sys.stdout.flush()
                    if not re.match(r'^\.', hmm):
                        os.system("hmmsearch --cpu " + str(args.cpu) +
                                  " --tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
                                  " -o " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                                  HMMdir + "/" + hmm + " " +
                                  binDir + "/" + i)
                        os.system("rm " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt")

            HMMresults = os.listdir(args.bin_dir + "/" + i + "-HMM")
            for result in HMMresults:
                if not re.match(r'^\.', result):
                    file = open(args.bin_dir + "/" + i + "-HMM/" + result)
                    for line in file:
                        line = line.rstrip()
                        if not re.match(r'^#', line):
                            ls = (delim(line))
                            orf = ls[0]
                            gene = ls[2]
                            evalue = ls[4]
                            bitscore = ls[5]
                            hmmFileName = (result.split("__")[1])
                            hmmFileName = hmmFileName.split(".hm")[0]
                            hmmName = deExt(result.split("__")[1])
                            threshold = bitDict[hmmFileName]["bit"]
                            process = bitDict[hmmFileName]["process"]
                            element = bitDict[hmmFileName]["element"]
                            if float(threshold) > 0:
                                if float(bitscore) > float(threshold):
                                    if orf not in HMMdict.keys():
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                        HMMdict[orf]["process"] = process
                                        HMMdict[orf]["element"] = element
                                    else:
                                        if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                            HMMdict[orf]["hmm"] = hmmName
                                            HMMdict[orf]["evalue"] = evalue
                                            HMMdict[orf]["bitscore"] = bitscore
                                            HMMdict[orf]["process"] = process
                                            HMMdict[orf]["element"] = element
                                        else:
                                            pass
                            else:
                                if float(evalue) <= float(1E-20):
                                    if orf not in HMMdict.keys():
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                        HMMdict[orf]["process"] = process
                                        HMMdict[orf]["element"] = element
                                    else:
                                        if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                            HMMdict[orf]["hmm"] = hmmName
                                            HMMdict[orf]["evalue"] = evalue
                                            HMMdict[orf]["bitscore"] = bitscore
                                            HMMdict[orf]["process"] = process
                                            HMMdict[orf]["element"] = element
                                        else:
                                            pass

            for key in HMMdict.keys():
                out.write(
                    i + "," + HMMdict[key]["hmm"] + "," +
                    "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key]["bitscore"] + "," + fastaFile[key]
                    + "\n")
            os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")

    out.close()

    os.system("mv %s/rosgenie.csv %s/rosgenie-summary.csv" % (args.out, args.out))


    # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    print("")
    print("....")
    print(".....")
    # COVERAGE-BASED ABUNDANCE
    if args.bams != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
        BAMmap = open(args.bams)
        for i in BAMmap:
            string = ''
            ls = i.rstrip().split("\t")
            cell = ls[0]
            for j in ls[1:]:
                string += " "
                string += j
            try:
                depth = open("%s/%s.depth" % (args.out, cell))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (args.out, cell, string))
                print("processing... " + cell)
                depth = open("%s/%s.depth" % (args.out, cell))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]

        os.system("mkdir %s/contigDepths" % args.out)
        os.system("mv %s/*depth %s/contigDepths/" % (args.out, args.out))

        cats = ["1-cysPeroxiredoxin_Cterminal", "Alkyl-hydroperoxide-reductase-ThiolSpecificAntioxidant", "catalase_Dyp_perox",
                "catalase_Glutathione_peroxidase", "Catalase", "catalase-peroxidase", "Catalase-rel-immune-response",
                "CCP_MauG", "Sod_CuZn", "SOD_FeMn_Cterminal", "SOD_FeMn_Nterminal"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/rosgenie-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[3]
                    contig = allButTheLast(orf, "_")
                    gene = ls[1]
                    Dict[cell][gene].append(float(depthDict[cell][contig]))

        outHeat = open("%s/rosgenie.readDepth.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            outHeat.write(i + ",")
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i + ",")
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
            outHeat.write("\n")

        outHeat.close()


    # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE
    elif args.bam != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

        try:
            depth = open("%s.depth" % (args.bam))
            for k in depth:
                LS = k.rstrip().split("\t")
                if LS[0] != "contigName":
                    depthDict[LS[0]] = LS[2]

        except FileNotFoundError:
            os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.bam, args.bam))
            depth = open("%s.depth" % (args.bam))
            for k in depth:
                LS = k.rstrip().split("\t")
                if LS[0] != "contigName":
                    depthDict[LS[0]] = LS[2]

        cats = ["1-cysPeroxiredoxin_Cterminal", "Alkyl-hydroperoxide-reductase-ThiolSpecificAntioxidant", "catalase_Dyp_perox",
                "catalase_Glutathione_peroxidase", "Catalase", "catalase-peroxidase", "Catalase-rel-immune-response",
                "CCP_MauG", "Sod_CuZn", "SOD_FeMn_Cterminal", "SOD_FeMn_Nterminal"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/rosgenie-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[3]
                    gene = ls[1]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[contig]))

        outHeat = open("%s/rosgenie.readDepth.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            outHeat.write(i + ",")
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i + ",")
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
            outHeat.write("\n")

        outHeat.close()
        print('......')
        print(".......")
        print("Finished!")


    # GENE COUNTS-BASED ABUNDANCE
    else:
        cats = ["1-cysPeroxiredoxin_Cterminal", "Alkyl-hydroperoxide-reductase-ThiolSpecificAntioxidant", "catalase_Dyp_perox",
                "catalase_Glutathione_peroxidase", "Catalase", "catalase-peroxidase", "Catalase-rel-immune-response",
                "CCP_MauG", "Sod_CuZn", "SOD_FeMn_Cterminal", "SOD_FeMn_Nterminal"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/rosgenie-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            print(ls)
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[3]
                    gene = ls[1]
                    Dict[cell][gene].append(gene)

        print("\n\n")
        normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.bin_dir):
            if lastItem(i.split(".")) == args.bin_ext:
                file = open("%s/%s" % (args.bin_dir, i), "r")
                file = fasta(file)
                print(i)
                normDict[i] = len(file.keys())

        print("\n\n")
        outHeat = open("%s/rosgenie.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            outHeat.write(i + ",")
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i + ",")
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    print(j)
                    print(normDict[j])
                    outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
            outHeat.write("\n")
        outHeat.close()

        print('......')
        print(".......")
        print("Finished!")


if __name__ == '__main__':
    main()

