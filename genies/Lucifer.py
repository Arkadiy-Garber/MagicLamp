#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD BITSCORE CUTOFF FOR EACH HMM HIT


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
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
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
        prog="Lucifer.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************
    
        Developed by Arkadiy Garber and Gustavo Ramírez;
        University of Southern California, Earth Sciences
        Please send comments and inquiries to arkadiyg@usc.edu
    
                                                        @                         
                                                        /@*     @                  
                                                   @     &     @@                  
                                                   @@    %    (%/                  
                                                   &%,   /#   %&                   
                                                   /#    .&   *.                   
                                                    ,   ..    ,                    
                                                   ,*   *,,   *.                   
                           ....                    ,,   ***/  ./                   
                        ...                        .,* ..... /.,(                  
                         ..             .            ,..../*,,                     
                          ...             ..           (.....                      
                           ....             .          /./..,                      
                          ., ...,*,     ......          ,*,.                       
                           ,.*%,,.,.........            ./*.                       
                            .,*.,,.....                  ,*                        
                            ..,. ,                        *.,                      
                           ,,.,.,,                      , ...                      
                          .....*,.*.  ,*,               , ,                        
                   ,.**,,....,,....,,**** ..         /(#* /.                       
                   ,**.,,,,,,.,.,,,*,.....,,,.     #####/./.                       
                 .. **,.*************..,,...,**,***..///*,/                        
               ..,***, ..,****.***,. .., % ,,,,.......   ./.                       
               ,,**.. ..,,*...,,.**,,..                  ./.                       
             ,*. .*   /.,..**,***.,*..                   ./.                       
            .,.,*,,     ,..***.**.,*.                     /.                       
              ,.,.#     .,.,,,.,*.,,                     ./.                       
               ,//(/,.  ,,*,***,.*,.                     ./.                       
               .*//##** *.,,,,,, ,,                      ./.                       
                  /.,,    ,,*,*,,.,                      ./.                       
                   ,,,. . ...... .                       ./.                       
                 ,.., ..  .......... .                   ./.                       
                      . .  ...... ... .                  ./.                       
                    , . . ..  ... ... ..                  /.                       
                      ..... ..  ..........               ./.                       
                   .. .................. . .             ./.                       
                   ........... ........ .... .           ./                        
                   .. .......... ....... ..... .         ./                        
                 ,.. . .......... ... ... .. ... ..       .                        
                ,,.. ................. ... ...........*, .* ,                      
              ,, . .. . ............ ...... ... ... ,.... /                        
             *,    .. .. ................... .... , .... ./                        
             ,#  . ...... .......  ............ , ......,./ ...                    
             *   . ... ... ....... ...... ....,.......... .,,* ..                  
             ,*  . .... ........... . .......... .. ... . / .. ., .                
              .* .......... ..............*....... .. ... / ..... , .              
                .*...... ... ......... *,.............(  */    . ,  , ,            
                 . .,,**  ... .. ,*,,.. ..   .... .....  //.     . ,**,*           
                   ..  ..         /   /  . ,,   ... ....../*        .*,(           
                 .**,.*                  ..,*.,,,. ...               ..            
                 ,.**.@                    ( , .                                                                                                                                             
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of genomes or metagenomes in FASTA format")
    parser.add_argument('-bin_ext', type=str, help="filename extension for genomes or metagenomes")
    parser.add_argument('-out', type=str, help="output directory (will be created if does not exist)",
                        default="genie_out")
    parser.add_argument('-contigs_source', type=str, help="By default Lucifer runs Prodigal to predict protein-coding regions. "
                                                          "Please specify whether you are providing a directory of genomes "
                                                          "or metagenomes (default=genomes)", default="genomes")
    parser.add_argument('--makeplots', type=str,
                        help="Would you like Lucifer to make some figures from your data? y = yes, n = no (default = n). "
                             "If so, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Warning: this part of the program is currently under beta-testing, and if there are any problems running Rscript, "
                             "or installing any of the required packages, you may get a bunch of error messages at the end. "
                             "This will not crash the program, however, and you can still expect to "
                             "see the main output (CSV files) from Genie.", default="n")
    parser.add_argument('-bams', type=str, help="a tab-delimeted file with two columns: first column has the genome or "
                                                "metagenome file names; second column has the corresponding BAM file "
                                                "(provide full path to the BAM file). BAM files are only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts.", default="NA")
    parser.add_argument('--d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=3)", default=3)

    # CHECKING FOR CONDA INSTALL
    os.system("echo ${lux_hmms}/bitscores.txt > HMMlib.txt")
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

    args = parser.parse_known_args()[0]

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
        print("Looks like you did not provide a directory of ORFs in FASTA format")
        print("Exiting")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print('Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
              ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit

    if args.makeplots == 'y':
        if conda == 0:
            if args.R != "NA":
                print(".")
            else:
                print('Looks like you told Lucifer to automatically generate R plots. '
                      'However, you have not provided the location of the directory that contains the R scripts '
                      '(as required of you because you did not go through the conda-based installation.')
                print("Exiting")
                raise SystemExit

    if conda == 1:
        os.system("echo ${lux_hmms} > HMMlib.txt")
        file = open("HMMlib.txt")
        for i in file:
            HMMdir = (i.rstrip())
            HMMdirLS = os.listdir(HMMdir)
        os.system("rm HMMlib.txt")
    else:
        HMMdir = args.hmm_dir
        HMMdirLS = os.listdir(args.hmm_dir)

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)
        print("To avoid overwriting potentially valuable files, FeGenie will now exit. "
              "Please delete or rename said directory and try running again.")
        print("Exiting")
        raise SystemExit
    except FileNotFoundError:
        print(".")
        os.system("mkdir %s" % args.out)
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)

    numHMMs = 0
    for i in HMMdirLS:
        if lastItem(i.split(".")) == "hmm":
            numHMMs += 1

    # STARTING MAIN ALGORITHM
    bits = open(HMMdir + "/bitscores.txt", "r")
    bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    print("\nreading in HMM bitscore cut-offs...")
    for i in bits:
        ls = i.rstrip().split("\t")
        bitDict[ls[0]]["bit"] = ls[1]
        bitDict[ls[0]]["desc"] = ls[2]
    print("...")

    count = 0
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    out = open("%s/lucifer.csv" % (args.out), "w")
    out.write(
        "bin" + "," + "gene" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "seq" + "\n")

    outTSV = open("%s/lucifer.tsv" % (args.out), "w")
    outTSV.write(
        "bin" + "\t" + "gene" + "\t" + "ORF" + "\t" + "evalue" + "\t" + "bitscore" + "\t" + "bitscore_cutoff" + "\t" + "seq" + "\n")

    for i in binDirLS:
        HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.bin_ext:
            print("")
            try:
                testFile = open("%s/%s-proteins.faa" % (binDir, i), "r")
                print("")
                print(".")
                print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))

            except FileNotFoundError:
                print("")
                print(".")
                print("Finding ORFs for " + i)
                if args.contigs_source == "genomes":
                    os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -q" % (
                    binDir, i, binDir, i, binDir, i))
                elif args.contigs_source == "metagenomes":
                    os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -p meta -q" % (
                    binDir, i, binDir, i, binDir, i))
                else:
                    print("WARNING: you did not specify whether the provided FASTA files are single genomes or "
                          "metagenome/metatranscriptome assemblies. By default, FeGenie is assuming that these are "
                          "single genomes, and running Prodigal accordingly. Just an FYI.")
                    os.system("prodigal -i %s%s -a %s%s-proteins.faa -o %s%s-prodigal.out -q" % (
                    binDir, i, binDir, i, binDir, i))

            fastaFile = open(args.bin_dir + "/" + i + "-proteins.faa", "r")
            fastaFile = fasta(fastaFile)
            os.system("mkdir " + args.bin_dir + "/" + i + "-HMM")
            count = 0
            for hmm in HMMdirLS:
                if lastItem(hmm.split(".")) == "hmm":
                    count += 1
                    perc = (count / numHMMs) * 100 - 1
                    sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc + 1))
                    sys.stdout.flush()
                    if not re.match(r'^\.', hmm):
                        os.system("hmmsearch "
                                  "--tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
                                  " -o " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                                  HMMdir + "/" + hmm + " " +
                                  binDir + "/" + i + "-proteins.faa")
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

            for key in HMMdict.keys():
                out.write(i + "," + HMMdict[key]["hmm"] + "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key][
                    "bitscore"] + "," + bitDict[HMMdict[key]["hmm"]]["bit"] + "," + fastaFile[key] + "\n")

                outTSV.write(
                    i + "\t" + HMMdict[key]["hmm"] + "\t" + key + "\t" + HMMdict[key]["evalue"] + "\t" + HMMdict[key][
                        "bitscore"] + "\t" + bitDict[HMMdict[key]["hmm"]]["bit"] + "\t" + fastaFile[key] + "\n")

            os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")
    out.close()
    outTSV.close()


    print(".")
    summaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    summary = open("%s/lucifer.csv" % (args.out), "r")
    for i in summary:
        ls = i.rstrip().split(",")
        if ls != "bin":
            genome = ls[0]
            gene = ls[1]
            orf = ls[2]
            evalue = ls[3]
            bitscore = ls[4]
            bitcut = ls[5]
            sequence = ls[6]
            summaryDict[genome][orf]["gene"] = gene
            summaryDict[genome][orf]["evalue"] = evalue
            summaryDict[genome][orf]["bitscore"] = bitscore
            summaryDict[genome][orf]["sequence"] = sequence
            summaryDict[genome][orf]["bitcut"] = bitcut

    # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
    print(".")
    print("Identifying genomic proximities and putative operons")
    CoordDict = defaultdict(lambda: defaultdict(list))
    for i in summaryDict.keys():
        if i != "bin":
            for j in summaryDict[i]:
                contig = allButTheLast(j, "_")
                numOrf = lastItem(j.split("_"))
                CoordDict[i][contig].append(int(numOrf))

    counter = 0
    print(".")
    print("Clustering ORFs...")
    print(".")
    out = open("%s/lucifer-2.csv" % (args.out), "w")
    outTSV = open("%s/lucifer-2.tsv" % (args.out), "w")
    for i in CoordDict.keys():
        for j in CoordDict[i]:
            LS = (CoordDict[i][j])
            clusters = (cluster(LS, args.d))
            for k in clusters:
                if len(RemoveDuplicates(k)) == 1:
                    orf = j + "_" + str(k[0])

                    out.write(summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf]["evalue"] +
                              "," + str(summaryDict[i][orf]["bitscore"]) + "," + str(
                        summaryDict[i][orf]["bitcut"]) + "," +
                              str(summaryDict[i][orf]["sequence"]) + "," + str(counter) + "\n")

                    out.write(
                        "###############################################" + "\n")
                    counter += 1

                    outTSV.write(summaryDict[i][orf]["gene"] + "\t" + i + "\t" + orf + "\t" + summaryDict[i][orf]["evalue"] +
                              "\t" + str(summaryDict[i][orf]["bitscore"]) + "\t" + str(
                        summaryDict[i][orf]["bitcut"]) + "\t" +
                              str(summaryDict[i][orf]["sequence"]) + "\t" + str(counter) + "\n")

                    outTSV.write("################################################" + "\n")

                else:
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)

                        out.write(summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf][
                                      "evalue"] +
                                  "," + str(summaryDict[i][orf]["bitscore"]) + "," + str(
                            summaryDict[i][orf]["bitcut"]) +
                                  "," + str(summaryDict[i][orf]["sequence"]) + "," + str(counter) + "\n")

                        outTSV.write(summaryDict[i][orf]["gene"] + "\t" + i + "\t" + orf + "\t" + summaryDict[i][orf][
                                      "evalue"] +
                                  "\t" + str(summaryDict[i][orf]["bitscore"]) + "\t" + str(
                            summaryDict[i][orf]["bitcut"]) +
                                  "\t" + str(summaryDict[i][orf]["sequence"]) + "\t" + str(counter) + "\n")

                    out.write(
                        "################################################" + "\n")
                    counter += 1
                    outTSV.write(
                        "################################################" + "\n")


    out.close()

    os.system("rm %s/lucifer.csv" % (args.out))
    os.system("rm %s/lucifer.tsv" % (args.out))
    os.system("mv %s/lucifer-2.csv %s/lucifer-summary.csv" % (args.out, args.out))
    os.system("mv %s/lucifer-2.tsv %s/lucifer-summary.tsv" % (args.out, args.out))


    # # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    # print("Writing heatmap-compatible CSV")
    cats = ["Proteorhodopsin", "Xantharhodopsin", "Heliorhodopsin", "Heliorhodopsin_Pfam", "Bac_rhodopsin",
            "GpcrRhopsn4", "Htr2", "7tm_1", "ASRT", "Photo_RC", "PsaA_PsaB", "PSII", "PCP", "Chloroa_b-bind",
            "Bac_chlorC", "BChl_A", "PsbH", "LuxAB", "LuxC", "LuxD", "LuxE"]


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
            os.system("jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (args.out, cell, string))
            depth = open("%s/%s.depth" % (args.out, cell))
            for k in depth:
                LS = k.rstrip().split("\t")
                if LS[0] != "contigName":
                    depthDict[cell][LS[0]] = LS[2]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/lucifer-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if not re.match(r'#', i):
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        gene = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        Dict[cell][gene].append(float(depthDict[cell][contig]))


        outHeat = open("%s/lucifer.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            if (not re.match(r'#', i) and i != "bin"):
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
    else:
        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/lucifer-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if not re.match(r'#', i):
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        gene = ls[0]
                        Dict[cell][gene].append(gene)

        normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.bin_dir):
            if lastItem(i.split(".")) == args.bin_ext:
                file = open("%s/%s-proteins.faa" % (args.bin_dir, i), "r")
                file = fasta(file)
                normDict[i] = len(file.keys())

        outHeat = open("%s/lucifer.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            if (not re.match(r'#', i) and i != "bin"):
                outHeat.write(i + ",")
        outHeat.write("\n")

        for i in cats:
            outHeat.write(bitDict[i]["desc"] + ",")
            for j in sorted(Dict.keys()):
                if (not re.match(r'#', j) and j != "bin"):
                    outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
            outHeat.write("\n")

        outHeat.close()
        print('......')
        print(".......")
        print("Finished!")

    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots == "y":
        if conda == 0:
            Rdir = args.R
        else:
            os.system("echo ${rscripts} > r.txt")
            file = open("r.txt")
            for i in file:
                Rdir = (i.rstrip())
            os.system("rm r.txt")

        os.system("Rscript -e 'install.packages(\"ggplot2\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"reshape\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"reshape2\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"tidyverse\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"argparse\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"ggdendro\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"ggpubr\", repos = \"http://cran.us.r-project.org\")\'")
        os.system("Rscript -e 'install.packages(\"grid\", repos = \"http://cran.us.r-project.org\")\'")

        os.system("Rscript --vanilla %s/DotPlot.R %s/lucifer.%s.heatmap.csv %s" % (
        Rdir, args.out, args.element, args.out))
        os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.%s.heatmap.csv %s" % (
        Rdir, args.out, args.element, args.out))


if __name__ == '__main__':
    main()

