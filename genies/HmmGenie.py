#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def main():
    def SUM(ls):
        count = 0
        for i in ls:
            count += float(i)
        return count


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
                    seq = ''
                else:
                    header = i[1:]
                    seq = ''
            else:
                seq += i
        Dict[header] = seq
        # print(count)
        return Dict

    def fastaRename(fasta_file):
        counter = 0
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
                    counter += 1
                    header = header + "_" + str(counter)
                    seq = ''
                else:
                    header = i[1:]
                    header = header.split(" ")[0]
                    counter += 1
                    header = header + "_" + str(counter)
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


    def compare(queryList, targetList):
        counter = 0
        for i in targetList:
            if i in queryList:
                counter += 1
        if counter == len(targetList):
            return True
        else:
            return False


    parser = argparse.ArgumentParser(
        prog="MagicLamp.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************

        Developed by Arkadiy Garber and Nancy Merino;
        University of Southern California, Earth Sciences
        Please send comments and inquiries to arkadiyg@usc.edu

      )\  .'(   )\   )\   )\   )\     )\.-.    )\.---.   )\  )\  .'(   )\.---.  
     ,') \  ) (  ',/ /  (  ',/ /   ,' ,-,_)  (   ,-._( (  \, /  \  ) (   ,-._( 
    (  '-' (   )    (    )    (   (  .   __   \  '-,    ) \ (   ) (   \  '-,   
     ) .-.  ) (  \(\ \  (  \(\ \   ) '._\ _)   ) ,-`   ( ( \ \  \  )   ) ,-`   
    (  ,  ) \  `.) /  )  `.) /  ) (  ,   (    (  ``-.   `.)/  )  ) \  (  ``-.  
     )/    )/      '.(       '.(   )/'._.'     )..-.(      '.(    )/   )..-.(                                                                              
                                  %(?/////////&//%                                                
              .,,.                   (%((&@@@#/*.                      .,,.        
              .,,.                     @(((/&@@@#///**                  ...        
                                         #&((///////////////*/@                                
                                                             #*@.                             
                                      ()                   * )//*
                                      <^^>             *     (/*   .
                                     .-""-.                  *)
                          .---.    ."-....-"-._     _...---''`/. '
                         ( (`\ \ .'            ``-''    _.-"'`
                          \ \ \ : :.                 .-'
                           `\`.\: `:.             _.'
                           (  .'`.`            _.'
                            ``    `-..______.-'
                                      ):.  (
                                    ."-....-".
                                  .':.        `.
                                  "-..______..-"

        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/
        https://ascii.co.uk/text
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

    parser.add_argument('-bin_ext', type=str, help="file name extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=5)", default=5)

    parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format",
                        default="NA")

    parser.add_argument('-out', type=str, help="name output directory (default=hmmgenie_out)",
                        default="hmmgenie_out")

    parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                        default=1000)

    parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                             "(default=1, max=16)", default=1)

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
                                               "corresponding to different genomes that you are providing, please use the \'-bams\' "
                                               "argument to provide a tab-delimited file that denotes which BAM file (or files) belongs "
                                               "with which genome", default="NA")

    parser.add_argument('-hmm_dir', type=str, help="directory of HMMs that you want HmmGenie to profile against your dataset", default="NA")

    parser.add_argument('-hmm_ext', type=str, help="filename extension for the HMM files (e.g. hmm, txt, default = hmm)", default="hmm")

    parser.add_argument('-rules', type=str, help="optional file containing rules for HMM detection/reporting. The template"
                                                 "for this file is provided in the main MagicLamp directory (called rules-template.csv). "
                                                 "If you provide this file, HmmGenie will ignore any arguments provided to -eval and -clu flags", default="NA")

    parser.add_argument('-eval', type=str, help="e-value cutoff for hmmsearch. Default = 1E-10.", default=float(1E-10))

    parser.add_argument('-bit', type=str, help="bit score cutoff for hmmsearch. Default = 20.", default=float(20))

    parser.add_argument('-clu', type=str, help="minimum size of a gene cluster/operon to be considered for reporting (default = 1).",
                        default=int(1))

    parser.add_argument('--tree', type=str,
                        help="Provide this flag if you want HmmGenie to make a phylogenetic tree. "
                             "In order to take advantage of this part of the pipeline, please make sure that for each HMM file,"
                             "you have a corresponding alignment file. Must be the same name with only the extension different. "
                             "This different extension should be provided via the -aln_ext argument",
                        const=True, nargs="?")

    parser.add_argument('-aln_ext', type=str, help="filename extension for the alignment files. "
                                                   "Use this flag only if you set --tree flag (e.g. fasta, fa, default = fa)."
                                                   "These files must be in the same directory as the HMM files (same file"
                                                   "name as the HMM files, but with a different extension, given by this argument)",
                        default="fa")

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--orfs', type=str,
                        help="include this flag if you are providing bins as open-reading frames or genes in FASTA amino-acid format",
                        const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('-translation_table', type=str, help="specify a translation table for Prodigal to use (default: 11)", default="11")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each iron gene category to be normalized to "
                             "the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, HmmGenie will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, HmmGenie will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--dot', type=str,
                        help="invoke this flag if you would like HmmGenie to make a dot plot based on the data", const=True, nargs="?")

    parser.add_argument('--dendro', type=str,
                        help="invoke this flag if you would like HmmGenie to make a heatmap-dendrogram plot based on the data",
                        const=True, nargs="?")

    parser.add_argument('--word', type=str,
                        help="invoke this flag if you would like HmmGenie to make a word cloud based on the data",
                        const=True, nargs="?")

    parser.add_argument('--phobius', type=str,
                        help="invoke this flag if you would like HmmGenie to use Phobius for signal peptide and transmembrane prediction",
                        const=True, nargs="?")

    # parser.add_argument('--makeplots', type=str,
    #                     help="include this flag if you would like HmmGenie to make some figures from your data?. "
    #                          "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
    #                          "Please be sure to install all the required R packages as instrcuted in the HmmGenie Wiki: "
    #                          "https://github.com/Arkadiy-Garber/HmmGenie/wiki/Installation. "
    #                          "If you see error or warning messages associated with Rscript, you can still expect to "
    #                          "see the main output (CSV files) from HmmGenie.", const=True, nargs="?")

    # CHECKING FOR CONDA INSTALL
    os.system("echo ${rscripts} > rscripts.txt")

    file = open("rscripts.txt")
    for i in file:
        rscriptDir = i.rstrip()

    rscript = rscriptDir + "/" + "DotPlot.R"
    try:
        test = open(rscript)

    except FileNotFoundError:
        os.system("which MagicLamp.py > mainDir.txt")

        file = open("mainDir.txt")
        for i in file:
            location = i.rstrip()
        location = allButTheLast(location, "/")

        rscriptDir = location + "/rscripts/"
        rscript = rscriptDir + "/" + "DotPlot.R"

        os.system("rm rscripts.txt mainDir.txt")
        try:
            test = open(rscript)
        except FileNotFoundError:
            print("MagicLamp script could not locate the required directories. Please run the setup.sh script if \n"
                  "you have Conda installed. Otherwise, please run the setupe-noconda.sh script and put MagicLamp.py \n"
                  "into your $PATH")
            raise SystemExit
        os.system("rm mainDir.txt")

    os.system("rm rscripts.txt")

    args = parser.parse_known_args()[0]

    # ************** Checking for the required arguments ******************* #
    cwd = os.getcwd()
    print("checking arguments")

    if args.bin_dir != "NA":
        binDir = args.bin_dir + "/"
        binDirLS = os.listdir(args.bin_dir)
        print(".")
    else:
        print("Looks like you did not provide a directory of genomes/bins or assemblies.")
        print("Exiting")
        raise SystemExit

    if args.bam != "NA" and args.bams != "NA":
        print("Please provide only one of the following flags: \'-bam\' or \'-bams\'.")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print(
            'Looks like you did not provide an extension for your genomes/bins or assemblies, so HmmGenie does not know'
            ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)

        answer = input("Would you like HmmGenie to proceed and potentially overwrite files in this directory? (y/n): ")
        if answer == "y":
            print("Ok, proceeding with analysis!")
            try:
                os.listdir(args.out + "/ORF_calls")
            except FileNotFoundError:
                os.system("mkdir %s/ORF_calls" % args.out)
        else:
            print("Exiting")
            raise SystemExit

    except FileNotFoundError:
        print(".")
        os.system("mkdir %s" % args.out)
        os.system("mkdir %s/ORF_calls" % args.out)

    if lastItem(args.out) == "/":
        outDirectory = "%s" % args.out[0:len(args.out)-1]
        outDirectoryLS = os.listdir(outDirectory)
    else:
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)

    print("All required arguments provided!")
    print("")

    # *************** MAKE NR A DIAMOND DB AND READ THE FILE INTO HASH MEMORY ************************ #
    if args.ref != "NA":
        try:
            testFile = open(args.ref + ".dmnd")

        except FileNotFoundError:
            print("Making diamond database out of provided reference file")
            os.system("diamond makedb --in %s -d %s" % (args.ref, args.ref))

    # *************** CALL ORFS FROM BINS AND READ THE ORFS INTO HASH MEMORY ************************ #
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in binDirLS:
        print(i)
        if lastItem(i.split(".")) == args.bin_ext:
            cell = i
            if not args.gbk:

                if args.orfs:
                    testFile = open("%s/%s" % (binDir, i), "r")
                    for line in testFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \|")
                                print(
                                    "Unfortunately, this is a problem for HmmGenie because it uses that character as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                else:
                    try:
                        testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                        print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                        for line in testFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print(
                                        "Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for HmmGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                    except FileNotFoundError:
                        binFile = open("%s/%s" % (binDir, i), "r")
                        for line in binFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for HmmGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                        print("Finding ORFs for " + cell)
                        if args.meta:
                            os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -p meta -q -g %s" % (
                                binDir, i, outDirectory, i, outDirectory, i, args.translation_table))
                        else:
                            os.system(
                                "prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -q -g %s" % (
                                    binDir, i, outDirectory, i, outDirectory, i, args.translation_table))
            else:
                os.system('gtt-genbank-to-AA-seqs -i %s/%s -o %s/%s.faa' % (binDir, i, outDirectory, i))

                faa = open("%s/%s.faa" % (binDir, i))
                faa = fasta(faa)

                gbkDict = defaultdict(list)
                counter = 0

                count = 0
                gbk = open("%s/%s" % (binDir, i))
                for gbkline in gbk:
                    ls = gbkline.rstrip()
                    if re.findall(r'/locus_tag', ls):
                        count += 1

                if count > 0:
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])
                            counter += 1

                        if counter > 0:
                            gbkDict[locus].append(locusTag)
                            counter = 0
                else:
                    # print(i)
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)
                            counter += 1

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])

                        if counter > 0:
                            gbkDict[locus].append(altContigName)
                            counter = 0

                idxOut = open("%s/ORF_calls/%s-proteins.idx" % (outDirectory, i), "w")
                faaOut = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "w")

                for gbkkey1 in gbkDict.keys():
                    counter = 0
                    for gbkey2 in gbkDict[gbkkey1]:
                        counter += 1
                        if len(faa[gbkey2]) > 0:
                            newOrf = gbkkey1 + "_" + str(counter)
                            idxOut.write(gbkey2 + "," + newOrf + "\n")
                            faaOut.write(">" + newOrf + "\n")
                            faaOut.write(str(faa[gbkey2]) + "\n")

                idxOut.close()
                faaOut.close()

            if args.orfs:
                file = open("%s/%s" % (binDir, i))
            else:
                file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
            file = fasta(file)
            for j in file.keys():
                orf = j.split(" # ")[0]
                BinDict[cell][orf] = file[j]

    # ******************** READ BITSCORE CUT-OFFS INTO HASH MEMORY ****************************** #
    if args.rules != "NA":
        meta = open(args.rules, "r")
        metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        operonDict = defaultdict(list)
        operonMetaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        operonMainDict = defaultdict(list)
        eDict = defaultdict(lambda: 'EMPTY')
        for i in meta:
            ls = i.rstrip().split(",")
            if ls[1] != "output gene name":
                hmm = ls[0]
                gene = ls[1]
                operon = ls[2]
                importance = ls[3]
                minHits = int(ls[4])
                bitcut = float(ls[5])
                evalue = float(ls[6])
                minlength = int(ls[7])
                maxlength = int(ls[8])
                metaDict[hmm]["gene"] = gene
                metaDict[hmm]["operon"] = operon
                metaDict[hmm]["minHits"] = minHits
                metaDict[hmm]["bitcut"] = bitcut
                metaDict[hmm]["evalue"] = evalue
                metaDict[hmm]["minlength"] = minlength
                metaDict[hmm]["maxlength"] = maxlength
                operonDict[operon].append(hmm)
                operonMetaDict[operon] = minHits
                if importance == "1":
                    operonMainDict[operon].append(hmm)
    else:
        metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.hmm_dir):
            if lastItem(i.split(".")) == remove(args.hmm_ext, ["."]):
                hmm = i
                metaDict[hmm]["gene"] = allButTheLast(hmm, ".")
                metaDict[hmm]["operon"] = "unnamed_operon"
                metaDict[hmm]["minHits"] = 1
                metaDict[hmm]["bitcut"] = float(args.bit)
                metaDict[hmm]["evalue"] = float(args.eval)
                metaDict[hmm]["minlength"] = 50
                metaDict[hmm]["maxlength"] = 10000

    # ******************* BEGINNING MAIN ALGORITHM **********************************))))
    print("\nStarting main pipeline...")
    HMMdirLS = os.listdir(args.hmm_dir)
    HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
    for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
        if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
            cell = i
            os.system(
                "mkdir -p " + outDirectory + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

            count = 0

            numHMMs = 0
            for hmm in HMMdirLS:
                if lastItem(hmm.split(".")) == args.hmm_ext:
                    numHMMs += 1

            for hmm in HMMdirLS:  # ITERATING THROUGH ALL THE HMM FILES IN THE HMM DIRECTORY
                if lastItem(hmm.split(".")) == remove(args.hmm_ext, ["."]):
                    if hmm not in metaDict.keys():
                        print("did not detect %s in the rules.csv file. Using the defaults for evalue and bit score cutoffs" % (hmm))
                        metaDict[hmm]["evalue"] = args.eval
                        metaDict[hmm]["bitcut"] = args.bit

                    count += 1
                    perc = (count / numHMMs) * 100
                    sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                    sys.stdout.flush()

                    if args.orfs:
                        os.system(
                            "hmmsearch -E %s --cpu %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/%s"
                            % (metaDict[hmm]["evalue"], int(args.t), outDirectory, i, hmm, outDirectory, i, hmm, args.hmm_dir, hmm, binDir, i))
                    else:
                        os.system(
                            "hmmsearch -E %s --cpu %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                            % (metaDict[hmm]["evalue"], int(args.t), outDirectory, i, hmm, outDirectory, i, hmm, args.hmm_dir, hmm, outDirectory, i))

                    # REMOVING THE STANDARD OUTPUT FILE
                    os.system("rm " + outDirectory + "/" + i + "-HMM/" + hmm + ".txt")

                    # READING IN THE HMMSEARCH RESULTS (TBLOUT) OUT FILE
                    hmmout = open(outDirectory + "/" + i + "-HMM/" + hmm + ".tblout", "r")

                    # COLLECTING SIGNIFICANT HMM HITS IN THE FILE
                    for line in hmmout:
                        if not re.match(r'#', line):
                            ls = delim(line)
                            evalue = float(ls[4])
                            bit = float(ls[5])
                            bitcut = float(metaDict[hmm]["bitcut"])
                            orf = ls[0]
                            seq = BinDict[cell][orf]
                            # print("")
                            # print(len(seq))
                            # print(metaDict[hmm]["minlength"])
                            # print(metaDict[hmm]["maxlength"])
                            # print("+")
                            # print(bit)
                            # print(bitcut)

                            if bit > bitcut and len(seq) >= metaDict[hmm]["minlength"] and len(seq) <= metaDict[hmm]["maxlength"]:
                                # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS
                                if orf not in HMMdict[i]:
                                    HMMdict[i][orf]["hmm"] = hmm
                                    HMMdict[i][orf]["evalue"] = evalue
                                    HMMdict[i][orf]["bit"] = bit
                                    HMMdict[i][orf]["seq"] = seq
                                    HMMdict[i][orf]["bitcut"] = metaDict[hmm]["bitcut"]
                                else:
                                    # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                                    if evalue > HMMdict[i][orf]["evalue"]:
                                        HMMdict[i][orf]["hmm"] = hmm
                                        HMMdict[i][orf]["evalue"] = evalue
                                        HMMdict[i][orf]["bit"] = bit
                                        HMMdict[i][orf]["seq"] = seq
                                        HMMdict[i][orf]["bitcut"] = metaDict[hmm]["bitcut"]
        print("")

    out = open("%s/summary.csv" % (outDirectory), "w")
    out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "seq" + "\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                      str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) +
                      "," + str(HMMdict[key][j]["bitcut"]) + "," + str(HMMdict[key][j]["seq"]) + "\n")

    out.close()
    # ****************************************** DEREPLICATION *********************************************************
    summary = open(outDirectory + "/summary.csv", "r")
    SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] != "category" and ls[0] != "HmmGenie":
            if len(ls) > 0:
                category = ls[0]
                cell = ls[0]
                orf = ls[1]
                hmm = ls[2]
                evalue = ls[3]
                hmmBit = ls[4]
                bitcut = ls[5]
                seq = ls[6]
                if cell not in SummaryDict.keys():
                    SummaryDict[cell][orf]["hmm"] = hmm
                    SummaryDict[cell][orf]["hmmBit"] = hmmBit
                    SummaryDict[cell][orf]["category"] = category
                    SummaryDict[cell][orf]["e"] = evalue
                    SummaryDict[cell][orf]["bitcut"] = bitcut
                    SummaryDict[cell][orf]["seq"] = seq

                else:
                    if orf not in SummaryDict[cell]:
                        SummaryDict[cell][orf]["hmm"] = hmm
                        SummaryDict[cell][orf]["hmmBit"] = hmmBit
                        SummaryDict[cell][orf]["category"] = category
                        SummaryDict[cell][orf]["e"] = evalue
                        SummaryDict[cell][orf]["bitcut"] = bitcut
                        SummaryDict[cell][orf]["seq"] = seq

                    else:
                        if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                            SummaryDict[cell][orf]["hmm"] = hmm
                            SummaryDict[cell][orf]["hmmBit"] = hmmBit
                            SummaryDict[cell][orf]["category"] = category
                            SummaryDict[cell][orf]["e"] = evalue
                            SummaryDict[cell][orf]["bitcut"] = bitcut
                            SummaryDict[cell][orf]["seq"] = seq

    # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
    print("\nIdentifying genomic proximities and putative operons")
    CoordDict = defaultdict(lambda: defaultdict(list))
    for i in SummaryDict.keys():
        if i != "cell":
            for j in SummaryDict[i]:
                contig = allButTheLast(j, "_")
                numOrf = lastItem(j.split("_"))
                CoordDict[i][contig].append(int(numOrf))

    counter = 0
    print("Clustering ORFs...")
    print("")
    out = open(outDirectory + "/summary-2.csv", "w")
    for i in CoordDict.keys():
        # print(".")
        for j in CoordDict[i]:
            LS = (CoordDict[i][j])
            clusters = (cluster(LS, args.d))
            for k in clusters:
                if len(RemoveDuplicates(k)) >= int(args.clu):
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)
                        out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                    out.write("###############################################\n")
                    counter += 1
    out.close()

    os.system("rm %s/summary.csv" % (args.out))

    os.system("mkdir -p %s/HMM_results" % outDirectory)
    # os.system("rm -f %s/ORF_calls/*-prodigal.out" % outDirectory)
    os.system("rm -rf %s/HMM_results/*-HMM" % outDirectory)
    os.system("mv -f %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

    # ****************************** CUSTOM-RULE-BASED FILTERING ************************************************
    if args.rules != "NA":
        clusterDict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        summary = open("%s/summary-2.csv" % args.out)
        out = open("%s/summary-3.csv" % args.out, "w")
        out.write("genome/metagenome,gene_call,hmm_file,gene,e_value,bit_score,bit_score_cutoff,cluster_id,aa_seq\n")
        for i in summary:
            ls = i.rstrip().split(",")

            if not re.match(r'#', i.rstrip()):
                hmmFile = ls[2]
                clusterNum = ls[6]
                clusterDict[ls[0]][clusterNum]["ls"].append(ls)
                clusterDict[ls[0]][clusterNum]["hmms"].append(hmmFile)
                if hmmFile not in clusterDict[ls[0]][clusterNum]["unqHmms"]:
                    clusterDict[ls[0]][clusterNum]["unqHmms"].append(hmmFile)

        masterDict = defaultdict(lambda: defaultdict(list))
        for i in clusterDict.keys():
            # print(i)
            for k in clusterDict[i]:
                # print(k)
                clusterNum = k
                ls = clusterDict[i][k]["ls"]
                hmms = (clusterDict[i][k]["hmms"])
                unqHmms = (clusterDict[i][k]["unqHmms"])
                operons = []
                for j in ls:
                    operon = metaDict[j[2]]["operon"]
                    if operon not in operons:
                        operons.append(operon)
                passDict = defaultdict(list)
                for j in operons:
                    operon = j
                    minHits = operonMetaDict[operon]
                    if len(unqHmms) >= int(minHits) and compare(unqHmms, operonMainDict[operon]):
                        passDict[operon].append(ls[0][6])
                for j in ls:
                    operon = metaDict[j[2]]["operon"]
                    if operon in passDict:
                        if j[6] in passDict[operon]:
                            masterDict[j[0]][clusterNum].append(j)

        for i in masterDict.keys():
            for j in masterDict[i]:
                for k in masterDict[i][j]:
                    # out.write(",".join(k) + "\n")
                    out.write(k[0] + "," + k[1] + "," + k[2] + "," + str(metaDict[k[2]]["gene"]) + "," + k[3] + "," + k[4] + "," + k[5] + "," + k[6] + "," + k[7] + "\n")
                out.write("#########################################\n")
            out.write("#********************************************************************\n")
            out.write("#********************************************************************\n")
            out.write("#********************************************************************\n")
        out.close()
        os.system("mv %s/summary-3.csv %s/genie-summary-rulesFiltered.csv" % (args.out, args.out))
        os.system("mv %s/summary-2.csv %s/genie-summary-allResults.csv" % (args.out, args.out))

    else:
        os.system("mv %s/summary-2.csv %s/genie-summary-allResults.csv" % (args.out, args.out))
        os.system("cp %s/genie-summary-allResults.csv %s/genie-summary-rulesFiltered.csv" % (args.out, args.out))


    # *********************************** CREATING A PHYLOGENETIC TREE ******************************************
    if args.tree:
        print("\nWorking on the phylogenetic tree:")
        # checking for presence of correct files
        alnDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in sorted(HMMdirLS):

            if lastItem(i.split(".")) == args.hmm_ext:
                base = allButTheLast(i, ".")
                alnDict[base]["hmm"] = i

            if lastItem(i.split(".")) == args.aln_ext:
                base = allButTheLast(i, ".")
                alnDict[base]["fa"] = i

        alnDict2 = defaultdict(list)
        for i in alnDict.keys():
            if len(alnDict[i]) == 2:
                print(i + " has an associated HMM and alignment file")
                alnDict2[i] = alnDict[i]
            else:
                print(i + " does not have an associated HMM and alignment file and will be included as a phylo tree")

        if len(alnDict2.keys()) == 0:
            print("HmmGenie did not detect any alignment files. Aborting tree construction.")

        else:
            phyloDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            summary = open("%s/genie-summary-rulesFiltered.csv" % args.out)
            for i in summary:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    if ls[1] != "gene_call":
                        base = allButTheLast(ls[2], ".")
                        if base in alnDict2.keys():
                            phyloDict[base][ls[3] + "-" + ls[0] + "-" + ls[1]] = ls[8]

            for i in HMMdirLS:
                if lastItem(i.split(".")) == args.aln_ext:
                    base = allButTheLast(i, ".")
                    if base in alnDict2.keys():
                        alnFile = open("%s/%s.%s" % (args.hmm_dir, base, args.aln_ext))
                        alnFile = fasta(alnFile)
                        for j in alnFile.keys():
                            phyloDict[base][j] = alnFile[j]

            for i in phyloDict.keys():
                out = open("%s/%s.faa" % (args.out, i), "w")
                for j in phyloDict[i]:
                    out.write(">" + j + "\n")
                    out.write(phyloDict[i][j] + "\n")
                out.close()
                print("\nAligning %s HMM hits to seed sequences..." % i)
                os.system("muscle -in %s/%s.faa -out %s/%s.fa > /dev/null 2>&1" % (args.out, i, args.out, i))
                print("Building phylogenetic tree for %s" % i)
                os.system("fasttree %s/%s.fa > %s/%s.tre > /dev/null 2>&1" % (args.out, i, args.out, i))

            os.system("mkdir -p %s/trees" % args.out)
            os.system("mv %s/*fa %s/trees/" % (args.out, args.out))
            os.system("mv %s/*tre %s/trees/" % (args.out, args.out))

    if args.word:

        Rdir = "purgatory"
        os.system("echo ${rscripts} > r.txt")
        Rfile = open("r.txt")
        for i in Rfile:
            Rdir = (i.rstrip())
        os.system("rm r.txt")

        try:
            test = open(Rdir + "/wordcloud.R")
            word = 1

        except FileNotFoundError:
            os.system("which MagicLamp.py > r.txt")
            Rfile = open("r.txt")
            for i in Rfile:
                Rdir = (i.rstrip())
            Rdir = allButTheLast(Rdir, "/")
            os.system("rm r.txt")
            try:
                test = open(Rdir + "/wordcloud.R")
                word = 1
            except FileNotFoundError:
                print("You have directed HmmGenie to make a word-cloud. However, HmmGenie cannot seem to find the required"
                      "R script. This script can be found in the main MagicLamp directory, titled wordcloud.R. ")
                answer = input(
                    "If you would like HmmGenie to generate the word-clouds, you will need to provide the full path to this script. "
                    "Would you like to provide the full path? (y/n): ")
                if answer == "y":
                    path = input("Please provide the full path to wordcloud.R: ")
                    try:
                        test = open(path)
                        word = 1
                        Rdir = allButTheLast(path, "/")
                    except FileNotFoundError:
                        print("Hmm, HmmGenie still cannot seem to locate the R script.")
                        answer = input("Please check your path and try again? (y/n): ")
                        if answer == "y":
                            path = input("Please provide the full path to wordcloud.R: ")
                            try:
                                test = open(path)
                                word = 1
                                Rdir = allButTheLast(path, "/")
                            except FileNotFoundError:
                                print("Hmm, HmmGenie still cannot seem to locate the R script. Moving on without the word-cloud"
                                      "for now. If you feel this is in error, please start on Issue on MagicLamp's GitHub repository.")
                                word = 0
                        else:
                            print("Very well. Moving on without word-clouds")
                            word = 0

                else:
                    print("Very well. Moving on without word-clouds")
                    word = 0

        if word == 1:
            print("Working on the word clouds")
            wordDict = defaultdict(list)
            summary = open("%s/genie-summary-rulesFiltered.csv" % args.out)
            for i in summary:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    if ls[1] != "gene_call":
                        wordDict[ls[0]].append(ls[3])

            for i in wordDict.keys():
                out = open("%s/%s.words.csv" % (args.out, allButTheLast(i, ".")), "w")
                for j in wordDict[i]:
                    out.write(j + "\n")
                out.close()

                os.system("Rscript --vanilla %s/wordcloud.R %s/%s.words.csv %s/%s.words.tiff > /dev/null 2>&1" % (
                Rdir, args.out, allButTheLast(i, "."), args.out, allButTheLast(i, ".")))

            # print("Cleaning up")
            os.system("mkdir -p %s/wordClouds" % args.out)
            os.system("mv %s/*words* %s/wordClouds/" % (args.out, args.out))

    if args.phobius:
        out = open("%s/genie-seqs.faa" % args.out, "w")
        summary = open("%s/genie-summary-rulesFiltered.csv" % args.out)
        for i in summary:
            if not re.match(r'#', i):
                ls = i.rstrip().split(",")
                if ls[1] != "gene_call":
                    out.write(">" + ls[0] + "-" + ls[1] + "\n")
                    out.write(remove(ls[8], ["*"]) + "\n")

        os.system("phobius.pl -short %s/genie-seqs.faa > %s/genie-seqs.phobius" % (args.out, args.out))
        phobiusDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        phobius = open("%s/genie-seqs.phobius" % args.out)
        for i in phobius:
            ls = delim(i.rstrip())
            phobiusDict[ls[0]]["tm"] = ls[1]
            phobiusDict[ls[0]]["sp"] = ls[2]

        summary = open("%s/genie-summary-rulesFiltered.csv" % args.out)
        out = open("%s/genie-summary-rulesFiltered-2.csv" % args.out, "w")
        for i in summary:
            if not re.match(r'#', i):
                ls = i.rstrip().split(",")
                if ls[1] != "gene_call":
                    out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + ls[5] + "," + ls[6] + "," + ls[7] + "," + str(phobiusDict[ls[0] + "-" + ls[1]]["tm"]) + "," + str(phobiusDict[ls[0] + "-" + ls[1]]["sp"]) + "," + ls[8] + "\n")
                else:
                    out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + ls[5] + "," + ls[6] + "," + ls[7] + "," + str("TM_domains") + "," + str("singal_peptide") + "," + ls[8] + "\n")
            else:
                out.write(i.rstrip() + "\n")
        out.close()

        os.system("mv %s/genie-summary-rulesFiltered-2.csv %s/genie-summary-rulesFiltered.csv" % (args.out, args.out))

    # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    cats = []
    for hmm in HMMdirLS:
        if lastItem(hmm.split(".")) == args.hmm_ext:
            cats.append(hmm)
    print("Working on the heatmap-format CSV")

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

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/genie-summary-rulesFiltered.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[cell][contig]))

        outHeat = open("%s/hmmgenie.readDepth.heatmap.csv" % (args.out), "w")
        outHeat.write("X")
        for i in sorted(Dict.keys()):
            outHeat.write("," + i)
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i)
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    outHeat.write("," + str(SUM(Dict[j][i])))
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

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/genie-summary-rulesFiltered.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[contig]))

        outHeat = open("%s/hmmgenie.readDepth.heatmap.csv" % (args.out), "w")
        outHeat.write("X")
        for i in sorted(Dict.keys()):
            outHeat.write("," + i)
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i)
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    outHeat.write("," + str(SUM(Dict[j][i])))
            outHeat.write("\n")

        outHeat.close()
        print('......')
        print(".......")
        print("Finished!")

    # GENE COUNTS-BASED ABUNDANCE
    else:
        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/genie-summary-rulesFiltered.csv" % (args.out), "r")
        for i in final:
            if not re.match(r'#', i):
                ls = (i.rstrip().split(","))
                print(ls)
                if ls[3] != "gene":
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    Dict[cell][gene].append(gene)

        normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.bin_dir):
            if lastItem(i.split(".")) == args.bin_ext:
                if args.orfs:
                    file = open("%s/%s" % (binDir, i), "r")
                    file = fasta(file)
                else:
                    file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    file = fasta(file)
                normDict[i] = len(file.keys())

        outHeat = open("%s/genie.heatmap.csv" % (outDirectory), "w")
        outHeat.write("X")
        for i in sorted(Dict.keys()):
            outHeat.write("," + i)
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i)
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    if args.norm:
                        outHeat.write("," + str((len(Dict[j][i]) / int(normDict[j])) * float(100)))
                    else:
                        outHeat.write("," + str((len(Dict[j][i]))))
            outHeat.write("\n")
        outHeat.close()

        if args.rules == "NA":
            os.system("rm %s/genie-summary-rulesFiltered.csv" % args.out)

        print('......')
        print(".......")
        print("Finished!")


    # # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************

    if args.bam == "NA" and args.bams == "NA":
        if args.norm:
            if args.dot:
                os.system("Rscript --vanilla %s/DotPlot.R %s/genie.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            if args.dendro:
                os.system("Rscript --vanilla %s/dendro-heatmap.R %s/genie.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
        else:
            if args.dot:
                os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/genie.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            if args.dendro:
                os.system("Rscript --vanilla %s/dendro-heatmap.R %s/genie.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
    else:
        if args.dot:
            os.system("Rscript --vanilla %s/DotPlot.R %s/genie.readDepth.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
        if args.dendro:
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/genie.readDepth.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))

    print("")
    print("Results are written to %s/genie-summary.csv and %s/genie.heatmap.csv" % (args.out, args.out))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()
