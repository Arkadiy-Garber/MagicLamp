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
        prog="LithoGenie.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************
    
        Developed by Arkadiy Garber and Gustavo Ramírez;
        University of Southern California, Earth Sciences
        Please send comments and inquiries to arkadiyg@usc.edu
        *******************************************************
        '''))

    parser.add_argument('-bins', type=str, help="directory of bins")

    parser.add_argument('-binx', type=str, help="filename extension for bins (do not include the period)")

    parser.add_argument('--meta', type=str,
                        help="use this flag if you are providing metagenome assemblies, rather than genomes. "
                             "If you are providing ORFs via the \'-orfs\' argument, ignore this flag.",
                        const=True, nargs="?")

    parser.add_argument('--orfs', type=str,
                        help="use this flag if you are providing ORFs (i.e. previously-made gene calls/predictions) instead of contigs ",
                        const=True, nargs="?")

    parser.add_argument('-hmms', type=str, help="directory of HMM files")

    parser.add_argument('-hmmx', type=str, help="filename extension for HMMs (do not include the period)")

    parser.add_argument('-bits', type=str, help="TSV file with bitscore cut-offs for each HMM in the HMM directory "
                                                "(provided via the \'-hmms\' argument). The file should contain two columns: "
                                                "the first should list the filename of each HMM, "
                                                "and the second column needs to have the trusted bitscores cut-offs for each "
                                                "HMM. If this is not provided, Genie will attempt to use the profile's "
                                                "trusted cutoff using the \'--cut_tc\' flag for hmmsearch.")

    parser.add_argument('-out', type=str, help="output directory (will be created if does not exist)",
                        default="genie_out")


    parser.add_argument('--norm', type=str,
                        help="normalize gene counts to total number of predicted ORFs in each genome or metagenome. "
                             "Ignore this flag if you are providing BAM files ",
                        const=True, nargs="?")

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

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                              "This number should be an integer and should reflect the maximum number of "
                                              "genes in between putative iron-related genes identified by the HMM database "
                                              "(default=3)", default=3)


    parser.add_argument('--skip', type=str, help="subvert all other functions of this program and only make the"
                                                      " heatmap. You can choose this flag if you would only like to change "
                                                 "how the data are normalized. Please provide the output directory using "
                                                 "the \'-out\' flag; this already must "
                                                 "contain the output summary file from a previously completed run",
                        const=True, nargs="?")

    parser.add_argument('--cpu', type=int, help="number of threads to allow for hmmsearch (default = 1)", default=1)

    args = parser.parse_known_args()[0]

    # CHECKING ARGUMENTS AND PATHS
    cwd = os.getcwd()
    print("checking arguments")

    if args.hmms != "NA":
        print(".")
        HMMdir = args.hmms
        HMMdirLS = os.listdir(args.hmms)
    else:
        print("You have not provided the location of the HMM library via the -hmm_dir argument. Please do so, and try "
              "again. The HMM library is found within the same directory as the FeGenie executable.")
        print("Exiting")
        raise SystemExit
    ##########################
    if args.bins != "NA":
        binDir = args.bins + "/"
        binDirLS = os.listdir(args.bins)
        print(".")
    else:
        print("Looks like you did not provide a directory of genomes/bins or assemblies.")
        print("Exiting")
        raise SystemExit
    ##########################
    if args.binx != "NA":
        print(".")
    else:
        print('Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
              ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit
    #########################
    if not args.skip:
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
    else:
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)
    #########################
    if args.bits:
        bits = open(args.hmms)
        bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in bits:
            ls = i.rstrip().split("\t")
            bitDict[ls[0]] = ls[1]
        print(".")
    else:
        bits = "cut_tc"
        print(".")
    #########################


    # STARTING MAIN ALGORITHM
    if not args.skip:
        print("Starting HMMER")
        count = 0
        BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        out = open("%s/summary.csv" % (args.out), "w")
        out.write(
            "bin" + "," + "hmm" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "sequence" + "\n")
        for i in binDirLS:
            HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.binx:

                if args.orfs:
                    seqdb = args.bins + "/" + i
                    fastaFile = open(args.bins + "/" + i, "r")
                    fastaFile = fasta(fastaFile)
                else:

                    if args.meta:
                        os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -p meta -q" % (
                            binDir, i, binDir, i, binDir, i))
                    else:
                        os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -q" % (
                            binDir, i, binDir, i, binDir, i))

                    seqdb = args.bins + "/" + i + "-proteins.faa"
                    fastaFile = open(args.bins + "/" + i + "-proteins.faa", "r")
                    fastaFile = fasta(fastaFile)

                os.system("mkdir " + args.bins + "/" + i + "-HMM")
                count = 0
                print("")
                for hmm in HMMdirLS:
                    if lastItem(hmm.split(".")) == args.hmmx:
                        count += 1
                        perc = (count / len(HMMdirLS)) * 100
                        if perc < 100:
                            sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc + 1))
                            sys.stdout.flush()
                        else:
                            sys.stdout.write("analyzing " + i + ": %d%%   \r" % (100))
                            sys.stdout.flush()
                        if args.bits:
                            os.system("hmmsearch -T " + str(bitDict[hmm]) + " --cpu " + str(args.cpu) +
                                      " --tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
                                      " -o " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                                      HMMdir + "/" + hmm + " " +
                                      seqdb)

                        else:
                            os.system("hmmsearch --cut_tc --cpu " + str(args.cpu) +
                                      " --tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
                                      " -o " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                                      HMMdir + "/" + hmm + " " +
                                      seqdb)

                        os.system("rm " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt")

                HMMresults = os.listdir(args.bins + "/" + i + "-HMM")
                for result in HMMresults:
                    if not re.match(r'^\.', result):
                        file = open(args.bins + "/" + i + "-HMM/" + result)
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

                                if orf not in HMMdict.keys():
                                    HMMdict[orf]["hmm"] = hmmName
                                    HMMdict[orf]["evalue"] = evalue
                                    HMMdict[orf]["bitscore"] = bitscore
                                else:
                                    if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                    else:
                                        pass

                for key in HMMdict.keys():
                    out.write(i + "," + HMMdict[key]["hmm"] + "," + key + "," + HMMdict[key]["evalue"] + "," +
                        HMMdict[key]["bitscore"] + "," + fastaFile[key] + "\n")

                os.system("rm -r " + args.bins + "/" + i + "-HMM")

        out.close()

        print("Finished HMMER")
        summaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
        summary = open("%s/summary.csv" % args.out, "r")
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[0] != "bin":
                genome = ls[0]
                gene = ls[1]
                orf = ls[2]
                evalue = ls[3]
                bitscore = ls[4]
                sequence = ls[5]
                summaryDict[genome][orf]["gene"] = gene
                summaryDict[genome][orf]["evalue"] = evalue
                summaryDict[genome][orf]["bitscore"] = bitscore
                summaryDict[genome][orf]["sequence"] = sequence

        print("Identifying genomic proximities and putative operons")
        CoordDict = defaultdict(lambda: defaultdict(list))
        for i in summaryDict.keys():
            if i != "bin":
                for j in summaryDict[i]:
                    contig = allButTheLast(j, "_")
                    numOrf = lastItem(j.split("_"))
                    CoordDict[i][contig].append(int(numOrf))

        counter = 0
        print("Clustering ORFs...")
        out = open("%s/summary-2.csv" % args.out, "w")
        out.write("hmm" + "," + "bin" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "cluster_id" + "," + "sequence" + "\n")
        for i in CoordDict.keys():
            for j in CoordDict[i]:
                LS = (CoordDict[i][j])
                clusters = (cluster(LS, args.d))
                for k in clusters:
                    if len(RemoveDuplicates(k)) == 1:
                        orf = j + "_" + str(k[0])

                        out.write(summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf]["evalue"] +
                                  "," + str(summaryDict[i][orf]["bitscore"]) + "," +  str(counter) + "," +
                                  str(summaryDict[i][orf]["sequence"]) + "\n")

                        out.write("#################################################################################\n")
                        counter += 1

                    else:
                        for l in RemoveDuplicates(k):
                            orf = j + "_" + str(l)

                            out.write(summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf]["evalue"] +
                                      "," + str(summaryDict[i][orf]["bitscore"]) + "," + str(counter) +
                                      "," + str(summaryDict[i][orf]["sequence"]) + "\n")

                            out.write("#################################################################################\n")
                        counter += 1

        out.close()

        os.system("mv %s/summary-2.csv %s/summary.csv" % (args.out, args.out))

        # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
        print("....")
        print(".....")
        # COVERAGE-BASED ABUNDANCE
        if args.bams != "NA":
            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            BAMmap = open(args.bams)
            normDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            for i in BAMmap:
                string = ''
                ls = i.rstrip().split("\t")
                cell = ls[0]
                for j in ls[1:]:
                    string += " "
                    string += j

                try:
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total / 1000000

                except FileNotFoundError:
                    os.system("mkdir -p %s/contigDepths" % args.out)
                    os.system("jgi_summarize_bam_contig_depths --outputDepth %s/contigDepths/%s.depth%s" % (args.out, cell, string))
                    print("processing... " + cell)
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total / 1000000

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(float(depthDict[cell][contig]))

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            outHeat = open("%s/readDepth.heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
                outHeat.write("\n")

            outHeat.close()

        # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE
        elif args.bam != "NA":
            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

            try:
                total = 0
                depth = open("%s/contigDepths/%s.depth" % (args.out, args.bam))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = float(LS[2])
                        total += float(LS[2])

            except FileNotFoundError:
                os.system("mkdir -p %s/contigDepths" % args.out)
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s/contigDepths/%s.depth %s" % (args.out, args.bam, args.bam))
                depth = open("%s/contigDepths/%s.depth" % (args.out, args.bam))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = float(LS[2])
                        total += float(LS[2])

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(float(depthDict[contig]))

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            outHeat = open("%s/readDepth.heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
                outHeat.write("\n")

            outHeat.close()

        # GENE COUNTS-BASED ABUNDANCE
        else:
            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(gene)

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in binDirLS:
                if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.binx:

                    if args.orfs:
                        seqdb = args.bins + "/" + i
                        fastaFile = open(args.bins + "/" + i, "r")
                        fastaFile = fasta(fastaFile)
                    else:

                        seqdb = args.bins + "/" + i + "-proteins.faa"
                        fastaFile = open(args.bins + "/" + i + "-proteins.faa", "r")
                        fastaFile = fasta(fastaFile)

                    normDict[i] = len(fastaFile.keys())

            outHeat = open("%s/heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):

                    if args.norm == "y":
                        outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
                    else:
                        outHeat.write(str((len(Dict[j][i]))) + ",")
                outHeat.write("\n")
            outHeat.close()

            print('......')
            print(".......")
            print("Finished!")


    else:
        print("Skipping HMMER and using output in " + args.out)
        # COVERAGE-BASED ABUNDANCE
        if args.bams != "NA":
            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            BAMmap = open(args.bams)
            normDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            for i in BAMmap:
                string = ''
                ls = i.rstrip().split("\t")
                cell = ls[0]
                for j in ls[1:]:
                    string += " "
                    string += j

                try:
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total / 1000000

                except FileNotFoundError:
                    os.system("mkdir -p %s/contigDepths" % args.out)
                    os.system("jgi_summarize_bam_contig_depths --outputDepth %s/contigDepths/%s.depth%s" % (
                    args.out, cell, string))
                    print("processing... " + cell)
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total / 1000000

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(float(depthDict[cell][contig]))

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            outHeat = open("%s/readDepth.heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
                outHeat.write("\n")

            outHeat.close()

        # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE
        elif args.bam != "NA":
            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

            try:
                total = 0
                depth = open("%s/contigDepths/%s.depth" % (args.out, args.bam))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = float(LS[2])
                        total += float(LS[2])

            except FileNotFoundError:
                os.system("mkdir -p %s/contigDepths" % args.out)
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s/contigDepths/%s.depth %s" % (
                args.out, args.bam, args.bam))
                depth = open("%s/contigDepths/%s.depth" % (args.out, args.bam))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = LS[2]
                        total += float(LS[2])

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(float(depthDict[contig]))

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            outHeat = open("%s/readDepth.heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    outHeat.write(str(SUM(Dict[j][i])) + ",")
                outHeat.write("\n")

            outHeat.close()

        # GENE COUNTS-BASED ABUNDANCE
        else:
            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/summary.csv" % args.out, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "hmm":
                    if not re.match(r'#', i):
                        cell = ls[1]
                        orf = ls[2]
                        contig = allButTheLast(orf, "_")
                        gene = ls[0]
                        Dict[cell][gene].append(gene)

            genes = []
            for i in Dict.keys():
                for j in Dict[i]:
                    if j not in genes:
                        genes.append(j)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in binDirLS:
                if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.binx:

                    if args.orfs:
                        seqdb = args.bins + "/" + i
                        fastaFile = open(args.bins + "/" + i, "r")
                        fastaFile = fasta(fastaFile)
                    else:

                        seqdb = args.bins + "/" + i + "-proteins.faa"
                        fastaFile = open(args.bins + "/" + i + "-proteins.faa", "r")
                        fastaFile = fasta(fastaFile)

                    normDict[i] = len(fastaFile.keys())

            outHeat = open("%s/heatmap.csv" % args.out, "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in sorted(genes):
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):

                    if args.norm:
                        outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
                    else:
                        outHeat.write(str((len(Dict[j][i]))) + ",")
                outHeat.write("\n")
            outHeat.close()

            print('......')
            print(".......")
            print("Finished!")


if __name__ == '__main__':
    main()



