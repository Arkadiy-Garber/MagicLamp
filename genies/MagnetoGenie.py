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


def checkMam(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if memoryDict[hmm] == "magnetosome_formation":
                count += 1
    return count


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
        prog="MagnetoGenie.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************
    
        Developed by Arkadiy Garber;
        University of Southern California, Earth Sciences
        University of Montana, Biological Sciences
        Please send comments and inquiries to arkadiyg@usc.edu
    
    
        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/ 
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins")
    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)")

    parser.add_argument('-out', type=str, help="name output directory (default=magnetogenie_out)",
                        default="magnetogenie_out")

    parser.add_argument('--contigs_source', type=str, help="are the provided contigs from a single organism (single)"
                                                           "or are you providing this program with metagenomic/metatranscriptomic assemblies (meta)? "
                                                           "(default=single)", default="single")
    parser.add_argument('-bams', type=str, help="a tab-delimeted file with two columns: first column has the genome or "
                                                "metagenome file names; second column has the corresponding BAM file "
                                                "(provide full path to the BAM file). BAM files are only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts.", default="NA")

    parser.add_argument('-bam', type=str, help="BAM file. This option is only required if you would like to create "
                                               "a heatmap that summarizes the abundance of a certain gene that is based on "
                                               "read coverage, rather than gene counts. If you have more than one BAM file"
                                               "corresponding to different genomes that you are providing, please use the \'-bams\' "
                                               "argument to provide a tab-delimited file that denotes which BAM file (or files) belongs "
                                               "with which genome", default="NA")

    parser.add_argument('--d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                              "This number should be an integer and should reflect the maximum number of "
                                              "genes in between putative iron-related genes identified by the HMM database "
                                              "(default=10)", default=10)

    parser.add_argument('--cpu', type=int, help="number of threads to allow for hmmsearch (default = 1)", default=1)


    # CHECKING FOR CONDA INSTALL
    os.system("echo ${magneto_hmms}/hmm-meta.txt > HMMlib.txt")
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
        os.system("echo ${magneto_hmms} > HMMlib.txt")
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

    os.system("mkdir -p %s" % args.out)
    count = 0
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    out = open("%s/magnetogenie.csv" % (args.out), "w")
    out.write(
        "bin" + "," + "gene" + "," + "process" + "," + "substrate" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "sequence" + "\n")
    for i in binDirLS:
        HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        if not re.match(r'^\.', i) and lastItem(i.split(".")) == args.bin_ext:

            try:
                testFile = open("%s/%s-proteins.faa" % (binDir, i), "r")
                print("")
                print(".")
                print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))

            except FileNotFoundError:
                print("")
                print(".")
                print("Finding ORFs for " + i)
                if args.contigs_source == "single":
                    os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -q" % (
                    binDir, i, binDir, i, binDir, i))
                elif args.contigs_source == "meta":
                    os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -p meta -q" % (
                    binDir, i, binDir, i, binDir, i))
                else:
                    print("WARNING: you did not specify whether the provided FASTA files are single genomes or "
                          "metagenome/metatranscriptome assemblies. By default, FeGenie is assuming that these are "
                          "single genomes, and running Prodigal accordingly. Just an FYI.")
                    os.system("prodigal -i %s%s -a %s%s-proteins.faa -o %s%s-prodigal.out -q" % (
                    binDir, i, binDir, i, binDir, i))

            # print("")
            # print("analyzing " + i)
            fastaFile = open(args.bin_dir + "/" + i + "-proteins.faa", "r")
            fastaFile = fasta(fastaFile)
            os.system("mkdir " + args.bin_dir + "/" + i + "-HMM")
            count = 0
            for hmm in HMMdirLS:
                if hmm != "hmm-meta.txt":
                    count += 1
                    perc = (count / len(HMMdirLS)) * 100
                    # print(str(perc) + "%")
                    # print("%.2f" % perc + "% done")
                    sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc + 1))
                    sys.stdout.flush()
                    if not re.match(r'^\.', hmm):
                        os.system("hmmsearch --cpu " + str(args.cpu) +
                                  " --tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
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
                    i + "," + HMMdict[key]["hmm"] + "," + HMMdict[key]["process"] + "," + HMMdict[key]["element"] +
                    "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key]["bitscore"] + "," + fastaFile[key]
                    + "\n")
            os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")

    out.close()

    summaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    summary = open("%s/magnetogenie.csv" % (args.out), "r")
    for i in summary:
        ls = i.rstrip().split(",")
        if ls != "bin":
            genome = ls[0]
            gene = ls[1]
            process = ls[2]
            substrate = ls[3]
            orf = ls[4]
            evalue = ls[5]
            bitscore = ls[6]
            sequence = ls[7]
            bitcut = bitDict[ls[1]]["bit"]
            summaryDict[genome][orf]["gene"] = gene
            summaryDict[genome][orf]["process"] = process
            summaryDict[genome][orf]["substrate"] = substrate
            summaryDict[genome][orf]["evalue"] = evalue
            summaryDict[genome][orf]["bitscore"] = bitscore
            summaryDict[genome][orf]["sequence"] = sequence
            summaryDict[genome][orf]["bitcut"] = bitcut

    # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
    print("")
    print("")
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
    out = open("%s/magnetogenie-2.csv" % (args.out), "w")
    for i in CoordDict.keys():
        for j in CoordDict[i]:
            LS = (CoordDict[i][j])
            clusters = (cluster(LS, args.d))
            for k in clusters:
                if len(RemoveDuplicates(k)) == 1:
                    orf = j + "_" + str(k[0])

                    out.write(summaryDict[i][orf]["substrate"] + "," + summaryDict[i][orf]["process"] + "," +
                              summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf]["evalue"] +
                              "," + str(summaryDict[i][orf]["bitscore"]) + "," + str(
                        summaryDict[i][orf]["bitcut"]) + "," +
                              str(summaryDict[i][orf]["sequence"]) + "," + str(counter) + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1

                else:
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)

                        out.write(summaryDict[i][orf]["substrate"] + "," + summaryDict[i][orf]["process"] + "," +
                                  summaryDict[i][orf]["gene"] + "," + i + "," + orf + "," + summaryDict[i][orf][
                                      "evalue"] +
                                  "," + str(summaryDict[i][orf]["bitscore"]) + "," + str(
                            summaryDict[i][orf]["bitcut"]) +
                                  "," + str(summaryDict[i][orf]["sequence"]) + "," + str(counter) + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1

    out.close()

    # ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************
    clusterDict = defaultdict(lambda: defaultdict(list))
    summary = open("%s/magnetogenie-2.csv" % (args.out), "r")
    for i in summary:
        print(i.rstrip())
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            clu = int(ls[9])
            clusterDict[clu]["line"].append(ls)
            clusterDict[clu]["gene"].append(ls[2])

    print("..")
    print("...")
    out = open("%s/magnetogenie-3.csv" % (args.out), "w")
    out.write("file" + "," + "gene" + "," "orf" + "," "eval" + "," "bit_score" + "," "bit_score_cutoff" + "," "AA_seq" + "," + "cluster_id" + "\n")
    for i in sorted(clusterDict.keys()):
        ls = (clusterDict[i]["gene"])


        ############################################################################################################
        if "MamA" in ls or "MamB" in ls or "MamE" in ls or "MamI" in ls or "MamK" in ls or "MamL" in ls or "MamM" in ls or "MamO" \
                in ls or "MamP" in ls or "MamQ" in ls:
            mam = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

            if unique(ls, mam) < 5:
                pass

            else:
                for j in clusterDict[i]["line"]:
                    out.write(j[3] + "," + j[2] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                out.write("####################################################" + "\n")

    out.close()

    os.system("rm %s/magnetogenie.csv" % (args.out))
    os.system("rm %s/magnetogenie-2.csv" % (args.out))
    os.system("mv %s/magnetogenie-3.csv %s/magnetogenie-summary.csv" % (args.out, args.out))

    print("....")
    print(".....")
    print('......')
    print(".......")
    print("Finished!")
    print("results are in %s/magnetogenie-summary.csv" % (args.out))



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

        cats = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/magnetogenie-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[2]
                    gene = ls[1]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[cell][contig]))

        outHeat = open("%s/magnetogenie.readDepth.heatmap.csv" % (args.out), "w")
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

        cats = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/magnetogenie-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[2]
                    gene = ls[1]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[contig]))

        outHeat = open("%s/magnetogenie.readDepth.heatmap.csv" % (args.out), "w")
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
        cats = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/magnetogenie-summary.csv" % (args.out), "r")
        for i in final:
            if not re.match(r'#', i):
                ls = (i.rstrip().split(","))
                if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                    cell = ls[0]
                    orf = ls[2]
                    gene = ls[1]
                    Dict[cell][gene].append(gene)

        print("\n\n")
        normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.bin_dir):
            if lastItem(i.split(".")) == args.bin_ext:
                file = open("%s/%s-proteins.faa" % (args.bin_dir, i), "r")
                file = fasta(file)
                for j in file.keys():
                    print(j)
                normDict[i] = len(file.keys())

        print("\n\n")
        outHeat = open("%s/magnetogenie.heatmap.csv" % (args.out), "w")
        outHeat.write("X" + ',')
        for i in sorted(Dict.keys()):
            outHeat.write(i + ",")
        outHeat.write("\n")

        for i in cats:
            outHeat.write(i + ",")
            for j in sorted(Dict.keys()):
                if not re.match(r'#', j):
                    print(len(Dict[j][i]))
                    print(int(normDict[j]))
                    print("")
                    outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
            outHeat.write("\n")
        outHeat.close()

        print('......')
        print(".......")
        print("Finished!")



if __name__ == '__main__':
    main()
















