#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


def main():
    def SUM(ls):
        count = 0
        for i in ls:
            count += float(i)
        return count

    def unique(ls, ls2):
        unqlist = []
        for i in ls:
            if i not in unqlist and i in ls2:
                unqlist.append(i)
        return len(unqlist)

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

    def remove(stringOrlist, list):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                pass
        outString = "".join(emptyList)
        return outString

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

    parser = argparse.ArgumentParser(
        prog="Genie Template",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************

        Developed by Arkadiy Garber;
        Middle Author Bioinformatics LLC
        Please send comments and inquiries to ark@midauthorbio.com
                                  .-=-.
                             /  ! ) )
                          __ \__/__/
                         / _<( ^.^ )   Your wish is my command...
                        / /   \ c /O
                        \ \_.-./=\.-._     _
                         `-._  `~`    `-,./_<
                             `\' \'\`'----'
                           *   \  . \          *
                                `-~~~\   .
                           .      `-._`-._   *
                                 *    `~~~-,      *
                       ()                   * )
                      <^^>             *     (   .
                     .-""-.                    )
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

    parser.add_argument('-genie', type=str, help="name of genie", default="NA")

    parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=5)", default=5)

    parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format",
                        default="NA")

    parser.add_argument('-out', type=str, help="name output directory (default=genie_out)",
                        default="genie_out")

    parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                        default=1000)

    parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each iron gene category to be normalized to "
                             "the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, this genie will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, this genie will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--makeplots', type=str,
                        help="include this flag if you would like this genie to make some figures from your data?. "
                             "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Please be sure to install all the required R packages as instrcuted in the MagicLamp repo home: "
                             "https://github.com/Arkadiy-Garber/MagicLamp. "
                             "If you see error or warning messages associated with Rscript, you can still expect to "
                             "see the main output (CSV files) from this genie.", const=True, nargs="?")

    args = parser.parse_known_args()[0]
    # CHECKING FOR CONDA INSTALL
    genie = args.genie
    os.system("echo ${%s_hmms} > HMMlib.txt" % genie)
    os.system("echo ${rscripts} > rscripts.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = i.rstrip()

    rscriptDir = "~/MagicLamp/rscripts/"
    try:
        file = open("rscripts.txt")
        for i in file:
            rscriptDir = i.rstrip()

    except FileNotFoundError:
        os.system("which MagicLamp.py > mainDir.txt")

        file = open("mainDir.txt")
        location = "~/MagicLamp/MagicLamp.py"
        for i in file:
            location = i.rstrip()
        location = allButTheLast(location, "/")

        HMMdir = location + "/hmms/%s/" % genie
        rscriptDir = location + "/rscripts/"

        try:
            file = open("rscripts.txt")
            for i in file:
                rscriptDir = i.rstrip()
        except FileNotFoundError:
            print("MagicLamp script could not locate the required directories. Please run the setup.sh script if \n"
                  "you have Conda installed. Otherwise, please run the setupe-noconda.sh script and put MagicLamp.py \n"
                  "into your $PATH")
            raise SystemExit

    os.system("rm -rf HMMlib.txt rscripts.txt mainDir.txt")

    # ************** Checking for the required arguments ******************* #
    cwd = os.getcwd()
    print("checking arguments")
    binDir = args.bin_dir + "/"
    binDirLS = os.listdir(args.bin_dir)
    os.system("mkdir -p %s" % args.out)
    os.system("mkdir -p %s/ORF_calls" % args.out)

    if lastItem(args.out) == "/":
        outDirectory = "%s" % args.out[0:len(args.out) - 1]
        outDirectoryLS = os.listdir(outDirectory)
    else:
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)

    # *************** CALL ORFS FROM BINS AND READ THE ORFS INTO HASH MEMORY ************************ #
    BinDict = defaultdict(lambda: defaultdict(lambda: ''))
    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext:
            cell = i

            # os.system("GB-or-FA.py %s/%s type.txt" % (binDir, i))
            # file = open("type.txt")
            fileType = "fa"
            # for j in file:
            #     fileType = j.rstrip()
            # os.system("rm type.txt")

            if fileType == "fa":

                try:
                    testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))

                except FileNotFoundError:
                    binFile = open("%s/%s" % (binDir, i), "r")
                    print("Finding ORFs for " + i)
                    if args.meta:
                        os.system(
                            "prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -p meta -q" % (
                                binDir, i, outDirectory, i, outDirectory, i))
                    else:
                        os.system(
                            "prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -q" % (
                                binDir, i, outDirectory, i, outDirectory, i))
            elif fileType == "gbk":
                os.system('gb2faa.py %s/%s %s/ORF_calls/%s-proteins.faa type.txt --protein-only' % (binDir, i, outDirectory, i))

                file = open("type.txt")
                fileType = "contigs"
                for j in file:
                    fileType = j.rstrip()
                os.system("rm type.txt")

                if fileType == "proteins":
                    pass

                else:
                    if args.meta:
                        os.system(
                            "prodigal -i %s/%s-proteins.faa -a %s/ORF_calls/%s-proteins2.faa -o %s/ORF_calls/%s-prodigal.out -p meta -q" % (
                                binDir, i, outDirectory, i, outDirectory, i))
                    else:
                        os.system(
                            "prodigal -i %s/%s-proteins.faa -a %s/ORF_calls/%s-proteins2.faa -o %s/ORF_calls/%s-prodigal.out -q" % (
                                binDir, i, outDirectory, i, outDirectory, i))

                    os.system("mv %s/ORF_calls/%s-proteins2.faa %s/ORF_calls/%s-proteins.faa" % (outDirectory, i, outDirectory, i))

            else:
                print("File type not recognized. Please provide a GenBank or FASTA file")
                raise SystemExit

            file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
            file = fasta(file)
            for j in file.keys():
                orf = j.split(" ")[0]
                BinDict[cell][orf] = str(file[j])

    # ******************* BEGINNING MAIN ALGORITHM **********************************))))
    print("starting main pipeline...")
    # HMMdirLS = os.listdir(HMMdir)
    # HMMdir = HMMdir + "/"
    HMMdirLS = [f for f in os.listdir(HMMdir) if f.lower().endswith('.hmm')]
    HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
        if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
            os.system(
                "mkdir -p " + outDirectory + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

            count = 0
            for hmm in HMMdirLS:  # ITERATING THROUGH ALL THE HMM FILES IN THE HMM DIRECTORY
                count += 1
                perc = (count / len(HMMdirLS)) * 100
                sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                sys.stdout.flush()

                os.system(
                    "hmmsearch --cpu %d --cut_tc --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                    % (int(args.t), outDirectory, i, hmm, outDirectory, i, hmm, HMMdir, hmm,
                       outDirectory, i)
                )

                # REMOVING THE STANDARD OUTPUT FILE
                os.system(
                    "rm " + outDirectory + "/" + i + "-HMM/" + hmm + ".txt"
                )

                # READING IN THE HMMSEARCH RESULTS (TBLOUT) OUT FILE
                hmmout = open(outDirectory + "/" + i + "-HMM/" + hmm + ".tblout", "r")

                # COLLECTING SIGNIFICANT HMM HITS IN THE FILE
                for line in hmmout:
                    if not re.match(r'#', line):
                        ls = delim(line)
                        evalue = float(ls[4])
                        bit = float(ls[5])
                        orf = ls[0]
                        # if evalue < float(1E-1):  # FILTERING OUT BACKGROUND NOISE
                            # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS

                        if orf not in HMMdict[i]:
                            # HMMdict[i][orf]["hmm"] = hmm.split(".hm")[0]
                            HMMdict[i][orf]["hmm"] = ls[2]
                            HMMdict[i][orf]["evalue"] = evalue
                            HMMdict[i][orf]["bit"] = bit
                        else:
                            # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                            if bit > HMMdict[i][orf]["bit"]:
                                HMMdict[i][orf]["hmm"] = ls[2]
                                HMMdict[i][orf]["evalue"] = evalue
                                HMMdict[i][orf]["bit"] = bit

            print("")

    out = open("%s/summary.csv" % (outDirectory), "w")
    out.write("organism" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            out.write(key + "," + j + "," + str(HMMdict[key][j]["hmm"]) + "," +
                      str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "\n")

    out.close()
    # ****************************************** DEREPLICATION *********************************************************
    summary = open(outDirectory + "/summary.csv", "r")
    SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] != "organism":
            if len(ls) > 0:
                cell = ls[0]
                orf = str(ls[1])
                hmm = ls[2]
                evalue = ls[3]
                hmmBit = ls[4]
                seq = BinDict[cell][orf]

                if cell not in SummaryDict.keys():
                    SummaryDict[cell][orf]["hmm"] = hmm
                    SummaryDict[cell][orf]["hmmBit"] = hmmBit
                    SummaryDict[cell][orf]["e"] = evalue
                    SummaryDict[cell][orf]["seq"] = seq

                else:
                    if orf not in SummaryDict[cell]:
                        SummaryDict[cell][orf]["hmm"] = hmm
                        SummaryDict[cell][orf]["hmmBit"] = hmmBit
                        SummaryDict[cell][orf]["e"] = evalue
                        SummaryDict[cell][orf]["seq"] = seq

                    else:
                        if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                            SummaryDict[cell][orf]["hmm"] = hmm
                            SummaryDict[cell][orf]["hmmBit"] = hmmBit
                            SummaryDict[cell][orf]["e"] = evalue
                            SummaryDict[cell][orf]["seq"] = seq

    # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
    print("Identifying genomic proximities and putative operons")
    CoordDict = defaultdict(lambda: defaultdict(list))
    # locusDict = defaultdict(lambda: '-')
    numDict = defaultdict(lambda: defaultdict(lambda: '-'))
    for i in SummaryDict.keys():
        if i != "organism":
            for j in SummaryDict[i]:
                # print(i + "\t" + j + "\t" + str(SummaryDict[i][j]["hmm"]) + "\t" + str(SummaryDict[i][j]["e"]) + "\t" + str(SummaryDict[i][j]["hmmBit"]))
                # print(j)
                # print(j.split(";"))
                contig = allButTheLast(j, "_")
                numOrf = lastItem(j.split("_"))
                # # locusDict[j] = contig + "_" + str(numOrf)
                # locusDict[contig + "_" + str(numOrf)] = j.split("_")[1]
                CoordDict[i][contig].append(int(numOrf))
                numDict[contig][int(numOrf)] = j

    counter = 0
    print("Clustering ORFs...")
    out = open(outDirectory + "/summary-2.csv", "w")
    out.write("organism" + "," + "locus" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "," + "clusterID" + "," + "ORF_sequence\n")
    for i in CoordDict.keys():
        print(".")
        for j in CoordDict[i]:
            LS = (CoordDict[i][j])
            clusters = (cluster(LS, args.d))
            for k in clusters:
                if len(RemoveDuplicates(k)) == 1:
                    # orf = j + "_" + str(k[0])
                    locus = numDict[j][k[0]]


                    out.write(
                        i + "," + locus + "," + SummaryDict[i][locus]["hmm"] + "," + SummaryDict[i][locus]["e"] + "," + str(
                            SummaryDict[i][locus]["hmmBit"]) + "," +
                        str(counter) + "," + str(SummaryDict[i][locus]["seq"]) + "\n")
                    out.write("################\n")
                    counter += 1

                else:
                    for l in RemoveDuplicates(k):
                        # orf = j + "_" + str(l)
                        locus = numDict[j][l]

                        out.write(i + "," + locus + "," + SummaryDict[i][locus]["hmm"] + "," + SummaryDict[i][locus][
                            "e"] + "," + str(SummaryDict[i][locus]["hmmBit"]) +
                                  "," + str(counter) + "," + str(SummaryDict[i][locus]["seq"]) + "\n")
                    out.write("################\n")
                    counter += 1
    out.close()

    os.system("rm %s/summary.csv" % (args.out))
    os.system("mv %s/summary-2.csv %s/%sgenie-summary.csv" % (args.out, args.out, genie))

    os.system("mkdir -p %s/HMM_results" % outDirectory)
    os.system("rm -rf %s/HMM_results/*-HMM" % outDirectory)
    os.system("mv -f %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

    # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    # GENE COUNTS-BASED ABUNDANCE
    cats = []
    Dict = defaultdict(lambda: defaultdict(list))
    final = open("%s/%sgenie-summary.csv" % (args.out, genie), "r")
    for i in final:
        if not re.match(r'#', i):
            ls = (i.rstrip().split(","))
            if ls[0] != "organism":
                cell = ls[0]
                gene = ls[2]
                if gene not in cats:
                    cats.append(gene)
                Dict[cell][gene].append(gene)

    normDict = defaultdict(lambda: 0)
    for i in os.listdir(args.bin_dir):
        if lastItem(i.split(".")) == args.bin_ext:
            file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
            file = fasta(file)
            normDict[i] = len(file.keys())

    outHeat = open("%s/%sgenie.heatmap.csv" % (outDirectory, genie), "w")
    outHeat.write("X" + ',')
    for i in sorted(Dict.keys()):
        outHeat.write(i + ",")
    outHeat.write("\n")

    for i in cats:
        outHeat.write(i + ",")
        for j in sorted(Dict.keys()):
            if not re.match(r'#', j):
                if args.norm:
                    outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
                else:
                    outHeat.write(str((len(Dict[j][i]))) + ",")
        outHeat.write("\n")
    outHeat.close()

    print("Finished!")

    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots:
        print("Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
              "This will not affect any of the output data that was already created. If you see plots generated, great! "
              "If not, you can plot the data as you wish on your own, or start an issue on MagicLamp's GitHub repository\n")

        if args.norm:
            os.system(
                "Rscript --vanilla %s/DotPlot.R %s/%sgenie.heatmap.csv %s/" % (rscriptDir, outDirectory, genie, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))
        else:
            os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))

    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots:
        print("Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
              "This will not affect any of the output data that was already created. If you see plots generated, great! "
              "If not, you can plot the data as you wish on your own, or start an issue on MagicLamp's GitHub repository\n")

        if args.norm:
            os.system("Rscript --vanilla %s/DotPlot.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))
        else:
            os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%sgenie.heatmap.csv %s/" % (
            rscriptDir, outDirectory, genie, outDirectory))

    print("Results are written to %s/%sgenie-summary.csv and %s/%sgenie.heatmap.csv" % (args.out, genie, args.out, genie))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()

