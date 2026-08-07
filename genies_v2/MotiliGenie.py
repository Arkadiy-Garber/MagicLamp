#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# MotiliGenie: HMM-based identification and subcategorization of motility & taxis gene profiler
# in genomes / proteomes.
#
# Part of the MagicLamp suite (https://github.com/Arkadiy-Garber/MagicLamp).
# Modeled on PortGenie.py / ATPGenie.py for stylistic consistency.
#
# HMMs are run with hmmsearch --cut_tc (trusted cutoffs built into each HMM).
# Each HMM is hard-coded to one of the subcategories below; the subcategory
# table was compiled from the curated PGL gene survey (Perez-Rodriguez,
# Mukherjee, Mehta, Viney) attached to the MagicLamp publication.


def main():

    # =================================================================
    # Helper functions (mirrored from PortGenie.py for stylistic consistency)
    # =================================================================

    def cluster(data, maxgap):
        # Arrange data into groups where successive elements
        # differ by no more than *maxgap*
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
                    header = header.split(" ")[0]
                    seq = ''
                else:
                    header = i[1:]
                    header = header.split(" ")[0]
                    seq = ''
            else:
                seq += i
        Dict[header] = seq
        return Dict

    def filt(list, items):
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
        ls = filt(ls, [""])
        return ls

    # =================================================================
    # Argparse
    # =================================================================

    parser = argparse.ArgumentParser(
        prog="MagicLamp.py MotiliGenie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r"""
        *******************************************************

        MotiliGenie - part of the MagicLamp suite.
        Developed by Arkadiy Garber, with curated gene categories from
        Ileana Perez-Rodriguez, Abhijit Mukherjee, Maya Mehta, and
        Isabella Viney (Penn Geomicrobiology Laboratory).
        Please send comments and inquiries to ark@midauthorbio.com
                      __  __       _   _ _ _  ____            _
                     |  \/  | ___ | |_(_) (_)/ ___| ___ _ __ (_) ___
                     | |\/| |/ _ \| __| | | | |  _ / _ \ '_ \| |/ _ \
                     | |  | | (_) | |_| | | | |_| |  __/ | | | |  __/
                     |_|  |_|\___/ \__|_|_|_|\____|\___|_| |_|_|\___|
                  Motility & taxis gene profiler

        Single-pass analysis:
          Profiles genomes/proteomes with TIGRFAM HMMs for genes
          involved in flagellar motility, chemotaxis, thermotaxis,
          and pili-based motility.

        Each hit is assigned to a functional subcategory defined in
        the curated gene survey, and counts are aggregated per
        subcategory in the heatmap output.

        *******************************************************
        """))

    parser.add_argument('-bin_dir', type=str, help="directory of bins/genomes/assemblies", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic 'cluster'."
                                             " This number should be an integer and should reflect the maximum number "
                                             "of genes in between putative motility & taxis gene profiler identified by the HMM "
                                             "database (default=5)", default=5)

    parser.add_argument('-out', type=str, help="name output directory (default=motiligenie_out)",
                        default="motiligenie_out")

    parser.add_argument('-t', type=int, help="number of threads to use for HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each subcategory to "
                             "be normalized to the number of predicted ORFs in each genome or metagenome.",
                        const=True, nargs="?")

    # =================================================================
    # Locate the HMM library
    #   Layout: $$motili_hmms points at a directory containing the .hmm
    #   files directly (flat layout, matching PortGenie convention).
    # =================================================================

    os.system("echo ${motili_hmms} > HMMlib.txt")
    file = open("HMMlib.txt")
    HMMrootdir = ""
    for i in file:
        HMMrootdir = i.rstrip()

    def have_lib(root):
        if not root:
            return False
        return os.path.isdir(root) and any(f.endswith(".hmm") for f in os.listdir(root))

    if not have_lib(HMMrootdir):
        os.system("which MagicLamp.py > mainDir.txt")
        try:
            file = open("mainDir.txt")
            location = ""
            for i in file:
                location = i.rstrip()
            location = allButTheLast(location, "/")
            HMMrootdir = location + "/hmms/motili_hmms"
        except FileNotFoundError:
            pass

    if not have_lib(HMMrootdir):
        print("MotiliGenie could not locate its HMM library. Please run the setup.sh script if\n"
              "you have Conda installed. Otherwise, please run setup-noconda.sh and put MagicLamp.py\n"
              "into your $PATH. Set the $motili_hmms environment variable to point at a directory\n"
              "containing the .hmm files for this genie.")
        os.system("rm -f HMMlib.txt mainDir.txt")
        raise SystemExit

    os.system("rm -f HMMlib.txt mainDir.txt")

    HMMdir = HMMrootdir

    args = parser.parse_known_args()[0]

    # =================================================================
    # HMM metadata
    #   acc (e.g. TIGR01181) -> (subcategory, gene_symbol, description)
    #
    # Subcategories defined in the curated PGL gene survey.
    # =================================================================
    HMM_META = {
        "TIGR00205": ('flagella', 'fliE', 'fliE'),
        "TIGR00206": ('flagella', 'fliF', 'fliF'),
        "TIGR00207": ('flagella', 'fliG', 'fliG'),
        "TIGR00208": ('flagella', 'fliS', 'fliS'),
        "TIGR00229": ('chemotaxis', 'cheBR', 'cheBR'),
        "TIGR00261": ('pili', 'traB', 'traB'),
        "TIGR00328": ('flagella', 'flhB', 'flhB'),
        "TIGR00575": ('miscellaneous', 'ligA1', 'ligA1'),
        "TIGR01103": ('flagella', 'fliP', 'fliP'),
        "TIGR01175": ('pili', 'pilM', 'pilM'),
        "TIGR01395": ('flagella', 'flgC', 'flgC'),
        "TIGR01396": ('flagella', 'flgB', 'flgB'),
        "TIGR01397": ('flagella', 'fliM', 'fliM'),
        "TIGR01398": ('flagella', 'flhA', 'flhA'),
        "TIGR01400": ('flagella', 'fliR', 'fliR'),
        "TIGR01402": ('flagella', 'fliQ', 'fliQ'),
        "TIGR01420": ('pili', 'pilT', 'pilT'),
        "TIGR01707": ('miscellaneous', 'gspI', 'gspI'),
        "TIGR01708": ('miscellaneous', 'gspH', 'gspH'),
        "TIGR01710": ('miscellaneous', 'gspG', 'gspG'),
        "TIGR01818": ('miscellaneous', 'glnG', 'glnG'),
        "TIGR02120": ('miscellaneous', 'gspF', 'gspF'),
        "TIGR02154": ('chemotaxis', 'phoB', 'phoB'),
        "TIGR02350": ('thermotaxis', 'dnaK', 'dnaK'),
        "TIGR02473": ('flagella', 'fliJ', 'fliJ'),
        "TIGR02479": ('flagella', 'fliA', 'fliA'),
        "TIGR02480": ('flagella', 'fliNY', 'fliNY'),
        "TIGR02492": ('flagella', 'flgK', 'flgK'),
        "TIGR02515": ('pili', 'pilQ', 'pilQ'),
        "TIGR02517": ('miscellaneous', 'gspD', 'gspD'),
        "TIGR02520": ('pili', 'pilN', 'pilN'),
        "TIGR02522": ('pili', 'cpaD', 'cpaD'),
        "TIGR02532": ('pili', 'pilA', 'pilA'),
        "TIGR02533": ('flagella', 'pilB1', 'pilB1'),
        "TIGR02537": ('flagella', 'flaG', 'flaG'),
        "TIGR02538": ('pili', 'pilB', 'pilB'),
        "TIGR02541": ('flagella', 'flgJ', 'flgJ'),
        "TIGR02550": ('flagella', 'flgL', 'flgL'),
        "TIGR02742": ('pili', 'trbC', 'trbC'),
        "TIGR02743": ('pili', 'traW', 'traW'),
        "TIGR02744": ('pili', 'trbI', 'trbI'),
        "TIGR02746": ('pili', 'traC', 'traC'),
        "TIGR02747": ('pili', 'traV', 'traV'),
        "TIGR02756": ('pili', 'traK', 'traK'),
        "TIGR02761": ('pili', 'traE', 'traE'),
        "TIGR02762": ('pili', 'traL', 'traL'),
        "TIGR02788": ('miscellaneous', 'virB11', 'virB11'),
        "TIGR02937": ('chemotaxis', 'sigD', 'sigD'),
        "TIGR03021": ('pili', 'pilP', 'pilP'),
        "TIGR03170": ('flagella', 'flgA', 'flgA'),
        "TIGR03177": ('pili', 'cpaB', 'cpaB'),
        "TIGR03300": ('miscellaneous', 'bamB', 'bamB'),
        "TIGR03310": ('miscellaneous', 'mocA', 'mocA'),
        "TIGR03496": ('flagella', 'fliI', 'fliI'),
        "TIGR03499": ('flagella', 'flhF', 'flhF'),
        "TIGR03500": ('flagella', 'fliOZ', 'fliOZ'),
        "TIGR03506": ('flagella', 'FlgEFG', 'FlgEFG'),
        "TIGR03787": ('miscellaneous', 'chvl', 'chvl'),
        "TIGR03818": ('flagella', 'motA', 'motA'),
        "TIGR03824": ('flagella', 'flgM', 'flgM'),
        "TIGR03825": ('flagella', 'fliH', 'fliH'),
    }

    SUBCATEGORIES = [
        'chemotaxis',
        'flagella',
        'miscellaneous',
        'pili',
        'thermotaxis',
    ]

    # =================================================================
    # Argument validation
    # =================================================================

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

    if args.bin_ext != "NA":
        print(".")
    else:
        print("Looks like you did not provide an extension for your genomes/bins or assemblies, so MotiliGenie does "
              "not know which files in the provided directory are FASTA files that you would like analyzed.")
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)
        # answer = input("Would you like MotiliGenie to proceed and potentially overwrite files in this directory? (y/n): ")
        answer = "y"
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
        outDirectory = "%s" % args.out[0:len(args.out) - 1]
    else:
        outDirectory = "%s" % args.out

    print("All required arguments provided!")
    print("")

    # =================================================================
    # ORF calling (Prodigal) or GenBank protein extraction
    # =================================================================

    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext:
            cell = i
            if not args.gbk:
                try:
                    testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    print("ORFs for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                    for line in testFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \\|")
                                print("Unfortunately, this is a problem for MotiliGenie because it uses that character "
                                      "as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                except FileNotFoundError:
                    binFile = open("%s/%s" % (binDir, i), "r")
                    for line in binFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \\|")
                                print("Unfortunately, this is a problem for MotiliGenie because it uses that character "
                                      "as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                    print("Finding ORFs for " + cell)
                    if args.meta:
                        os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa "
                                  "-o %s/ORF_calls/%s-prodigal.out -p meta -q" % (
                                      binDir, i, outDirectory, i, outDirectory, i))
                    else:
                        os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa "
                                  "-o %s/ORF_calls/%s-prodigal.out -q" % (
                                      binDir, i, outDirectory, i, outDirectory, i))
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
                            locus = ls.split("       ")[1].split(" ")[0]
                        if re.findall(r'/locus_tag', ls):
                            locusTag = ls.split("=")[1]
                            locusTag = remove(locusTag, ["\""])
                            counter += 1
                        if counter > 0:
                            gbkDict[locus].append(locusTag)
                            counter = 0
                else:
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = ls.split("       ")[1].split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = ls.split("            ")[1]
                            start = gene.split("..")[0]
                            end = gene.split("..")[1]
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)
                            counter += 1
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

            file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
            file = fasta(file)
            for j in file.keys():
                orf = j.split(" # ")[0]
                BinDict[cell][orf] = file[j]

    # =================================================================
    # HMM profiling (single pass with --cut_tc)
    # =================================================================
    print("")
    print("starting main pipeline...")
    print("Profiling motility & taxis gene profiler with TIGRFAM HMMs (hmmsearch --cut_tc)")

    HMMdirLS = [h for h in os.listdir(HMMdir) if h.endswith(".hmm")]

    HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))

    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext:
            os.system("mkdir -p " + outDirectory + "/" + i + "-HMM")
            count = 0
            for hmm in HMMdirLS:
                count += 1
                perc = (count / len(HMMdirLS)) * 100
                sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                sys.stdout.flush()

                os.system(
                    "hmmsearch --cpu %d --cut_tc --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                    % (int(args.t), outDirectory, i, hmm, outDirectory, i, hmm, HMMdir, hmm,
                       outDirectory, i)
                )
                os.system("rm " + outDirectory + "/" + i + "-HMM/" + hmm + ".txt")

                hmmout = open(outDirectory + "/" + i + "-HMM/" + hmm + ".tblout", "r")
                for line in hmmout:
                    if not re.match(r'#', line):
                        ls = delim(line)
                        evalue = float(ls[4])
                        bitscore = float(ls[5])
                        orf = ls[0]
                        # Track the best-scoring HMM per ORF
                        if orf not in HMMdict[i] or bitscore > HMMdict[i][orf]["bit"]:
                            HMMdict[i][orf]["hmm"] = hmm
                            HMMdict[i][orf]["evalue"] = evalue
                            HMMdict[i][orf]["bit"] = bitscore
            print("")

    # =================================================================
    # Initial summary CSV
    # =================================================================
    out = open("%s/summary.csv" % outDirectory, "w")
    out.write("cell,ORF,HMM,gene_symbol,description,subcategory,evalue,bitscore\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            hmm = HMMdict[key][j]["hmm"]
            acc = hmm.split(".")[0]
            meta = HMM_META.get(acc, ("unknown", "unknown", "unknown"))
            sub, gene_sym, desc = meta
            out.write(",".join([
                key, j, hmm, gene_sym, "\"" + desc + "\"", sub,
                str(HMMdict[key][j]["evalue"]), str(HMMdict[key][j]["bit"])
            ]) + "\n")
    out.close()

    # =================================================================
    # Dereplicate (keep best-bitscore HMM per ORF)
    # =================================================================
    summary = open(outDirectory + "/summary.csv", "r")
    SummaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for line in summary:
        ls = line.rstrip().split(",")
        if ls[0] in ("cell", "MotiliGenie"):
            continue
        if len(ls) < 8:
            continue
        cell = ls[0]
        orf = ls[1]
        hmm = ls[2]
        gene_sym = ls[3]
        # description was quoted to handle commas; rejoin if needed
        if ls[4].startswith("\"") and not ls[4].endswith("\""):
            j = 4
            desc_parts = [ls[4]]
            while j + 1 < len(ls) and not ls[j].endswith("\""):
                j += 1
                desc_parts.append(ls[j])
            desc = ",".join(desc_parts).strip("\"")
            sub = ls[j + 1]
            evalue = ls[j + 2]
            bitscore = ls[j + 3]
        else:
            desc = ls[4].strip("\"")
            sub = ls[5]
            evalue = ls[6]
            bitscore = ls[7]
        seq = BinDict[cell][orf]

        record = {
            "hmm": hmm, "gene_sym": gene_sym, "desc": desc, "sub": sub,
            "evalue": evalue, "bit": bitscore, "seq": seq,
        }

        if orf not in SummaryDict[cell]:
            SummaryDict[cell][orf] = record
        else:
            try:
                if float(bitscore) > float(SummaryDict[cell][orf]["bit"]):
                    SummaryDict[cell][orf] = record
            except (ValueError, TypeError):
                pass

    # =================================================================
    # Cluster ORFs by genomic proximity (putative operons)
    # =================================================================
    print("Identifying genomic proximities and putative operons")
    CoordDict = defaultdict(lambda: defaultdict(list))
    for i in SummaryDict.keys():
        if i != "cell":
            for j in SummaryDict[i]:
                contig = allButTheLast(j, "_")
                numOrf = lastItem(j.split("_"))
                try:
                    CoordDict[i][contig].append(int(numOrf))
                except ValueError:
                    pass

    counter = 0
    print("Clustering ORFs...")
    print("")
    out = open(outDirectory + "/summary-2.csv", "w")
    for i in CoordDict.keys():
        print(".")
        for j in CoordDict[i]:
            LS = CoordDict[i][j]
            if not LS:
                continue
            clusters = cluster(LS, args.d)
            for k in clusters:
                for l in RemoveDuplicates(k):
                    orf = j + "_" + str(l)
                    rec = SummaryDict[i][orf]
                    if rec == "EMPTY" or rec.get("hmm", "EMPTY") == "EMPTY":
                        continue
                    out.write(",".join([
                        i, orf, rec["hmm"], rec["gene_sym"], "\"" + rec["desc"] + "\"",
                        rec["sub"], str(rec["evalue"]), str(rec["bit"]),
                        str(counter), str(rec["seq"])
                    ]) + "\n")
                out.write(",".join(["#"] * 10) + "\n")
                counter += 1
    out.close()

    # =================================================================
    # Final summary CSV
    # =================================================================
    out = open("%s/motiligenie-summary.csv" % args.out, "w")
    out.write("file,ORF,gene,gene_symbol,description,subcategory,evalue,bit_score,cluster_id,seq\n")
    summary = open("%s/summary-2.csv" % args.out, "r")
    for line in summary:
        if not re.match(r'#', line):
            out.write(line)
        else:
            out.write("####################################################\n")
    out.close()

    os.system("rm %s/summary.csv" % args.out)
    os.system("rm %s/summary-2.csv" % args.out)

    os.system("mkdir -p %s/HMM_results" % outDirectory)
    os.system("rm -f %s/ORF_calls/*-prodigal.out" % outDirectory)
    os.system("rm -rf %s/HMM_results/*-HMM" % outDirectory)
    os.system("mv -f %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

    # =================================================================
    # Heatmap-compatible CSV
    #   Two row blocks:
    #     1. Per-HMM detection rows (one row per HMM in the library)
    #     2. Per-subcategory aggregate rows (TOTAL_<subcategory>)
    # =================================================================
    print("....")
    print(".....")

    hmm_rows = []
    for acc in sorted(HMM_META.keys()):
        sub, gene_sym, _ = HMM_META[acc]
        # Sanitize subcategory for use in a row label: strip spaces, slashes
        sub_clean = sub.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
        hmm_rows.append("%s_%s_%s" % (acc, gene_sym, sub_clean))

    sub_rows = []
    for s in SUBCATEGORIES:
        sub_clean = s.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
        sub_rows.append("TOTAL_" + sub_clean)

    cats = hmm_rows + sub_rows

    Dict = defaultdict(lambda: defaultdict(list))
    final = open("%s/motiligenie-summary.csv" % args.out, "r")
    for line in final:
        if re.match(r'#', line):
            continue
        ls = line.rstrip().split(",")
        if len(ls) < 9:
            continue
        if ls[0] in ("file", "bin", "assembly", "genome"):
            continue
        cell = ls[0]
        hmm_acc = ls[2].split(".")[0]
        meta = HMM_META.get(hmm_acc)
        if not meta:
            continue
        sub, gene_sym, _ = meta
        sub_clean = sub.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")

        # Per-HMM row
        hmm_row = "%s_%s_%s" % (hmm_acc, gene_sym, sub_clean)
        Dict[cell][hmm_row].append(hmm_acc)
        # Subcategory aggregate row
        Dict[cell]["TOTAL_" + sub_clean].append(hmm_acc)

    normDict = defaultdict(lambda: 'EMPTY')
    for i in os.listdir(args.bin_dir):
        if lastItem(i.split(".")) == args.bin_ext:
            try:
                file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                file = fasta(file)
                normDict[i] = len(file.keys())
            except FileNotFoundError:
                normDict[i] = 0

    outHeat = open("%s/motiligenie.heatmap.csv" % outDirectory, "w")
    outHeat.write("X,")
    for i in sorted(Dict.keys()):
        outHeat.write(i + ",")
    outHeat.write("\n")
    for c in cats:
        outHeat.write(c + ",")
        for j in sorted(Dict.keys()):
            if re.match(r'#', j):
                continue
            if args.norm:
                if normDict[j] and int(normDict[j]) > 0:
                    outHeat.write(str((len(Dict[j][c]) / int(normDict[j])) * float(100)) + ",")
                else:
                    outHeat.write("0,")
            else:
                outHeat.write(str(len(Dict[j][c])) + ",")
        outHeat.write("\n")
    outHeat.close()

    print("......")
    print(".......")
    print("Finished!")
    print("")
    print("Results are written to %s/motiligenie-summary.csv and %s/motiligenie.heatmap.csv" % (args.out, args.out))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()
