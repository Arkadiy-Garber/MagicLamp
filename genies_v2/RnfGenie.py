#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# RnfGenie: HMM-based identification of the Rnf (Rhodobacter nitrogen
# fixation) electron transport complex, a membrane-bound ferredoxin:NAD+
# oxidoreductase. In many anaerobes (e.g. Acetobacterium woodii, Clostridium
# spp.) Rnf couples electron transfer to Na+ or H+ translocation across the
# membrane, generating an ion gradient that drives ATP synthesis via the
# corresponding F-type or V/A-type ATP synthase.
#
# The complex consists of six subunits encoded by rnfABCDGE. Detection of
# all six in a tight genomic cluster is strong evidence for a functional
# Rnf complex.
#
# Part of the MagicLamp suite (https://github.com/Arkadiy-Garber/MagicLamp).
# Modeled on Lucifer.py / ATPGenie.py / PortNaGenie.py for stylistic
# consistency. All HMMs run with hmmsearch --cut_tc.


def main():

    # =================================================================
    # Helper functions
    # =================================================================

    def cluster(data, maxgap):
        '''Arrange data into groups where successive elements
           differ by no more than *maxgap*'''
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
        prog="MagicLamp.py RnfGenie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r'''
        *******************************************************

        RnfGenie - part of the MagicLamp suite.
        Developed by Arkadiy Garber, Daniel Naphas, and collaborators.
        Please send comments and inquiries to ark@midauthorbio.com

           ____        __ ____            _
          |  _ \ _ __ / _|  _ \ ___ _ __ (_) ___
          | |_) | '_ \ |_| | | / _ \ '_ \| |/ _ \
          |  _ <| | | |  _| |_| |  __/ | | | |  __/
          |_| \_\_| |_|_| |____/ \___|_| |_|_|\___|

              Rnf electron transport complex profiler
                 (rnfABCDGE / Na+- or H+-pumping
                ferredoxin:NAD+ oxidoreductase)

        Single-pass analysis:
          Profiles genomes/proteomes with TIGRFAM HMMs for the six
          Rnf subunits (rnfA, rnfB, rnfC, rnfD, rnfG, rnfE) and reports
          per-genome subunit detection plus genomic-proximity clustering
          to flag putative complete operons.

        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins/genomes/assemblies", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic 'cluster'."
                                             " This number should be an integer and should reflect the maximum number "
                                             "of genes in between putative Rnf-related genes identified by the HMM "
                                             "database (default=5)", default=5)

    parser.add_argument('-out', type=str, help="name output directory (default=rnfgenie_out)",
                        default="rnfgenie_out")

    parser.add_argument('-t', type=int, help="number of threads to use for HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each Rnf subunit to "
                             "be normalized to the number of predicted ORFs in each genome or metagenome.",
                        const=True, nargs="?")

    # =================================================================
    # Locate the HMM library
    # =================================================================

    os.system("echo ${rnf_hmms} > HMMlib.txt")
    file = open("HMMlib.txt")
    HMMrootdir = ""
    for i in file:
        HMMrootdir = i.rstrip()
    print(HMMrootdir)

    def have_lib(root):
        if not root:
            return False
        return root

    if not have_lib(HMMrootdir):
        os.system("which MagicLamp.py > mainDir.txt")
        try:
            file = open("mainDir.txt")
            location = ""
            for i in file:
                location = i.rstrip()
            location = allButTheLast(location, "/")
            HMMrootdir = location + "/hmms/rnf"
        except FileNotFoundError:
            pass

    if not have_lib(HMMrootdir):
        print("RnfGenie could not locate its HMM library. Please run the setup.sh script if\n"
              "you have Conda installed. Otherwise, please run setup-noconda.sh and put MagicLamp.py\n"
              "into your $PATH. The library must contain a 'rnf/' subdirectory with the six rnf HMMs.")
        os.system("rm -f HMMlib.txt mainDir.txt")
        raise SystemExit

    os.system("rm -f HMMlib.txt mainDir.txt")

    HMMdir = HMMrootdir

    args = parser.parse_known_args()[0]

    # =================================================================
    # HMM metadata: acc -> (gene_symbol, subunit_label, description)
    # =================================================================
    HMM_META = {
        "TIGR01943": ("rnfA", "RnfA", "Rnf electron transport complex subunit A (membrane)"),
        "TIGR01944": ("rnfB", "RnfB", "Rnf electron transport complex subunit B (cytoplasmic FeS)"),
        "TIGR01945": ("rnfC", "RnfC", "Rnf electron transport complex subunit C (NADH-binding, FeS)"),
        "TIGR01946": ("rnfD", "RnfD", "Rnf electron transport complex subunit D (membrane, FMN)"),
        "TIGR01947": ("rnfG", "RnfG", "Rnf electron transport complex subunit G (membrane, FMN)"),
        "TIGR01948": ("rnfE", "RnfE", "Rnf electron transport complex subunit E (membrane)"),
    }

    SUBUNIT_ORDER = ["TIGR01943", "TIGR01944", "TIGR01945", "TIGR01946", "TIGR01947", "TIGR01948"]

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
        print("Looks like you did not provide an extension for your genomes/bins or assemblies, so RnfGenie does "
              "not know which files in the provided directory are FASTA files that you would like analyzed.")
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)
        answer = input("Would you like RnfGenie to proceed and potentially overwrite files in this directory? (y/n): ")
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
                                print("Unfortunately, this is a problem for RnfGenie because it uses that character "
                                      "as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                except FileNotFoundError:
                    binFile = open("%s/%s" % (binDir, i), "r")
                    for line in binFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \\|")
                                print("Unfortunately, this is a problem for RnfGenie because it uses that character "
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
    print("Profiling Rnf subunits with TIGRFAM HMMs (hmmsearch --cut_tc)")

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
                        if orf not in HMMdict[i] or bitscore > HMMdict[i][orf]["bit"]:
                            HMMdict[i][orf]["hmm"] = hmm
                            HMMdict[i][orf]["evalue"] = evalue
                            HMMdict[i][orf]["bit"] = bitscore
            print("")

    # =================================================================
    # Initial summary CSV
    # =================================================================
    out = open("%s/summary.csv" % outDirectory, "w")
    out.write("cell,ORF,HMM,gene_symbol,subunit,description,evalue,bitscore\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            hmm = HMMdict[key][j]["hmm"]
            acc = hmm.split(".")[0]
            meta = HMM_META.get(acc, ("unknown", "unknown", "unknown"))
            gene_sym, sub_label, desc = meta
            out.write(",".join([
                key, j, hmm, gene_sym, sub_label, "\"" + desc + "\"",
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
        if ls[0] in ("cell", "RnfGenie"):
            continue
        if len(ls) < 8:
            continue
        cell = ls[0]
        orf = ls[1]
        hmm = ls[2]
        gene_sym = ls[3]
        sub_label = ls[4]
        # description was quoted to preserve commas; reassemble if it spans fields
        if ls[5].startswith("\"") and not ls[5].endswith("\""):
            j = 5
            desc_parts = [ls[5]]
            while j + 1 < len(ls) and not ls[j].endswith("\""):
                j += 1
                desc_parts.append(ls[j])
            desc = ",".join(desc_parts).strip("\"")
            evalue = ls[j + 1]
            bitscore = ls[j + 2]
        else:
            desc = ls[5].strip("\"")
            evalue = ls[6]
            bitscore = ls[7]
        seq = BinDict[cell][orf]

        record = {
            "hmm": hmm, "gene_sym": gene_sym, "sub": sub_label, "desc": desc,
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
    # Cluster ORFs by genomic proximity (putative rnfABCDGE operons)
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
                cluster_rows = []
                cluster_subs = []
                cluster_orfs = []
                for l in RemoveDuplicates(k):
                    orf = j + "_" + str(l)
                    rec = SummaryDict[i][orf]
                    if rec == "EMPTY" or rec.get("hmm", "EMPTY") == "EMPTY":
                        continue
                    cluster_rows.append(rec)
                    cluster_subs.append(rec["sub"])
                    cluster_orfs.append(orf)
                if not cluster_rows:
                    continue
                # Cluster completeness: how many of the 6 Rnf subunits are present
                unique_subs = RemoveDuplicates(cluster_subs)
                completeness = "%d/6" % len(unique_subs)
                for orf, rec in zip(cluster_orfs, cluster_rows):
                    out.write(",".join([
                        i, orf, rec["hmm"], rec["gene_sym"], rec["sub"],
                        "\"" + rec["desc"] + "\"", str(rec["evalue"]), str(rec["bit"]),
                        str(counter), completeness, str(rec["seq"])
                    ]) + "\n")
                out.write(",".join(["#"] * 11) + "\n")
                counter += 1
    out.close()

    # =================================================================
    # Final summary CSV
    # =================================================================
    out = open("%s/rnfgenie-summary.csv" % args.out, "w")
    out.write("file,ORF,gene,gene_symbol,subunit,description,evalue,bit_score,cluster_id,cluster_completeness,seq\n")
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
    #   Rows: one per Rnf subunit (six rows) + one summary row for
    #   maximum-completeness cluster per genome (the best operon found).
    # =================================================================
    print("....")
    print(".....")

    sub_rows = []
    for acc in SUBUNIT_ORDER:
        gene_sym, sub_label, _ = HMM_META[acc]
        sub_rows.append("%s_%s" % (acc, sub_label))
    cats = sub_rows + ["best_cluster_completeness"]

    # Per-genome subunit counts
    Dict = defaultdict(lambda: defaultdict(list))
    # Per-genome best cluster completeness (max number of distinct subunits found in any cluster)
    BestCluster = defaultdict(int)

    final = open("%s/rnfgenie-summary.csv" % args.out, "r")
    cluster_id_subs = defaultdict(lambda: defaultdict(set))   # cell -> cluster_id -> set of subunits
    for line in final:
        if re.match(r'#', line):
            continue
        ls = line.rstrip().split(",")
        if len(ls) < 11:
            continue
        if ls[0] in ("file", "bin", "assembly", "genome"):
            continue
        cell = ls[0]
        hmm_acc = ls[2].split(".")[0]
        meta = HMM_META.get(hmm_acc)
        if not meta:
            continue
        gene_sym, sub_label, _ = meta
        # subunit row
        row_label = "%s_%s" % (hmm_acc, sub_label)
        Dict[cell][row_label].append(hmm_acc)
        # track cluster membership
        try:
            cid = int(ls[8])
            cluster_id_subs[cell][cid].add(sub_label)
        except (ValueError, IndexError):
            pass

    for cell, clu in cluster_id_subs.items():
        for cid, subs in clu.items():
            if len(subs) > BestCluster[cell]:
                BestCluster[cell] = len(subs)

    normDict = defaultdict(lambda: 'EMPTY')
    for i in os.listdir(args.bin_dir):
        if lastItem(i.split(".")) == args.bin_ext:
            try:
                file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                file = fasta(file)
                normDict[i] = len(file.keys())
            except FileNotFoundError:
                normDict[i] = 0

    outHeat = open("%s/rnfgenie.heatmap.csv" % outDirectory, "w")
    outHeat.write("X,")
    for i in sorted(Dict.keys()):
        outHeat.write(i + ",")
    outHeat.write("\n")
    for c in cats:
        outHeat.write(c + ",")
        for j in sorted(Dict.keys()):
            if re.match(r'#', j):
                continue
            if c == "best_cluster_completeness":
                # raw count out of 6 -- not normalized even with --norm,
                # since this is a categorical operon-quality metric
                outHeat.write(str(BestCluster[j]) + ",")
            elif args.norm:
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
    print("Results are written to %s/rnfgenie-summary.csv and %s/rnfgenie.heatmap.csv" % (args.out, args.out))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()