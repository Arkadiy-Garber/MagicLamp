#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# ATPGenie: HMM-based identification and categorization of ATP synthase / ATPase
# subunits (F-type and V/A-type), plus classification of the membrane proteolipid
# rotor subunit (atpE / K) as Na+ or H+ pumping using custom 2TM-hairpin HMMs.
#
# Part of the MagicLamp suite (https://github.com/Arkadiy-Garber/MagicLamp).
# Modeled on Lucifer.py (Garber & Merino, USC).
#
# Subunit HMMs (TIGRFAM + Pfam) are run with hmmsearch --cut_tc, using the
# trusted cutoffs built into each HMM. The four ion-preference HMMs
# (F-type_Na_pumping, F-type_H_pumping, V-type_Na_binding, V-type_H_pumping)
# are run with no bitscore cutoff; the best-scoring of the four per
# proteolipid ORF determines the ion call.


def main():

    # =================================================================
    # Helper functions (mirrored from Lucifer.py for stylistic consistency)
    # =================================================================

    def unique(ls, ls2):
        unqlist = []
        for i in ls:
            if i not in unqlist and i in ls2:
                unqlist.append(i)
        return len(unqlist)

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
        prog="MagicLamp.py ATPGenie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r'''
        *******************************************************

        ATPGenie - part of the MagicLamp suite.
        Developed by Arkadiy Garber, Paul Carini.
        Please send comments and inquiries to ark@midauthorbio.com

                    _  _____ ____   ____            _
             /\    |_||_   _|  _ \ / ___| ___ _ __ (_) ___
            /  \    _   | | | |_) | |  _ / _ \ '_ \| |/ _ \
           / /\ \  | |  | | |  __/| |_| |  __/ | | | |  __/
          / ____ \ | |  | | | |    \____|\___|_| |_|_|\___|
         /_/    \_\|_|  |_| |_|

                       ATP synthase / ATPase profiler
                          + Na+/H+ ion classifier

        Two-stage analysis:
          1. Profile genomes/proteomes with TIGRFAM + Pfam HMMs
             (run with hmmsearch --cut_tc) to identify F-type and
             V/A-type ATPase subunits.
          2. Classify the membrane proteolipid rotor (atpE / K) as
             Na+ or H+ pumping using four custom 2TM-hairpin HMMs
             (best-of-four wins; no bitscore cutoff).

        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins/genomes/assemblies", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic 'cluster'."
                                             " This number should be an integer and should reflect the maximum number "
                                             "of genes in between putative ATPase-related genes identified by the HMM "
                                             "database (default=5)", default=5)

    parser.add_argument('-out', type=str, help="name output directory (default=atpgenie_out)",
                        default="atpgenie_out")

    parser.add_argument('-t', type=int, help="number of threads to use for HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each ATPase subunit category to "
                             "be normalized to the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, ATPGenie will create a heatmap-compatible CSV output with raw gene "
                             "counts. With normalization, ATPGenie will create a heatmap-compatible CSV with "
                             "'normalized gene abundances'.", const=True, nargs="?")

    # =================================================================
    # Locate the HMM library (mirrors Lucifer's discovery logic)
    # =================================================================

    os.system("echo ${atp_hmms} > HMMlib.txt")
    file = open("HMMlib.txt")
    HMMrootdir = ""
    for i in file:
        HMMrootdir = i.rstrip()

    def have_lib(root):
        if not root:
            return False
        for sub in ["f-type", "v-a-type", "ion-preference"]:
            if not os.path.isdir(os.path.join(root, sub)):
                return False
        return True

    if not have_lib(HMMrootdir):
        os.system("which MagicLamp.py > mainDir.txt")
        try:
            file = open("mainDir.txt")
            location = ""
            for i in file:
                location = i.rstrip()
            location = allButTheLast(location, "/")
            HMMrootdir = location + "/hmms/atp"
        except FileNotFoundError:
            pass

    if not have_lib(HMMrootdir):
        print("ATPGenie could not locate its HMM library. Please run the setup.sh script if\n"
              "you have Conda installed. Otherwise, please run setup-noconda.sh and put MagicLamp.py\n"
              "into your $PATH. The library must contain three subdirectories: f-type/, v-a-type/, ion-preference/.")
        os.system("rm -f HMMlib.txt mainDir.txt")
        raise SystemExit

    os.system("rm -f HMMlib.txt mainDir.txt")

    # Three sub-libraries
    HMMdir_f   = os.path.join(HMMrootdir, "f-type")
    HMMdir_v   = os.path.join(HMMrootdir, "v-a-type")
    HMMdir_ion = os.path.join(HMMrootdir, "ion-preference")

    args = parser.parse_known_args()[0]

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
        print("Looks like you did not provide an extension for your genomes/bins or assemblies, so ATPGenie does not "
              "know which files in the provided directory are FASTA files that you would like analyzed.")
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)
        # answer = input("Would you like ATPGenie to proceed and potentially overwrite files in this directory? (y/n): ")
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
                                print("Unfortunately, this is a problem for ATPGenie because it uses that character "
                                      "as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                except FileNotFoundError:
                    binFile = open("%s/%s" % (binDir, i), "r")
                    for line in binFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \\|")
                                print("Unfortunately, this is a problem for ATPGenie because it uses that character "
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
    # Subunit metadata (HMM accession -> complex_type, gene symbol, label)
    # =================================================================
    SUBUNIT_META = {
        # F-type
        "TIGR00962": ("F-type",   "atpA",  "F1_alpha"),
        "TIGR01039": ("F-type",   "atpD",  "F1_beta"),
        "TIGR01146": ("F-type",   "atpG",  "F1_gamma"),
        "TIGR01145": ("F-type",   "atpH",  "F1_delta"),
        "TIGR01216": ("F-type",   "atpC",  "F1_epsilon"),
        "TIGR01131": ("F-type",   "atpB",  "F0_a"),
        "TIGR01144": ("F-type",   "atpF",  "F0_b"),
        "TIGR01260": ("F-type",   "atpE",  "F0_c"),       # <- proteolipid (ion-preference target)
        # V/A-type
        "TIGR01042": ("V-A-type", "ntpA",  "V1_A"),
        "TIGR01043": ("V-A-type", "ntpA",  "A-type_A"),
        "TIGR01040": ("V-A-type", "ntpB",  "V1_B"),
        "TIGR01041": ("V-A-type", "ntpB",  "A-type_B"),
        "TIGR00309": ("V-A-type", "ntpD",  "V1_D"),
        "TIGR01100": ("V-A-type", "ntpC",  "V1_C"),
        "TIGR01101": ("V-A-type", "ntpF",  "V1_F"),
        "TIGR01147": ("V-A-type", "ntpG",  "V1_G"),
        "PF00137":   ("V-A-type", "ntpK",  "V0_K"),       # <- proteolipid (ion-preference target)
    }

    # The membrane proteolipid HMMs whose hits get re-scanned with ion HMMs
    PROTEOLIPID_HMMS = {"TIGR01260", "PF00137"}

    print("")
    # =================================================================
    # PASS 1 : profile each genome with f-type + v-a-type HMMs (--cut_tc)
    # =================================================================
    print("starting main pipeline...")
    print("Pass 1: profiling ATPase subunits with TIGRFAM + Pfam HMMs (hmmsearch --cut_tc)")

    HMMdirLS_f = [h for h in os.listdir(HMMdir_f) if h.endswith(".hmm")]
    HMMdirLS_v = [h for h in os.listdir(HMMdir_v) if h.endswith(".hmm")]
    all_pass1 = [(HMMdir_f, h, "F-type") for h in HMMdirLS_f] + \
                [(HMMdir_v, h, "V-A-type") for h in HMMdirLS_v]

    HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))

    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext:
            os.system("mkdir -p " + outDirectory + "/" + i + "-HMM")
            count = 0
            for hmmDir, hmm, complexType in all_pass1:
                count += 1
                perc = (count / len(all_pass1)) * 100
                sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                sys.stdout.flush()

                os.system(
                    "hmmsearch --cpu %d --cut_tc --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                    % (int(args.t), outDirectory, i, hmm, outDirectory, i, hmm, hmmDir, hmm,
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
                        # --cut_tc has already filtered; no extra threshold needed
                        if orf not in HMMdict[i]:
                            HMMdict[i][orf]["hmm"] = hmm
                            HMMdict[i][orf]["evalue"] = evalue
                            HMMdict[i][orf]["bit"] = bitscore
                            HMMdict[i][orf]["complex"] = complexType
                        else:
                            if bitscore > HMMdict[i][orf]["bit"]:
                                HMMdict[i][orf]["hmm"] = hmm
                                HMMdict[i][orf]["evalue"] = evalue
                                HMMdict[i][orf]["bit"] = bitscore
                                HMMdict[i][orf]["complex"] = complexType
            print("")

    # =================================================================
    # PASS 2 : ion-preference classification for proteolipid hits
    #   Best-of-four wins. No bitscore or E-value cutoff.
    # =================================================================
    print("Pass 2: ion-preference classification of proteolipid hits (atpE / K)")

    ionDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
    HMMdirLS_ion = [h for h in os.listdir(HMMdir_ion) if h.endswith(".hmm")]

    def ion_hmm_complex(acc):
        if acc.startswith("F-type"):
            return "F-type"
        if acc.startswith("V-type"):
            return "V-A-type"
        return "unknown"

    def confidence_label(delta):
        # Interpretation of within-complex bitscore difference (delta).
        # Each bit doubles the likelihood ratio: delta=10 -> ~1000x,
        # delta=20 -> ~1e6x. We label coarsely:
        try:
            d = float(delta)
        except (TypeError, ValueError):
            return "-"
        if d >= 20:
            return "high"
        if d >= 10:
            return "medium"
        if d >= 3:
            return "low"
        return "ambiguous"

    for i in binDirLS:
        if lastItem(i.split(".")) != args.bin_ext:
            continue

        # Collect proteolipid ORFs (best Pass-1 hit was atpE or K).
        candidates = []
        for orf in HMMdict[i]:
            best = HMMdict[i][orf]["hmm"]
            if best == "EMPTY":
                continue
            acc_p1 = best.split(".")[0]
            if acc_p1 in PROTEOLIPID_HMMS:
                candidates.append(orf)

        if not candidates:
            continue

        candFaa = "%s/ORF_calls/%s-proteolipid_candidates.faa" % (outDirectory, i)
        cf = open(candFaa, "w")
        for orf in candidates:
            seq = BinDict[i][orf]
            if seq != "EMPTY" and len(seq) > 0:
                cf.write(">" + orf + "\n" + str(seq) + "\n")
        cf.close()

        # Run all four ion HMMs against the candidates with no cutoff.
        # Collect per-ORF bitscores against every ion HMM, then derive:
        #   - winner = highest-bitscore HMM (its prefix sets ion_complex)
        #   - preference score = winner - same-complex-alternative
        # The within-complex delta is what tells you Na+ vs H+ confidence;
        # comparing across F vs V would conflate ion vs complex separation.
        perOrfScores = defaultdict(dict)  # orf -> {ion_hmm_acc: bitscore}
        perOrfEvalues = defaultdict(dict)
        for hmm in HMMdirLS_ion:
            acc = hmm.split(".")[0]
            tbl = "%s/%s-HMM/%s.tblout" % (outDirectory, i, hmm)
            os.system(
                "hmmsearch --cpu %d --tblout %s -o /dev/null %s/%s %s"
                % (int(args.t), tbl, HMMdir_ion, hmm, candFaa)
            )
            try:
                hmmout = open(tbl, "r")
            except FileNotFoundError:
                continue
            for line in hmmout:
                if not re.match(r'#', line):
                    ls = delim(line)
                    evalue = float(ls[4])
                    bitscore = float(ls[5])
                    orf = ls[0]
                    # Track best (per ion HMM) score for this ORF
                    if acc not in perOrfScores[orf] or bitscore > perOrfScores[orf][acc]:
                        perOrfScores[orf][acc] = bitscore
                        perOrfEvalues[orf][acc] = evalue

        # Resolve winner + within-complex preference score for each ORF
        for orf, scores in perOrfScores.items():
            if not scores:
                continue
            # Global best-of-four winner picks the complex
            winner_acc = max(scores, key=scores.get)
            winner_bit = scores[winner_acc]
            ion_complex = ion_hmm_complex(winner_acc)
            if "_Na_" in winner_acc:
                ion = "Na+"
                same_complex_alt = winner_acc.replace("_Na_pumping", "_H_pumping") \
                                             .replace("_Na_binding", "_H_pumping")
            elif "_H_" in winner_acc:
                ion = "H+"
                if winner_acc.startswith("F-type"):
                    same_complex_alt = "F-type_Na_pumping"
                else:
                    same_complex_alt = "V-type_Na_binding"
            else:
                ion = "ambiguous"
                same_complex_alt = None

            alt_bit = scores.get(same_complex_alt) if same_complex_alt else None
            if alt_bit is None:
                # No same-complex-alternative hit at all -> winner is the
                # only signal. Treat the alternative as scoring 0 bits,
                # so preference = full winner bitscore.
                pref = winner_bit
            else:
                pref = winner_bit - alt_bit

            ionDict[i][orf]["hmm"] = winner_acc + ".hmm"
            ionDict[i][orf]["bit"] = winner_bit
            ionDict[i][orf]["evalue"] = perOrfEvalues[orf][winner_acc]
            ionDict[i][orf]["ion"] = ion
            ionDict[i][orf]["ion_complex"] = ion_complex
            ionDict[i][orf]["alt_bit"] = alt_bit if alt_bit is not None else "-"
            ionDict[i][orf]["pref"] = round(pref, 2)
            ionDict[i][orf]["conf"] = confidence_label(pref)

        # Override Pass-1 complex_type for proteolipid ORFs whose ion-HMM
        # winner disagrees -- the ion HMM is the higher-resolution model
        # for the c/K rotor specifically. We DO NOT rewrite the gene/HMM
        # name or its bitscore/E-value, because those came from the actual
        # Pass-1 hit (PF00137 or TIGR01260) and rewriting them would create
        # a CSV row whose score doesn't match its labelled HMM. Instead we
        # only flip complex_type, and the heatmap categorization function
        # downstream uses (Pass-1 HMM, ion call, complex_type) jointly.
        for orf in candidates:
            ionRec = ionDict[i].get(orf)
            if not ionRec:
                continue
            ion_complex = ionRec.get("ion_complex", "unknown")
            if ion_complex in ("F-type", "V-A-type") and ion_complex != HMMdict[i][orf]["complex"]:
                HMMdict[i][orf]["complex"] = ion_complex

    # =================================================================
    # Initial summary CSV (all pass-1 hits + ion calls where applicable)
    # =================================================================
    out = open("%s/summary.csv" % outDirectory, "w")
    out.write("cell,ORF,HMM,evalue,bitscore,complex_type,ion_preference,ion_HMM,ion_winner_bit,ion_alt_bit,ion_preference_score,ion_confidence\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            ionRec = ionDict[key].get(j)
            if ionRec and ionRec.get("ion", "EMPTY") != "EMPTY":
                ion_pref = ionRec["ion"]
                ion_hmm  = ionRec["hmm"]
                ion_bit  = ionRec["bit"]
                ion_alt  = ionRec.get("alt_bit", "-")
                ion_pref_score = ionRec.get("pref", "-")
                ion_conf = ionRec.get("conf", "-")
            else:
                ion_pref = "-"; ion_hmm = "-"; ion_bit = "-"
                ion_alt = "-"; ion_pref_score = "-"; ion_conf = "-"
            out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                      str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "," +
                      HMMdict[key][j]["complex"] + "," + str(ion_pref) + "," + str(ion_hmm) + "," +
                      str(ion_bit) + "," + str(ion_alt) + "," + str(ion_pref_score) + "," +
                      str(ion_conf) + "\n")
    out.close()

    # =================================================================
    # Dereplicate (keep best-bitscore HMM per ORF across all libraries)
    # =================================================================
    summary = open(outDirectory + "/summary.csv", "r")
    SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] in ("cell", "ATPGenie"):
            continue
        if len(ls) < 12:
            continue
        cell = ls[0]
        orf = ls[1]
        hmm = ls[2]
        evalue = ls[3]
        hmmBit = ls[4]
        complexType = ls[5]
        ion_pref = ls[6]
        ion_hmm = ls[7]
        ion_bit = ls[8]
        ion_alt = ls[9]
        ion_pref_score = ls[10]
        ion_conf = ls[11]
        seq = BinDict[cell][orf]

        record = {
            "hmm": hmm, "hmmBit": hmmBit, "e": evalue, "seq": seq,
            "complex": complexType, "ion": ion_pref, "ion_hmm": ion_hmm, "ion_bit": ion_bit,
            "ion_alt": ion_alt, "ion_pref_score": ion_pref_score, "ion_conf": ion_conf,
        }

        if orf not in SummaryDict[cell]:
            SummaryDict[cell][orf] = record
        else:
            if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                SummaryDict[cell][orf] = record

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
                if len(RemoveDuplicates(k)) == 1:
                    orf = j + "_" + str(k[0])
                    rec = SummaryDict[i][orf]
                    if rec == "EMPTY" or rec.get("hmm", "EMPTY") == "EMPTY":
                        continue
                    out.write(",".join([i, orf, rec["hmm"], rec["e"], str(rec["hmmBit"]),
                                        rec["complex"], rec["ion"], rec["ion_hmm"], str(rec["ion_bit"]),
                                        str(rec["ion_alt"]), str(rec["ion_pref_score"]), str(rec["ion_conf"]),
                                        str(counter), str(rec["seq"])]) + "\n")
                    out.write(",".join(["#"] * 14) + "\n")
                    counter += 1
                else:
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)
                        rec = SummaryDict[i][orf]
                        if rec == "EMPTY" or rec.get("hmm", "EMPTY") == "EMPTY":
                            continue
                        out.write(",".join([i, orf, rec["hmm"], rec["e"], str(rec["hmmBit"]),
                                            rec["complex"], rec["ion"], rec["ion_hmm"], str(rec["ion_bit"]),
                                            str(rec["ion_alt"]), str(rec["ion_pref_score"]), str(rec["ion_conf"]),
                                            str(counter), str(rec["seq"])]) + "\n")
                    out.write(",".join(["#"] * 14) + "\n")
                    counter += 1
    out.close()

    # =================================================================
    # Filter out singleton proteolipid hits with no ion call
    #   (small membrane peptides are the most error-prone subunit;
    #    F1/V1 catalytic singletons are kept)
    # =================================================================
    clusterDict = defaultdict(lambda: defaultdict(list))
    summary = open("%s/summary-2.csv" % args.out, "r")
    for i in summary:
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            if len(ls) < 14:
                continue
            try:
                clu = int(ls[12])
            except ValueError:
                continue
            clusterDict[clu]["line"].append(ls)
            clusterDict[clu]["gene"].append(ls[2].split(".")[0])

    print("..")
    print("...")
    out = open("%s/summary-3.csv" % args.out, "w")
    out.write("file,ORF,gene,evalue,bit_score,complex_type,ion_preference,"
              "ion_HMM,ion_winner_bit,ion_alt_bit,ion_preference_score,ion_confidence,"
              "cluster_id,seq\n")
    for i in sorted(clusterDict.keys()):
        ls = clusterDict[i]["gene"]
        keep = True
        if len(RemoveDuplicates(ls)) == 1 and ls[0] in ("TIGR01260", "PF00137"):
            row = clusterDict[i]["line"][0]
            ion_pref = row[6]
            if ion_pref in ("-", "ambiguous", ""):
                keep = False
        if not keep:
            continue
        for j in clusterDict[i]["line"]:
            out.write(",".join(j) + "\n")
        out.write("####################################################\n")
    out.close()

    os.system("rm %s/summary.csv" % args.out)
    os.system("rm %s/summary-2.csv" % args.out)
    os.system("mv %s/summary-3.csv %s/atpgenie-summary.csv" % (args.out, args.out))

    os.system("mkdir -p %s/HMM_results" % outDirectory)
    os.system("rm -f %s/ORF_calls/*-prodigal.out" % outDirectory)
    os.system("rm -rf %s/HMM_results/*-HMM" % outDirectory)
    os.system("mv -f %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

    # =================================================================
    # Heatmap-compatible CSV
    #   Rows: each ATPase subunit category, with proteolipid split into
    #   Na+ / H+ / unclassified rows so the ion call shows up directly.
    # =================================================================
    print("....")
    print(".....")

    cats = [
        # F-type subunits
        "TIGR00962_atpA_F1alpha",
        "TIGR01039_atpD_F1beta",
        "TIGR01146_atpG_F1gamma",
        "TIGR01145_atpH_F1delta",
        "TIGR01216_atpC_F1epsilon",
        "TIGR01131_atpB_F0a",
        "TIGR01144_atpF_F0b",
        "TIGR01260_atpE_F0c_Na",
        "TIGR01260_atpE_F0c_H",
        "TIGR01260_atpE_F0c_unclassified",
        # V/A-type subunits
        "TIGR01042_ntpA_V1A",
        "TIGR01043_ntpA_AtypeA",
        "TIGR01040_ntpB_V1B",
        "TIGR01041_ntpB_AtypeB",
        "TIGR00309_ntpD_V1D",
        "TIGR01100_ntpC_V1C",
        "TIGR01101_ntpF_V1F",
        "TIGR01147_ntpG_V1G",
        "PF00137_ntpK_V0K_Na",
        "PF00137_ntpK_V0K_H",
        "PF00137_ntpK_V0K_unclassified",
    ]

    def category_for(hmm_acc, ion_pref, complex_type):
        # Proteolipid bucketing is driven by the corrected complex_type
        # (which the ion HMM may have flipped) rather than by which Pass-1
        # proteolipid HMM happened to win. This puts e.g. a PF00137 Pass-1
        # hit whose ion HMM identified it as F-type into the F0c row.
        if hmm_acc in ("TIGR01260", "PF00137"):
            if complex_type == "F-type":
                if ion_pref == "Na+": return "TIGR01260_atpE_F0c_Na"
                if ion_pref == "H+":  return "TIGR01260_atpE_F0c_H"
                return "TIGR01260_atpE_F0c_unclassified"
            if complex_type == "V-A-type":
                if ion_pref == "Na+": return "PF00137_ntpK_V0K_Na"
                if ion_pref == "H+":  return "PF00137_ntpK_V0K_H"
                return "PF00137_ntpK_V0K_unclassified"
        for c in cats:
            if c.startswith(hmm_acc + "_"):
                return c
        return None

    Dict = defaultdict(lambda: defaultdict(list))
    final = open("%s/atpgenie-summary.csv" % args.out, "r")
    for i in final:
        if re.match(r'#', i):
            continue
        ls = i.rstrip().split(",")
        if len(ls) < 14:
            continue
        if ls[0] in ("file", "bin", "assembly", "genome"):
            continue
        cell = ls[0]
        hmm_acc = ls[2].split(".")[0]
        complex_type = ls[5]
        ion_pref = ls[6]
        cat = category_for(hmm_acc, ion_pref, complex_type)
        if cat is None:
            continue
        Dict[cell][cat].append(hmm_acc)

    normDict = defaultdict(lambda: 'EMPTY')
    for i in os.listdir(args.bin_dir):
        if lastItem(i.split(".")) == args.bin_ext:
            try:
                file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                file = fasta(file)
                normDict[i] = len(file.keys())
            except FileNotFoundError:
                normDict[i] = 0

    outHeat = open("%s/atpgenie.heatmap.csv" % outDirectory, "w")
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
    print("Results are written to %s/atpgenie-summary.csv and %s/atpgenie.heatmap.csv" % (args.out, args.out))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()