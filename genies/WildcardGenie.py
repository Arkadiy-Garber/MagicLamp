#!/usr/bin/env python3
import time
from collections import defaultdict
import re
import os
import argparse
import textwrap
import sys
import subprocess
import statistics


def reverseComplement(seq):
    out = []
    for i in range(len(seq) - 1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


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


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: 'EMPTY')
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

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


parser = argparse.ArgumentParser(
    prog="YfGenie.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber;
    Arizona State University
    Please send comments and inquiries to agarber4@asu.edu

    *******************************************************
    '''))

parser.add_argument('-a', type=str, help="GenBank or RefSeq assembly accession. "
                                         "If unavailable, provide the input files via the -c, -g, and -p arguments",
                    default="")

parser.add_argument('-c', type=str, help="contig sequences in FASTA format", default="")

parser.add_argument('-g', type=str, help="GFF file", default="")

parser.add_argument('-p', type=str, help="protein seqences in FASTA format", default="")

parser.add_argument('-n', type=str, help="max distance between genes, in bp, to consider them part of the "
                                         "same gene neighborhood (default = 10000)", default=10000)

parser.add_argument('-i', type=str, help="locus identifier in the attributes field of the GFF file (default=Name)",
                    default="Name")

parser.add_argument('-d', type=str,
                    help="HMM database (path to HMM_dir). This argument is required if the --hmm flag is invoked",
                    default="")

parser.add_argument('-m', type=str, help="optional file containing preferred gene and pathway names for each HMM",
                    default="")

parser.add_argument('-y', type=str, help="single-column file of gene names of interest", default="")

parser.add_argument('-o', type=str, help="output basename. Default is basename of input file", default="")

parser.add_argument('-t', type=str, help="number parallel threads to use for hmmsearch (default = 1)", default=1)

parser.add_argument('--clean', type=str, help="remove all intermediate files and keep only the .csv summary",
                    const=True, nargs="?")

parser.add_argument('--hmm', type=str, help="run hmmsearch given a set of HMMs", const=True, nargs="?")

parser.add_argument('--gc', type=str, help="calculate GC content and AA frequencies", const=True, nargs="?")

parser.add_argument('--gff', type=str, help="rely on gene names as listed in the .gff file", const=True, nargs="?")

parser.add_argument('--add_cds', type=str, help="don't ask", const=True, nargs="?")

parser.add_argument('--dot_cds', type=str, help="don't ask", const=True, nargs="?")

parser.add_argument('--ncbi', type=str, help="data has already been downloaded from NCBI", const=True, nargs="?")

parser.add_argument('--cut_nc', type=str, help="use the HMM's noise cutoff to set all thresholding, instead of the default trusted cutoff", const=True, nargs="?")

parser.add_argument('--e', type=str, help="override baked-in bit score cutoffs with an evalue of 1E-6", const=True, nargs="?")

parser.add_argument('--x', type=str, help="set no significance cutoffs", const=True, nargs="?")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

if not args.gc and not args.hmm and not args.gff:
    print("Please select either or both of the flags for an analysis: --hmm, --gc. Exiting for now.")
    raise SystemExit

os.system("echo ~/miniconda3/etc/profile.d > path2conda")
path2condaFile = open("path2conda")
for i in path2condaFile:
    path2conda = i.rstrip()
os.system("rm path2conda")

CONTIGS = args.c
GFF = args.g
PROTS = args.p
if "" in [args.c, args.g, args.p]:
    if args.a == "":
        print("One of the required input files is missing. Please check command. Exiting for now.")
        raise SystemExit
    else:
        if not args.ncbi:
            os.system("echo %s > accession" % args.a)
            cmd = '. %s/conda.sh && conda activate bit && ' \
                  'bit-dl-ncbi-assemblies -w accession -f gff -j %s && ' \
                  'bit-dl-ncbi-assemblies -w accession -f protein -j %s && ' \
                  'bit-dl-ncbi-assemblies -w accession -f fasta -j %s && ' \
                  'conda deactivate' % (path2conda, args.t, args.t, args.t)
            subprocess.call(cmd, shell=True, executable='/bin/bash')
            os.system("gzip -d -f %s*gz" % args.a)

        CONTIGS = args.a + ".fa"
        GFF = args.a + ".gff"
        PROTS = args.a + ".faa"

try:
    file = open(CONTIGS)
    file = open(GFF)
    file = open(PROTS)
except FileNotFoundError:
    print("One or more of the required input files is missing. Exiting for now.")
    raise SystemExit

output = args.o
if output == "":
    output = allButTheLast(CONTIGS, ".")
    print(output)

if args.hmm:
    if "" in [args.d]:
        print("Please provide the HMM database. Exiting for now.")
        raise SystemExit
    print("Quering HMMs using hmmsearch")
    os.system("mkdir -p hmmhits_%s" % output)

    geneDict = defaultdict(lambda: '-')
    pathDict = defaultdict(lambda: '-')

    if args.m != "":
        file = open("%s" % args.m)
        for i in file:
            ls = (i.rstrip().split(","))
            if ls[0] != "gene":
                geneDict[ls[2].split(".")[0]] = ls[0]
                pathDict[ls[2].split(".")[0]] = ls[1]

    prots = open(PROTS)
    prots = fasta(prots)
    contigs = open(CONTIGS)
    contigs = fasta(contigs)
    gffDict = defaultdict(lambda: defaultdict(list))
    gff = open(GFF)
    for i in gff:
        if re.match(r'##FASTA', i):
            break
        else:
            if not re.match(r'#', i):
                ls = i.rstrip().split("\t")
                if ls[2] == "CDS" and not re.findall(r'pseudo=true', ls[8]):
                    try:
                        attr = args.i
                        attr1 = attr.split(",")[0]
                        attr2 = "ID"
                        # if attr != "Name":
                        #     attr2 = attr.split(",")[1]

                        if args.add_cds:
                            ID = ls[8].split(attr1 + "=cds-")[1].split(";")[0]
                            # if attr != "Name":
                            #     ID2 = ls[8].split(attr2 + "=cds-")[1].split(";")[0]
                        elif args.dot_cds:
                            ID = ls[8].split(attr1 + "=")[1].split(";")[0].split(".cds")[0]
                            print(ID)
                        else:
                            ID = ls[8].split(attr1 + "=")[1].split(";")[0]
                            # if attr != "Name":
                            #     ID2 = ls[8].split(attr2 + "=")[1].split(";")[0]
                    except IndexError:
                        print("Detected an issue with you GFF file formatting. See below line. Exiting for now.")
                        raise SystemExit
                    if ls[2] == "CDS":
                        gffDict[ID]["contig"] = [(ls[0])]
                        gffDict[ID]["coords"] = [int(ls[3]), int(ls[4])]
                        gffDict[ID]["strand"] = [ls[6]]
                        # if attr != "Name":
                        #     gffDict[ID2]["contig"] = [(ls[0])]
                        #     gffDict[ID2]["coords"] = [int(ls[3]), int(ls[4])]
                        #     gffDict[ID2]["strand"] = [ls[6]]

    CPU = args.t
    prots = open(PROTS)
    prots = fasta(prots)
    HMMdb = os.listdir(args.d)
    for i in HMMdb:
        if lastItem(i.split(".")) in ["HMM", "hmm"]:
            HMM = "%s/" % args.d + i
            TBLOUT = "hmmhits_%s/%s.tblout" % (output, i.split(".")[0])

            if args.x:
                os.system("hmmsearch --cpu %s --tblout %s %s %s > /dev/null 2>&1" % (CPU, TBLOUT, HMM, PROTS))
                # os.system("hmmsearch --cpu %s --tblout %s %s %s" % (CPU, TBLOUT, HMM, PROTS))
            else:
                if args.e:
                    os.system("hmmsearch -E 1E-6 --cpu %s --tblout %s %s %s > /dev/null 2>&1" % (CPU, TBLOUT, HMM, PROTS))
                    # os.system("hmmsearch -E 1E-6 --cpu %s --tblout %s %s %s" % (CPU, TBLOUT, HMM, PROTS))
                else:
                    if args.cut_nc:
                        os.system("hmmsearch --cut_nc --cpu %s --tblout %s %s %s > /dev/null 2>&1" % (CPU, TBLOUT, HMM, PROTS))
                        # os.system("hmmsearch --cut_nc --cpu %s --tblout %s %s %s" % (CPU, TBLOUT, HMM, PROTS))
                    else:
                        os.system("hmmsearch --cut_tc --cpu %s --tblout %s %s %s > /dev/null 2>&1" % (CPU, TBLOUT, HMM, PROTS))
                        # os.system("hmmsearch --cut_tc --cpu %s --tblout %s %s %s" % (CPU, TBLOUT, HMM, PROTS))

    iDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    cluDict = defaultdict(lambda: 'EMPTY')
    idxDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    coordDict = defaultdict(list)
    results = os.listdir("hmmhits_" + output)
    for i in results:
        if lastItem(i.split(".")) == "tblout":
            tblout = open("hmmhits_%s/%s" % (output, i))
            for j in tblout:
                if not re.match(r'#', j):
                    ls = (delim(j.rstrip()))
                    contig = gffDict[ls[0]]["contig"][0]
                    start = gffDict[ls[0]]["coords"][0]
                    end = gffDict[ls[0]]["coords"][0]
                    ID = ls[0]
                    desc = ls[2]
                    hmm = ls[3].split(".")[0]
                    bit = ls[5]

                    if ID not in iDict.keys():
                        iDict[ID]["desc"] = desc
                        iDict[ID]["hmm"] = hmm
                        iDict[ID]["start"] = start
                        iDict[ID]["end"] = end
                        iDict[ID]["contig"] = contig
                        iDict[ID]["bit"] = bit

                        idxDict[contig][start] = ID
                        idxDict[contig][end] = ID

                    else:
                        if float(bit) > float(iDict[ID]["bit"]):
                            iDict[ID]["desc"] = desc
                            iDict[ID]["hmm"] = hmm
                            iDict[ID]["start"] = start
                            iDict[ID]["end"] = end
                            iDict[ID]["contig"] = contig
                            iDict[ID]["bit"] = bit

                            idxDict[contig][start] = ID
                            idxDict[contig][end] = ID
                        else:
                            pass

                    if start not in coordDict[contig]:
                        coordDict[contig].append(start)
                    if end not in coordDict[contig]:
                        coordDict[contig].append(end)

    out = open("%s_hmmout.csv" % output, "w")
    out.write(
        "contig,start,end,strand,HMM_accession,HMM_description,locus,gene_call,pathway,bit,protein,gene_seq,gene_gc,A,G,V,I,L,M,F,Y,W,S,T,N,Q,C,P,R,H,K,D,E\n")
    for i in coordDict.keys():
        clus = (cluster(coordDict[i], int(args.n)))
        for j in clus:
            for k in j:
                ID = idxDict[i][k]
                hmm = iDict[ID]["hmm"]
                bit = iDict[ID]["bit"]
                contig = i
                start = gffDict[ID]["coords"][0]
                end = gffDict[ID]["coords"][1]
                seq = contigs[contig][start - 1:end]
                strand = gffDict[ID]["strand"][0]
                prot = prots[ID]

                if strand == "-":
                    seq = reverseComplement(seq)

                GC = (seq.count("G") + seq.count("C")) / len(seq)
                A = round(prot.count("A") / len(prot), 4)
                G = round(prot.count("G") / len(prot), 4)
                V = round(prot.count("V") / len(prot), 4)
                I = round(prot.count("I") / len(prot), 4)
                L = round(prot.count("L") / len(prot), 4)
                M = round(prot.count("M") / len(prot), 4)
                F = round(prot.count("F") / len(prot), 4)
                Y = round(prot.count("Y") / len(prot), 4)
                W = round(prot.count("W") / len(prot), 4)
                S = round(prot.count("S") / len(prot), 4)
                T = round(prot.count("T") / len(prot), 4)
                N = round(prot.count("N") / len(prot), 4)
                Q = round(prot.count("Q") / len(prot), 4)
                C = round(prot.count("C") / len(prot), 4)
                P = round(prot.count("P") / len(prot), 4)
                R = round(prot.count("R") / len(prot), 4)
                H = round(prot.count("H") / len(prot), 4)
                K = round(prot.count("K") / len(prot), 4)
                D = round(prot.count("D") / len(prot), 4)
                E = round(prot.count("E") / len(prot), 4)

                out.write(str(i) + "," + str(start) + "," + str(end) + "," + str(strand[0]) + "," +
                          str(hmm) + "," + str(iDict[ID]["desc"]) + "," +
                          ID + "," + str(geneDict[hmm]) + "," + str(pathDict[hmm]) + "," + str(bit) +
                          "," + str(prot) + "," + seq + "," + str(GC) + "," +
                          str(A) + "," + str(G) + "," + str(V) + "," + str(I) + "," +
                          str(L) + "," + str(M) + "," + str(F) + "," + str(Y) + "," +
                          str(W) + "," + str(S) + "," + str(T) + "," + str(N) + "," +
                          str(Q) + "," + str(C) + "," + str(P) + "," + str(R) + "," +
                          str(H) + "," + str(K) + "," + str(D) + "," + str(E) + "\n")
            out.write(
                "###############################################################################################\n")
    out.close()

################################################################################################################

if args.gc:
    print("Calculating GC content and amino acid frequences")
    contigs = open(CONTIGS)
    contigs = fasta(contigs)
    faa = open(PROTS)
    faa = fasta(faa)

    GClist = []
    totalLength = 0
    for i in contigs.keys():
        seq = contigs[i]
        GC = (seq.count("G") + seq.count("C")) / len(seq)
        GClist.append(GC)
        totalLength += len(seq)
    header = i

    try:
        GCmean = round(statistics.mean(GClist), 4)
        GCvar = round(statistics.stdev(GClist), 4)
    except statistics.StatisticsError:
        GCmean = round(GClist[0], 4)
        GCvar = 0

    Dict = defaultdict(list)
    SEQ = ""
    for i in faa.keys():
        seq = faa[i]
        SEQ += seq

    A = round(SEQ.count("A") / len(SEQ), 4)
    G = round(SEQ.count("G") / len(SEQ), 4)
    V = round(SEQ.count("V") / len(SEQ), 4)
    I = round(SEQ.count("I") / len(SEQ), 4)
    L = round(SEQ.count("L") / len(SEQ), 4)
    M = round(SEQ.count("M") / len(SEQ), 4)
    F = round(SEQ.count("F") / len(SEQ), 4)
    Y = round(SEQ.count("Y") / len(SEQ), 4)
    W = round(SEQ.count("W") / len(SEQ), 4)
    S = round(SEQ.count("S") / len(SEQ), 4)
    T = round(SEQ.count("T") / len(SEQ), 4)
    N = round(SEQ.count("N") / len(SEQ), 4)
    Q = round(SEQ.count("Q") / len(SEQ), 4)
    C = round(SEQ.count("C") / len(SEQ), 4)
    P = round(SEQ.count("P") / len(SEQ), 4)
    R = round(SEQ.count("R") / len(SEQ), 4)
    H = round(SEQ.count("H") / len(SEQ), 4)
    K = round(SEQ.count("K") / len(SEQ), 4)
    D = round(SEQ.count("D") / len(SEQ), 4)
    E = round(SEQ.count("E") / len(SEQ), 4)

    out = open("%s_gc.tsv" % output, "w")
    out.write("accession\tGC\tlength\tA\tG\tV\tI\tL\tM\tF\tY\tW\tS\tT\tN\tQ\tC\tP\tR\tH\tK\tD\tE\theader\n")
    out.write(str(output) + "\t" + str(GCmean) + "\t" + str(totalLength) + "\t" + str(A) + "\t" + str(G) + "\t" + str(
        V) + "\t" +
              str(I) + "\t" + str(L) + "\t" + str(M) + "\t" + str(F) + "\t" + str(Y) + "\t" + str(W) + "\t" + str(
        S) + "\t" +
              str(T) + "\t" + str(N) + "\t" + str(Q) + "\t" + str(C) + "\t" + str(P) + "\t" + str(R) + "\t" +
              str(H) + "\t" + str(K) + "\t" + str(D) + "\t" + str(E) + "\t" + header + "\n")
    out.close()
time.sleep(1)

if args.gff:
    print("Scanning GFF file for genes-of-interest")
    prots = open(PROTS)
    prots = fasta(prots)
    contigs = open(CONTIGS)
    contigs = fasta(contigs)
    genes = open(args.y)
    geneList = []
    for i in genes:
        geneList.append(i.rstrip())

    geneDict = defaultdict(list)
    gff = open(GFF)
    for i in gff:
        if re.match(r'##FASTA', i):
            break
        else:
            if not re.match(r'#', i):
                ls = i.rstrip().split("\t")
                if ls[2] == "CDS" and not re.findall(r'pseudo=true', ls[8]):
                    try:
                        gene = ls[8].split("gene=")[1].split(";")[0]
                    except IndexError:
                        try:
                            gene = ls[8].split("product=")[1].split(";")[0]
                        except IndexError:
                            gene = "na"
                    if gene in geneList:
                        geneDict[gene] = ls

    out = open("%s_genes.csv" % output, "w")
    out.write("contig,start,end,strand,locus,gene,prot,gene_seq,gene_gc,A,G,V,I,L,M,F,Y,W,S,T,N,Q,C,P,R,H,K,D,E\n")
    for i in geneDict.keys():
        ls = geneDict[i]
        if args.add_cds:
            ID = ls[8].split(args.i + "=cds-")[1].split(";")[0]
        else:
            ID = ls[8].split(args.i + "=")[1].split(";")[0]

        contig = ls[0]
        start = int(ls[3])
        end = int(ls[4])
        strand = ls[6]

        prot = prots[ID]
        seq = contigs[contig][start - 1:end]
        if strand == "-":
            seq = reverseComplement(seq)

        GC = (seq.count("G") + seq.count("C")) / len(seq)
        A = round(prot.count("A") / len(prot), 4)
        G = round(prot.count("G") / len(prot), 4)
        V = round(prot.count("V") / len(prot), 4)
        I = round(prot.count("I") / len(prot), 4)
        L = round(prot.count("L") / len(prot), 4)
        M = round(prot.count("M") / len(prot), 4)
        F = round(prot.count("F") / len(prot), 4)
        Y = round(prot.count("Y") / len(prot), 4)
        W = round(prot.count("W") / len(prot), 4)
        S = round(prot.count("S") / len(prot), 4)
        T = round(prot.count("T") / len(prot), 4)
        N = round(prot.count("N") / len(prot), 4)
        Q = round(prot.count("Q") / len(prot), 4)
        C = round(prot.count("C") / len(prot), 4)
        P = round(prot.count("P") / len(prot), 4)
        R = round(prot.count("R") / len(prot), 4)
        H = round(prot.count("H") / len(prot), 4)
        K = round(prot.count("K") / len(prot), 4)
        D = round(prot.count("D") / len(prot), 4)
        E = round(prot.count("E") / len(prot), 4)

        out.write(str(contig) + "," + str(start) + "," + str(end) + "," + str(strand) + "," +
                  ID + "," + str(i) + "," + str(prot) + "," + str(seq) + "," + str(GC) + "," +
                  str(A) + "," + str(G) + "," + str(V) + "," + str(I) + "," +
                  str(L) + "," + str(M) + "," + str(F) + "," + str(Y) + "," +
                  str(W) + "," + str(S) + "," + str(T) + "," + str(N) + "," +
                  str(Q) + "," + str(C) + "," + str(P) + "," + str(R) + "," +
                  str(H) + "," + str(K) + "," + str(D) + "," + str(E) + "\n")
time.sleep(1)

if args.clean:
    os.system("rm -rf hmmhits_%s" % output)
    if args.a:
        os.system("rm -rf accession")
        os.system("rm %s.fa" % args.a)
        os.system("rm %s.gff" % args.a)
        os.system("rm %s.faa" % args.a)

print("DONE!")
