#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys

# TODO: ADD BITSCORE CUTOFF FOR EACH HMM HIT


def deExt(string):
    ls = string.split(".")
    ls2 = ls[0:len(ls)-1]
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
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
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
        prog="WspGenie.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************
    
        Developed by Arkadiy Garber;
        University of Montana, Biological Sciences
        Please send comments and inquiries to rkdgarber@gmail.com
        
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
    
        *******************************************************
        '''))


    parser.add_argument('-bin', type=str, help="FASTA format file")
    parser.add_argument('-format', type=str, help="is the input fasta file ORFs or contigs (orfs/contigs). "
                                                  "If contigs is chosen, then prodigal will be run."
                                                  "Default = contigs", default = "contigs")
    parser.add_argument('-out', type=str, help="name output directory (default=wspgenie_out)",
                        default="wspgenie_out")
    parser.add_argument('-hmm_dir', type=str, help='directory of HMMs. Provide the directory path (e.g. hmms/iron or hmms/lux) if you do not have conda installed and were not able to run the setup.sh script', default="NA")


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

    args = parser.parse_args()


    bits = open(args.hmm_dir + "/bitscores.txt")
    bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in bits:
        ls = i.rstrip().split("\t")
        bitDict[ls[0]]["gene"] = ls[1]
        bitDict[ls[0]]["bit"] = ls[2]

    if args.format == "contigs":
        os.system("prodigal -i %s -a %s-proteins.faa -o %s-prodigal.out -q" % (args.bin, args.bin, args.bin))
    else:
        os.system("mv %s %s-proteins.faa" % (args.bin, args.bin))

    os.system("mkdir " + args.out)

    hmms = os.listdir(args.hmm_dir)
    for i in hmms:
        if lastItem(i.split(".")) == "hmm":
            os.system(
                "hmmsearch --tblout %s/%s.tblout -o %s/%s.txt %s/%s %s-proteins.faa" % (args.out, i, args.out, i, args.hmm_dir, i, args.bin))


    results = os.listdir(args.out)
    resultsDict = defaultdict(lambda: defaultdict(list))
    for i in results:
        if lastItem(i.split(".")) == "tblout":
            result = open(args.out + "/" + i, "r")
            for line in result:
                if not re.match(r'#', line):
                    ls = delim(line)
                    evalue = float(ls[4])
                    bit = float(ls[5])
                    orf = ls[0]
                    query = ls[2]
                    threshold = bitDict[query]["bit"]
                    gene = bitDict[query]["gene"]
                    if float(bit) > float(threshold):
                        resultsDict[orf]["query"].append(query)
                        resultsDict[orf]["threshold"].append(threshold)
                        resultsDict[orf]["gene"].append(gene)
                        resultsDict[orf]["bit"].append(bit)

    out = open(args.out + "/wspgenie.csv", "w")
    for i in sorted(resultsDict.keys()):
        for j in range(0, len(resultsDict[i]["query"])):
            out.write(i + "," + str(resultsDict[i]["query"][j]) + "," + str(resultsDict[i]["gene"][j]) + "," + str(resultsDict[i]["bit"][j]) + "," + str(resultsDict[i]["threshold"][j]) + "\n")
    out.close()


    summaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    summary = open(args.out + "/wspgenie.csv", "r")
    for i in summary:
        ls = i.rstrip().split(",")
        orf = ls[0]
        contig = allButTheLast(orf, "_")
        orfNum = lastItem(orf.split("_"))
        domain = ls[1]
        gene = ls[2]
        bit = ls[3]
        threshold = ls[4]
        summaryDict[contig][orfNum]["domain"].append(domain)
        summaryDict[contig][orfNum]["gene"].append(gene)
        summaryDict[contig][orfNum]["bit"].append(bit)
        summaryDict[contig][orfNum]["threshold"].append(threshold)

    orfs = defaultdict(list)
    for i in summaryDict.keys():
        for j in summaryDict[i]:
            orfs[i].append(int(j))

    ORFS = defaultdict(list)
    for i in orfs.keys():
        clu = cluster(orfs[i], 2)
        for j in clu:
            if len(j) > 3:
                ORFS[i].append(j)

    out = open(args.out + "/wspgenie-2.csv", "w")
    out.write("orf" + "," + "gene" + "," + "domain" + "," + "domain_bitscore_ratio" + "," + "domain" + "," + "domain_bitscore_ratio" + "," + "domain" + "," + "domain_bitscore_ratio" + "\n")
    for i in ORFS:
        for j in ORFS[i]:
            for k in j:
                ORF = (i + "_" + str(k))
                out.write(ORF + "," + str(summaryDict[i][str(k)]["gene"][0]))
                for l in range(0, len(summaryDict[i][str(k)]["domain"])):
                    out.write("," + summaryDict[i][str(k)]["domain"][l] + "," + str(float(summaryDict[i][str(k)]["threshold"][l]) / float(summaryDict[i][str(k)]["bit"][l])))
                out.write("\n")
        out.write("#" + "\n")

    out.close()
    os.system("rm %s/wspgenie.csv" % args.out)
    os.system("mv %s/wspgenie-2.csv %s/wspgenie.csv" % (args.out, args.out))


    for i in results:
        if lastItem(i.split(".")) in ["txt", "tblout"]:
            os.system("rm %s/%s" % (args.out, i))


if __name__ == '__main__':
    main()


