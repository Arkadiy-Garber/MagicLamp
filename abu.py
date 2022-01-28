#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys
import time


def SUM(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count

def firstNum(string):
    outputNum = []
    for i in string:
        try:
            int(i)
            outputNum.append(i)
        except ValueError:
            break
    Num = "".join(outputNum)
    return Num


def Num(ls):
    outputNum = "0"
    for i in ls:
        try:
            int(i)
            outputNum = i
        except ValueError:
            pass
    return outputNum


def Strip(ls):
    outList = []
    for i in ls:
        gene = i.split("|")[0]
        outList.append(gene)
    return outList

def unique(ls, ls2):
    unqlist = []
    for i in ls:
        if i not in unqlist and i in ls2:
            unqlist.append(i)
    return len(unqlist)

def Unique(ls):
    unqList = []
    for i in ls:
        if i not in unqList:
            unqList.append(i)
    return unqList

def Unique2(ls):
    unqList = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in unqList:
            unqList.append(hmm)
    return unqList


def checkDFE1(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"]:
                count += 1
    return count

def checkDFE2(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"]:
                count += 1
    return count

def checkGACE(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"]:
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


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


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
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
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

parser = argparse.ArgumentParser(
    prog="abu.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-hmm', type=str, help="HMM file", default="NA")

parser.add_argument('-ref', type=str, help="reference protein sequences (e.g. nr or refseq database)", default="NA")

parser.add_argument('-out', type=str, help="name output directory", default="hmm_histograms")

parser.add_argument('-cpu', type=str, help="threads for hmmsearch (default = 4)", default="4")

parser.add_argument('-tru', type=str, help="comma-separated list of words and word-fragments in the annotations that coorrespond to true positives", default="NA")

parser.add_argument('-fal', type=str, help="comma-separated list of words and word-fragments in the annotations that coorrespond to false positives", default="NA")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

###################################################################
if "NA" in [args.hmm, args.ref]:
    print("Please provide all required files")
    raise SystemExit
else:
    pass

###################################################################
os.system("echo ${rscripts} > rscripts.txt")
file = open("rscripts.txt")
rscriptDir = ""
for i in file:
    rscriptDir = i.rstrip()

try:
    test = open(rscriptDir + "/Histogram.R")

except FileNotFoundError:
    os.system("which MagicLamp.py > mainDir.txt")

    file = open("mainDir.txt")
    location = ''
    for i in file:
        location = i.rstrip()
    location = allButTheLast(location, "/")

    rscriptDir = location + "/rscripts"

    try:
        test = open(rscriptDir + "/Histogram.R")
    except FileNotFoundError:
        print("FeGenie could not locate the required rscript directory. Please run the setup.sh script if "
              "you have Conda installed. Otherwise, please run the setupe-noconda.sh script and put MagicLamp.py into your $PATH")
        raise SystemExit

    os.system("rm mainDir.txt")

os.system("rm -f rscripts.txt")
###################################################################


os.system("hmmsearch -o %s/%s.tblout --noali --cpu %s %s %s" % (args.out, args.hmm, args.cpu, args.hmm, args.ref))

tblout = open("%s/%s.tblout" % (args.out, args.hmm))
out = open("%s/%s.mod.tblout" % (args.out, args.hmm), "w")
out.write("Fullseq_evalue|Fullseq_score|Fullseq_bias|Best1domain_evalue|Best1domain_score|Best1domain_bias|Numdom_exp|Numdom_N|Sequence|Desc\n")
for i in tblout:
    if not re.match(r'#', i):
        ls = delim(i.rstrip())
        if re.findall(r'>>', i):
            break
        else:
            if len(ls) > 8 and not re.match(r'-', ls[0]) and ls[0] != "E-value":
                lastField = (" ".join(ls[9:]))
                out.write(ls[0] + "|" + ls[1] + "|" + ls[2] + "|" + ls[3] + "|" + ls[4] + "|" + ls[5] + "|" + ls[6] + "|" + ls[7] + "|" + ls[8] + "|" + lastField + "\n")
out.close()

true = args.tru
trueList = true.split(",")
out = open("%s/true.csv" % (args.out), "w")
for i in trueList:
    out.write(i + "\n")
out.close()

false = args.fal
falseList = false.split(",")
out = open("%s/false.csv" % (args.out), "w")
for i in falseList:
    out.write(i + "\n")
out.close()

os.system("Rscript --vanilla %s/Histogram.R %s/ %s/%s.mod.tblout %s/true.csv %s/false.csv" % (rscriptDir, args.out, args.out, args.hmm, args.out, args.out))






