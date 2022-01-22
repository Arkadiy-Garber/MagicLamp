#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def SUM(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count

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


def compare(queryList, targetList):
    counter = 0
    for i in targetList:
        if i in queryList:
            counter += 1
    if counter == len(targetList):
        return True
    else:
        return False


meta = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/hmm_info.csv", "r")
metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
operonDict = defaultdict(list)
operonMetaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
operonMainDict = defaultdict(list)
eDict = defaultdict(lambda: 'EMPTY')
for i in meta:
    ls = i.rstrip().split("\t")
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
        # metaDict[hmm]["importance"] = importance
        metaDict[hmm]["minHits"] = minHits
        metaDict[hmm]["bitcut"] = bitcut
        metaDict[hmm]["evalue"] = evalue
        metaDict[hmm]["minlength"] = minlength
        metaDict[hmm]["maxlength"] = maxlength
        operonDict[operon].append(hmm)
        operonMetaDict[operon] = minHits
        if importance == "1":
            operonMainDict[operon].append(hmm)

# for i in metaDict.keys():
#     print(i)
#     print(metaDict[i])
#     print(operonDict[metaDict[i]["operon"]])
#     print("")

BinDict = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/Genomes/GCA_006384225.1_ASM638422v1_genomic.fna-proteins.faa")
BinDict = fasta(BinDict)

hmmDir = os.listdir("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/HMM")
HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
for i in hmmDir:
    if re.findall(r'tblout', i):
        tblout = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/HMM/%s" % i)
        for line in tblout:
            if not re.match(r'#', line):
                ls = delim(line)
                evalue = float(ls[4])
                bit = float(ls[5])
                hmm = allButTheLast(i, ".")
                bitcut = float(metaDict[hmm]["bitcut"])
                orf = ls[0]
                seq = BinDict[orf]
                genome = "GCA_006384225.1_ASM638422v1_genomic.fna-proteins.faa"
                if bit > bitcut and len(seq) > metaDict[hmm]["minlength"] and len(seq) < metaDict[hmm]["maxlength"]:
                    # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS
                    if orf not in HMMdict[i]:
                        HMMdict[genome][orf]["hmm"] = hmm
                        HMMdict[genome][orf]["evalue"] = evalue
                        HMMdict[genome][orf]["bit"] = bit
                        HMMdict[genome][orf]["seq"] = seq
                        HMMdict[genome][orf]["bitcut"] = metaDict[hmm]["bitcut"]
                    else:
                        # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                        if bit > HMMdict[i][orf]["bit"]:
                            HMMdict[genome][orf]["hmm"] = hmm
                            HMMdict[genome][orf]["evalue"] = evalue
                            HMMdict[genome][orf]["bit"] = bit
                            HMMdict[genome][orf]["seq"] = seq
                            HMMdict[genome][orf]["bitcut"] = metaDict[hmm]["bitcut"]


out = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/summary.csv", "w")
out.write(
    "cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "seq" + "\n")
for key in HMMdict.keys():
    for j in HMMdict[key]:
        out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                  str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) +
                  "," + str(HMMdict[key][j]["bitcut"]) + "," + str(HMMdict[key][j]["seq"]) + "\n")

out.close()


# ****************************************** DEREPLICATION *********************************************************
summary = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/summary.csv")
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
CoordDict = defaultdict(lambda: defaultdict(list))
for i in SummaryDict.keys():
    if i != "cell":
        for j in SummaryDict[i]:
            contig = allButTheLast(j, "_")
            numOrf = lastItem(j.split("_"))
            CoordDict[i][contig].append(int(numOrf))

counter = 0
out = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/summary-2.csv", "w")
for i in CoordDict.keys():
    for j in CoordDict[i]:
        LS = (CoordDict[i][j])
        clusters = (cluster(LS, 2))
        for k in clusters:
            if len(RemoveDuplicates(k)) >= int(0):
                for l in RemoveDuplicates(k):
                    orf = j + "_" + str(l)
                    # print(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf][
                    #     "e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(
                    #     SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(
                    #     SummaryDict[i][orf]["seq"]) + "\n")]
                    out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                out.write("###############################################\n")
                counter += 1
out.close()

os.system("rm /Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/summary.csv")
os.system("mv /Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/summary-2.csv /Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/genie-summary.csv")

# ****************************** CUSTOM-RULE-BASED FILTERING ************************************************
clusterDict = defaultdict(lambda: defaultdict(list))
summary = open("/Users/arkadiygarber/Downloads/Files_for_Arkadiy_211013/genie-summary.csv")
for i in summary:
    ls = i.rstrip().split(",")
    if not re.match(r'#', i.rstrip()):
        # genome = ls[0]
        # orfCall = ls[1]
        hmmFile = ls[2]
        # evalue = ls[3]
        # bitscore = ls[4]
        # bitcut = ls[5]
        # cluster = ls[6]
        # seq = ls[7]
        gene = metaDict[hmmFile]["gene"]
        minHits = metaDict[hmmFile]["minHits"]
        minlength = metaDict[hmmFile]["minlength"]
        maxlength = metaDict[hmmFile]["maxlength"]
        operon = metaDict[hmmFile]["operon"]

        clusterDict[cluster]["ls"].append(ls)
        clusterDict[cluster]["hmms"].append(hmmFile)
        if hmmFile not in clusterDict[cluster]["unqHmms"]:
            clusterDict[cluster]["unqHmms"].append(hmmFile)

masterDict = defaultdict(lambda: defaultdict(list))
for i in clusterDict.keys():
    clusterNum = i
    ls = clusterDict[i]["ls"]
    hmms = (clusterDict[i]["hmms"])
    # print(hmms)
    unqHmms = (clusterDict[i]["unqHmms"])
    # print(unqHmms)
    operons = []
    for j in ls:
        operon = metaDict[j[2]]["operon"]
        if operon not in operons:
            operons.append(operon)
    # print(operons)
    passDict = defaultdict(lambda: 'EMPTY')
    for j in operons:
        operon = j
        minHits = operonMetaDict[operon]
        if len(unqHmms) >= int(minHits) and compare(unqHmms, operonMainDict[operon]):
            passDict[operon] = ls[0][6]
            # print(operon + "\t" + str(operonMetaDict[operon]) + "\t" + str(operonMainDict[operon]) + "\t" + str(
            #     compare(unqHmms, operonMainDict[operon])))
        for j in ls:
            operon = metaDict[j[2]]["operon"]
            if operon in passDict:
                if j[6] == passDict[operon]:
                    masterDict[j[0]][clusterNum].append(j)
                    # print(j)
        # print("")


for i in masterDict.keys():
    print(i)
    for j in masterDict[i]:
        for k in masterDict[i][j]:
            print(k)
        print("")






















