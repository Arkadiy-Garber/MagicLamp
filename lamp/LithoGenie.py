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


def checkFleet(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if memoryDict[hmm] in FLEET:
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


parser = argparse.ArgumentParser(
    prog="LithoGenie.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber and Gustavo Ramírez;
    University of Southern California, Earth Sciences
    Please send comments and inquiries to arkadiyg@usc.edu


...............**/((#((((((/**,..    ...                                   ..,(##,.......,,,,***,,,
...........,*/(%%&%####((*,..             .....,/%@@&%#*....              ,&&@@@@@&(....,.....,***,
.....,,,.,*(#&&&&%##(/,...       .........,#%%###%&@@@&&%(.....           ,/%&%##&(...,,,....,,,,
,,,,,,,*/#%&&&&%%(/,..        ........,,,*@@@% ,,,.  ,(%&%#/,......      .%/,.,,*(%%(,.....,,..,,,,
,,,,,*/#%&&&&&%#/,.           .....,,****%@@%..**,.    /&@@@&%*%@@/      (/     .*((/,.......,,..,,
,***/(%&&&&&%(,.             ..,,********&@&, .**...     ,(&@@@((##(      . .,...,***.............,
**/(#%&&@&%(*.   .*(/,    ....,,,,*&&%//(%       ...    ((%@&***/#(      .,*..       ............
//#%%&@&%(*.  *&&&&&&*.....,,,,(%%###%%&&,  .,**///,.    /%#/*,,**/(&%,,%/*....,,,,,,............
(#%%&&%(*,.  (@&%#(//(#/*,.*/(***(%%%####&%  ..,**///,.   %%(   .,,,,,(#%&%/,.,**/(((#%%%%,  ......
##%&%(*,.    @&(*,.  .*((((#%@@&&%%##(##/  ........,,,/(((%#/.    *   ***/(####((/*,,***///,.   .,/
#%%#/*..    .%/..    .,,. .,,.   ,/(/,      .****/((/,..,,.   ,.  .#. ..       ..    .*///(#%%(.   
##(%&&%(*    .     .  ..,,*///(/,.   .               .*/*. ,  ..         ..,*////*,..  ,(#%&@@@. 
(/&@@*         .   ..,,**///(/*.         ... .    .*/,.. *.     *@(    ..,,*/((((((*,... *(&@@@& 
*,#&&&%*.      .,    ...,,**///,,.     ,      ..     ,*..   .,,  ,. .*/(##(*.,*//((/,....  (&@@#*
..*####/*,..,*/*     ...,,,***,..      ,*  .,    ,*. .,..   .,. .        .*#&@@#.//*,.....  (&%(,
. .,....,*#%&%#.         .....       .  /* ..   .,/,       .,, #%      .,,***,  /&@%.,.  ....  #(**
         /##(/. ..,***//(%%*         .   #,      ,*,      .  ,%      ..,**///((*. #@@,  .....  ,(*,
       ,***,   .,*//(/*,.   .(      .,.   (#.              ,,       ..,**//(((((/, ,%@, ....    (##
      .*** ..,,,,,*,.    .    .    ..,.    ,#%.           ......     .,**/(((((/*.  .&& ...     (&%
      ,,. .....   ,**/***,.               . .. *     .,,**///////*,...,*/((((//**.   /@.        &@@
 **,,,,  .       .,**/((//////*,.  *%.# (( &%.&   .,,*///((((((((/******//////*,..   ,@*       *@@@
&&&&%*         . ,///////*****,,,.*%.&,%& %%.&( .,***////////////**,,,,,,****,,,..   .&.       @@@@
@@/.    .,.   ..,**********,,.. #*.& &(.%/.&..,,,************,,,..      ,,,,...    ,(       (@@&&
/,    .,***.   .////***,,,,,... ,.,# #/.(.,%  .,,,,,**,,,,,,... .,*/(/*. ....      *       #@@%
#*. .. .,***,.   ..**//***,,.....    , ,/ ,. /.   .......,,,,,,,,,**////*.          ..     .&@@%((%
%, .,,..,****.    .,,,,........     ,                .........   .......                 ,@@@%/,,**
&(     ...,,,.    .                                                                    .&&&%#**,,**
@,     .,,.                   /% &( (.,#                                            (%#(/***%&%**
@@&(*     .......        ..  /# /&.@(,@.(%    ....                                 .(&%(*,,**&@@&/*
*#%#/,.      .......          ,  * .* *../     .........                        *&@@@@@@@&@@&(/**
.,/##(***                      .                     .,,,,.     .  ,         ,**/**(##(#%@@&&/..,,,
..,/##*,*#%&&%(//*,......,*(((####&%%                  ,,,   . ,,            ...  .,,/(##(*,.......
...,*/*,....,*/(((/***/**,,,,,*(%@@&&%/    .......       ,  .. ..      .,.        ..,..............
....,,*,,.             ..../*,,*&@@@(.   .,**,,..                   ,*,         .........       .
...  .,,,,..             ..**/////#%#((.   .,**/**,,,.... ..          .*, //       ......          
.... .....                 ..,#@@@@@@&( .,  .,,*******,,,,..     .    ..,*(&&      ....    .       
,.... ...,.              .. ,&&&&&@@&%#,  .(%%(/*,,....,,*(##(,        ,,,*//      ...   ,/.       
*........,,..             .,(%&%#((((,.   ,#%#/*,...     .,/((*     ,@@&&&(.     ...    /&(      .,
#/,.......,,..     ,,.   ./((%&%#/,....   /##(*,. ..*(#*.    ,/*,. &@@@@@@@@&   ....   #&&(.    .. 
@&%(/, .....,.,..,(/*    /%&@@@&%/,,,/(*  /,.  /&@@&@@@@@#/,   .*,**,*#@@@@&&*  ...   /(/*/(   .   
@@@@(//////((#%#(*..   ,*#&@@@&(,/((.      .%@@@&&&@@@@@&&&(/. .,   #@@@/. ..    ,(/,..(#/.    
(#&@@@@@@@@@@@(,,,......,*/%@@@@#**##/*. //#@@/*//#&@@@@@&.   .(#@@&%##/**/*    ,((/,  */.     
,***/((######(/,,,,,,,,...,**(&@@@@&%(/*//%&&@@%#/*,,,,(&@@@@%//**%@@@%/*,,,*##.   *#((.   .,     .
,,...        ...,,,,//,....****,,*/%&@@@@@@@%(///*,,,,**/%%%(/%&&@&%%%#/* ./%@#.  /*.         ..,
.,,....       ....,,(#,,,,.......,*((#&&&(//////*****((##//(((##/****,  *%@@%..*#/,.        .,,,.
*,,,.....      ....../%,,,,,,,,,,..,*/(#(((///////(###((((/*,,.        #%@@@%/,,,.        .,....,
//*,,,.....      ....,/&@&/,,,,,,,*#/**,,/#%#(/**/(((//((%#/*((....... ,*#(//****,        ..,...,**
*///**,..  ...    ..,,**/#&@@@@&&%#//,,,,*%@@@//((/////(###%&/,,,.*@///*,.,,,,.       .,,,,,,***,
,*////*,,.  .....  ....,,*/((((((/**,*,,,,*(%@@@@@@@@#(#%&&&&&%##//%%%///,...*,      ..,,****,*,,,,

    Image design: Kazuki Takahashi (1996);
    ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/ 
    *******************************************************
    '''))

parser.add_argument('-bin_dir', type=str, help="directory of bins")
parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)")
parser.add_argument('-outdir', type=str, help="output directory (will be created if does not exist)",
                    default="genie_out")
parser.add_argument('-out', type=str, help="basename of output file (default = out)", default="out")
parser.add_argument('--contigs_source', type=str, help="are the provided contigs from a single organism (single)"
                                                       "or are you providing this program with metagenomic/metatranscriptomic assemblies (meta)? "
                                                       "(default=single)", default="single")
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

parser.add_argument('--d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                          "This number should be an integer and should reflect the maximum number of "
                                          "genes in between putative iron-related genes identified by the HMM database "
                                          "(default=3)", default=3)
parser.add_argument('--makeplots', type=str,
                    help="Would you like LithoGenie to make some figures from your data? y = yes, n = no (default = n). "
                         "If so, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                         "Warning: this part of the program is currently under beta-testing, and if there are any problems running Rscript, "
                         "or installing any of the required packages, you may get a bunch of error messages at the end. "
                         "This will not crash the program, however, and you can still expect to "
                         "see the main output (CSV files) from Genie.", default="n")

parser.add_argument('--element', type=str,
                    help="which element do you want to highlight in the heatmap output? [sulfur, hydrogen, "
                         "methane, nitrogen, oxygen, carbon-monoxide, C1compounds, carbon, urea,"
                         "halogenetated-compounds, arsenic, selenium, nitriles, iron, ALL] (default = ALL)",
                    default="ALL")

parser.add_argument('--only_heat', type=str, help="subvert all other functions of this scpript and only make the"
                                                  " heatmap. Please provide the output directory using the \'-outdir\' flag; this already must "
                                                  "contain the output summary file from a previously completed run. y = yes, n = no (default = n)",
                    default="n")

parser.add_argument('--cpu', type=int, help="number of threads to allow for hmmsearch (default = 1)", default=1)

parser.add_argument('--norm', type=str, help="normalize gene counts to total number of predicted ORFs in each genome or metagenome (y = yes, n = no) (default = yes)", default="yes")

# CHECKING FOR CONDA INSTALL
os.system("echo ${HMM_dir}/hmm-meta.txt > HMMlib.txt")
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
    parser.add_argument('-hmm_dir', type=str, help='directory of HMMs', default="NA")

    parser.add_argument('--R', type=str,
                        help="location of R scripts directory (note: this optional argument requires Rscript to be "
                             "installed on your system). The R scripts directory is in the same directory as the "
                             "Genie.py code", default="NA")

args = parser.parse_args()

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

if args.makeplots == 'y':
    if conda == 0:
        if args.R != "NA":
            print(".")
        else:
            print('Looks like you told LithoGenie to automatically generate R plots. '
                  'However, you have not provided the location of the directory that contains the R scripts '
                  '(as required of you because you did not go through the conda-based installation.')
            print("Exiting")
            raise SystemExit

if conda == 1:
    os.system("echo ${HMM_dir} > HMMlib.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = (i.rstrip())
        HMMdirLS = os.listdir(HMMdir)
    os.system("rm HMMlib.txt")
else:
    HMMdir = args.hmm_dir
    HMMdirLS = os.listdir(args.hmm_dir)

if args.only_heat != "y":
    try:
        os.listdir(args.outdir)
        print("Looks like you already have a directory with the name: " + args.outdir)
        print("To avoid overwriting potentially valuable files, FeGenie will now exit. "
              "Please delete or rename said directory and try running again.")
        print("Exiting")
        raise SystemExit
    except FileNotFoundError:
        print(".")
        os.system("mkdir %s" % args.outdir)
        outDirectory = "%s" % args.outdir
        outDirectoryLS = os.listdir("%s" % args.outdir)
else:
    outDirectory = "%s" % args.outdir
    outDirectoryLS = os.listdir("%s" % args.outdir)

# STARTING MAIN ALGORITHM
bits = open(HMMdir + "/hmm-meta.txt", "r")
bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
print("\nreading in HMM bitscore cut-offs...")
for i in bits:
    ls = i.rstrip().split("\t")
    bitDict[ls[0]]["bit"] = ls[1]
    bitDict[ls[0]]["process"] = ls[2]
    bitDict[ls[0]]["element"] = ls[3]
print("...")

if args.only_heat == "n":
    count = 0
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    out = open("%s/%s.csv" % (args.outdir, args.out), "w")
    out.write(
        "bin" + "," + "gene" + "," + "process" + "," + "substrate" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "sequence" + "\n")
    for i in binDirLS:
        HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        if not re.match(r'^\.', i) and i != (args.out + ".csv") and lastItem(i.split(".")) == args.bin_ext:

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
    summary = open("%s/%s.csv" % (args.outdir, args.out), "r")
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
    out = open("%s/%s-2.csv" % (args.outdir, args.out), "w")
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
    summary = open("%s/%s-2.csv" % (args.outdir, args.out), "r")
    for i in summary:
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            clu = int(ls[9])  ###############################
            clusterDict[clu]["line"].append(ls)
            clusterDict[clu]["gene"].append(ls[2])
            clusterDict[clu]["category"].append(ls[0])

    print("..")
    print("...")
    out = open("%s/%s-3.csv" % (args.outdir, args.out), "w")
    for i in sorted(clusterDict.keys()):  ###########################
        ls = (clusterDict[i]["gene"])

        if "EetA" in ls or "EetB" in ls or "Ndh2" in ls or "FmnB" in ls or "FmnA" in ls or "DmkA" in ls or "DmkB" in ls or "PplA" in ls:
            fleet = ["EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"]

            if unique(ls, fleet) < 5:  # If there are less than 5 FLEET genes in the cluster
                if len(remove2(ls, fleet)) < 1:  # If FLEET genes are the only ones in the cluster
                    pass
                else:  # If there are other genes in the cluster that are not FLEET
                    for j in clusterDict[i]["line"]:
                        if j[2] not in fleet:  # avoiding the fleet genes
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:  # if there are 5 or more of the FLEET genes present within cluster
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "MtoA" in ls or "MtrA" in ls or "MtrC_TIGR03507" in ls or "MtrB_TIGR03509" in ls:
            if "MtoA" in ls and "MtrB_TIGR03509" in ls:
                for j in clusterDict[i]["line"]:
                    if j[2] in ["MtrB_TIGR03509", "MtoA"]:
                        out.write(
                            j[0] + "," + "iron-oxidation" + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    else:
                        out.write(
                            j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                            j[7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            if "MtrA" in ls and "MtrB_TIGR03509" in ls:
                for j in clusterDict[i]["line"]:
                    if j[2] in ["MtrA", "MtrB_TIGR03509"]:
                        out.write(
                            j[0] + "," + "iron-reduction" + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    else:
                        out.write(
                            j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                            j[7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            elif "MtrB_TIGR03509.hmm" not in ls:
                pass
        ############################################################################################################
        elif "FoxA" in ls or "FoxB" in ls or "FoxC" in ls:
            foxabc = ["FoxA", "FoxB", "FoxC"]

            if unique(ls, foxabc) < 2:
                if len(remove2(ls, foxabc)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in foxabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:

                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "FoxE" in ls or "FoxY" in ls or "FoxZ" in ls:
            foxeyz = ["FoxE", "FoxY", "FoxZ"]

            if "FoxE" not in ls:

                if len(remove2(ls, foxeyz)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in foxeyz:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "Cyc1" in ls:
            if ("Cyc2_repCluster3" not in ls and "Cyc2_repCluster2" not in ls and "Cyc2_repCluster1" not in ls and
                        "MtoA" not in ls and "MtrB_TIGR03509" not in ls and "MtrA" not in ls and "MtrC_TIGR03507" not in ls):
                pass
        ############################################################################################################
        elif "CymA" in ls:
            if ("MtrB_TIGR03509" not in ls and "MtrA" not in ls and "MtoA" not in ls and "MtrC_TIGR03507" not in ls):
                pass
            else:
                if "MtrC_TIGR03507" in ls:
                    out.write(j[0] + "," + "iron-reduction" + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")
                else:
                    out.write(j[0] + "," + "iron-oxidation" + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

        ############################################################################################################
        elif "thiosulfate_oxidation_soxB" in ls or "thiosulfate_oxidation_soxY" in ls or "thiosulfate_oxidation_soxC" in ls \
                or "soxA" in ls or "soxX" in ls or "soxZ" in ls:
            soxabcxyz = ["thiosulfate_oxidation_soxB", "thiosulfate_oxidation_soxY", "thiosulfate_oxidation_soxC",
                         "soxA", "soxX", "soxZ"]

            if unique(ls, soxabcxyz) < 3:
                if len(remove2(ls, soxabcxyz)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in soxabcxyz:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "anaerobic_sulfite_reductase_asrA" in ls or "anaerobic_sulfite_reductase_asrB" in ls or "anaerobic_sulfite_reductase_asrC" in ls:
            asrabc = ["anaerobic_sulfite_reductase_asrA", "anaerobic_sulfite_reductase_asrB",
                      "anaerobic_sulfite_reductase_asrC"]

            if unique(ls, asrabc) < 2:
                if len(remove2(ls, asrabc)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in asrabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "particulate_methane_monooxygenase_pmoA" in ls or "partiuclate_methane_monooxygenase_pmoB" in ls or "partiuclate_methane_monooxygenase_pmoC" in ls:
            pmoabc = ["particulate_methane_monooxygenase_pmoA", "partiuclate_methane_monooxygenase_pmoB",
                      "partiuclate_methane_monooxygenase_pmoC"]

            if unique(ls, pmoabc) < 2:
                if len(remove2(ls, pmoabc)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in pmoabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "soluble_methane_monooxygenase_mmoB" in ls or "soluble_methane_monooxygenase_mmoD" in ls:
            mmoBD = ["soluble_methane_monooxygenase_mmoB", "soluble_methane_monooxygenase_mmoD"]

            if unique(ls, mmoBD) < 2:
                if len(remove2(ls, mmoBD)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in mmoBD:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "methyl_coenzymeM_reductase_mcrA" in ls or "methyl_coenzymeM_reductase_mcrB" in ls or "methyl_coenzymeM_reductase_mcrG" in ls:
            mcrabc = ["methyl_coenzymeM_reductase_mcrA", "methyl_coenzymeM_reductase_mcrB",
                      "methyl_coenzymeM_reductase_mcrG"]

            if unique(ls, mcrabc) < 2:
                if len(remove2(ls, mcrabc)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in mcrabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "Fe_nitrogenase_alpha" in ls or "Fe_nitrogenase_beta" in ls or "Fe_nitrogenase_delta" in ls:
            fenitrogenase = ["Fe_nitrogenase_alpha", "Fe_nitrogenase_beta", "Fe_nitrogenase_delta"]

            if unique(ls, fenitrogenase) < 2:
                if len(remove2(ls, fenitrogenase)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in fenitrogenase:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "MoFe_nitrogenase_nifD" in ls or "MoFe_nitrogenase_nifK" in ls or "nitrogenase_nifH" in ls:
            mofenitrogenase = ["MoFe_nitrogenase_nifD", "MoFe_nitrogenase_nifK", "nitrogenase_nifH"]

            if unique(ls, mofenitrogenase) < 2:
                if len(remove2(ls, mofenitrogenase)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in mofenitrogenase:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "V_containing_nitrogenase_Alpha" in ls or "V_containing_nitrogenase_Beta" in ls or "V_containing_nitrogenase_Delta" in ls:
            vnitrogenase = ["V_containing_nitrogenase_Alpha", "V_containing_nitrogenase_Beta",
                            "V_containing_nitrogenase_Delta"]

            if unique(ls, vnitrogenase) < 2:
                if len(remove2(ls, vnitrogenase)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in vnitrogenase:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrite_oxidoreductase_nxrA" in ls or "nitrite_oxidoreductase_nxrB" in ls:
            nxrab = ["nitrite_oxidoreductase_nxrA", "nitrite_oxidoreductase_nxrB"]

            if unique(ls, nxrab) < 2:
                if len(remove2(ls, nxrab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nxrab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrate_reductase_napA" in ls or "nitrate_reductase_napB" in ls:
            napab = ["nitrate_reductase_napA", "nitrate_reductase_napB"]

            if unique(ls, napab) < 2:
                if len(remove2(ls, napab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in napab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrate_reductase_narG" in ls or "nitrate_reductase_narH" in ls:
            nargh = ["nitrate_reductase_narG", "nitrate_reductase_narH"]

            if unique(ls, nargh) < 2:
                if len(remove2(ls, nargh)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nargh:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrite_reductase_nrfH" in ls or "nitrite_reductase_nrfA" in ls or "nitrite_reductase_nrfD" in ls:
            nrfhad = ["nitrite_reductase_nrfH", "nitrite_reductase_nrfA", "nitrite_reductase_nrfD"]

            if unique(ls, nrfhad) < 2:
                if len(remove2(ls, nrfhad)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nrfhad:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrite_reductase_nirB" in ls or "nitrite_reductase_nirD" in ls or "nitrite_reductase_nirK" in ls or "nitrite_reductase_nirS" in ls:
            nirbdks = ["nitrite_reductase_nirB", "nitrite_reductase_nirD", "nitrite_reductase_nirK",
                       "nitrite_reductase_nirS"]

            if unique(ls, nirbdks) < 2:
                if len(remove2(ls, nirbdks)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nirbdks:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitrite_reductase_nirB" in ls or "nitrite_reductase_nirD" in ls or "nitrite_reductase_nirK" in ls or "nitrite_reductase_nirS" in ls:
            nirbdks = ["nitrite_reductase_nirB", "nitrite_reductase_nirD", "nitrite_reductase_nirK",
                       "nitrite_reductase_nirS"]

            if unique(ls, nirbdks) < 2:
                if len(remove2(ls, nirbdks)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nirbdks:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "nitric_oxide_reductase_norB" in ls or "nitric_oxide_reductase_norC" in ls:
            norbc = ["nitric_oxide_reductase_norB", "nitric_oxide_reductase_norC"]

            if unique(ls, norbc) < 2:
                if len(remove2(ls, norbc)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in norbc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################

        elif "nitrous_oxide_reductase_nosD" in ls or "nitrous_oxide_reductase_nosZ" in ls:
            nosdz = ["nitrous_oxide_reductase_nosD", "nitrous_oxide_reductase_nosZ"]

            if unique(ls, nosdz) < 2:
                if len(remove2(ls, nosdz)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nosdz:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "CytCoxidase_aa3_coxA" in ls or "CytCoxidase_aa3_coxB" in ls:
            coxab = ["CytCoxidase_aa3_coxA", "CytCoxidase_aa3_coxB"]

            if unique(ls, coxab) < 2:
                if len(remove2(ls, coxab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in coxab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "CytCoxidase_cbb3_ccoN" in ls or "CytCoxidase_cbb3_ccoO" in ls or "CytCoxidase_cbb3_ccoP" in ls:
            cconop = ["CytCoxidase_cbb3_ccoN", "CytCoxidase_cbb3_ccoO", "CytCoxidase_cbb3_ccoP"]

            if unique(ls, cconop) < 2:
                if len(remove2(ls, cconop)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in cconop:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "Cyt_bo3_ubiquinol_oxidase_CyoA" in ls or "Cyt_bo3_ubiquinol_oxidase_CyoD" in ls or "Cyt_bo3_ubiquinol_oxidase_CyoE" in ls:
            cyoade = ["Cyt_bo3_ubiquinol_oxidase_CyoA", "Cyt_bo3_ubiquinol_oxidase_CyoD",
                      "Cyt_bo3_ubiquinol_oxidase_CyoE"]

            if unique(ls, cyoade) < 2:
                if len(remove2(ls, cyoade)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in cyoade:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "Cyt_bd_oxidase_CydA" in ls or "Cyt_bd_oxidase_CydB" in ls:
            cydab = ["Cyt_bd_oxidase_CydA", "Cyt_bd_oxidase_CydB"]

            if unique(ls, cydab) < 2:
                if len(remove2(ls, cydab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in cydab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "Cyt_aa3_quinol_oxidase_QoxA" in ls or "Cyt_aa3_quinol_oxidase_QoxB" in ls:
            qoxab = ["Cyt_aa3_quinol_oxidase_QoxA", "Cyt_aa3_quinol_oxidase_QoxB"]

            if unique(ls, qoxab) < 2:
                if len(remove2(ls, qoxab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in qoxab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        ############################################################################################################
        elif "methylamine_dehydrogenase_mauA" in ls or "methylamine_dehydrogenase_mauB" in ls:
            mauab = ["methylamine_dehydrogenase_mauA", "methylamine_dehydrogenase_mauB"]

            if unique(ls, mauab) < 2:
                if len(remove2(ls, mauab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in mauab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "formate_dehydrogenase_alpha" in ls or "formate_dehydrogenase_beta" in ls or "formate_dehydrogenase_gamma" in ls:
            formatedehyd = ["formate_dehydrogenase_alpha", "formate_dehydrogenase_beta", "formate_dehydrogenase_gamma"]

            if unique(ls, formatedehyd) < 2:
                if len(remove2(ls, formatedehyd)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in formatedehyd:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "carbon_monoxide_dehydrogenase_coxS" in ls or "carbon_monoxide_dehydrogenase_coxM" in ls or "carbon_monoxide_dehydrogenase_coxL" in ls:
            coxsml = ["carbon_monoxide_dehydrogenase_coxS", "carbon_monoxide_dehydrogenase_coxM",
                      "carbon_monoxide_dehydrogenase_coxL"]

            if unique(ls, coxsml) < 2:
                if len(remove2(ls, coxsml)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in coxsml:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")


        ############################################################################################################
        elif "wood_ljungdahl_codhD" in ls or "wood_ljungdahl_codhC" in ls or "wood_ljungdahl_codh_catalytic" in ls:
            codhdccat = ["wood_ljungdahl_codhD", "wood_ljungdahl_codhC", "wood_ljungdahl_codh_catalytic"]

            if unique(ls, codhdccat) < 2:
                if len(remove2(ls, codhdccat)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in codhdccat:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")


        ############################################################################################################
        elif "acetate_citrate_lyase_rTCA_aclA" in ls or "acetate_citrate_lyase_rTCA_aclB" in ls:
            aclab = ["acetate_citrate_lyase_rTCA_aclA", "acetate_citrate_lyase_rTCA_aclB"]

            if unique(ls, aclab) < 2:
                if len(remove2(ls, aclab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in aclab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "urease_ureA" in ls or "urease_ureB" in ls or "urease_ureC" in ls:
            ureabc = ["urease_ureA", "urease_ureB", "urease_ureC"]

            if unique(ls, ureabc) < 2:
                if len(remove2(ls, ureabc)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in ureabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "arsenite_oxidase_small" in ls or "arsenite_oxidase_large" in ls:
            arsenite = ["arsenite_oxidase_small", "arsenite_oxidase_large"]

            if unique(ls, arsenite) < 2:
                if len(remove2(ls, arsenite)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in arsenite:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "DMSO_reductase_typeII_pcrA" in ls or "DMSO_reductase_typeII_pcrB" in ls:
            pcrab = ["DMSO_reductase_typeII_pcrA", "DMSO_reductase_typeII_pcrB"]

            if unique(ls, pcrab) < 2:
                if len(remove2(ls, pcrab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in pcrab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "selenate_reductase_ygfM" in ls or "selenate_reductase_MoBindingSubunit" in ls or "selenate_reductase_ygfK" in ls:
            ygfmk = ["urease_ureA", "urease_ureB", "urease_ureC"]

            if unique(ls, ygfmk) < 2:
                if len(remove2(ls, ygfmk)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in ygfmk:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "nitrile_hydratase_nthA" in ls or "nitrile_hydratase_nthB" in ls:
            nthab = ["nitrile_hydratase_nthA", "nitrile_hydratase_nthB"]

            if unique(ls, nthab) < 2:
                if len(remove2(ls, nthab)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in nthab:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        ############################################################################################################
        elif "adenylyl_sulfate_kinase_cysC_apsK" in ls or "Sulfate_andenyltransferase_CysN" in ls:
            cyscn = ["adenylyl_sulfate_kinase_cysC_apsK", "Sulfate_andenyltransferase_CysN"]

            if unique(ls, cyscn) < 2:
                if len(remove2(ls, cyscn)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[2] not in cyscn:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                    6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")


        ############################################################################################################
        elif "dissimilatory_sulfite_reductase_DsrD" in ls or "dissimilatory_sulfite_reductase-sulfur_oxidation_dsrA" in ls \
                or "dissimilatory_sulfite_reductase-sulfur_oxidation_dsrB" in ls or "dsrK" in ls or "dsrM" in ls \
                or "dsrE" in ls or "dsrF" in ls or "dsrH" in ls or "dsrC" in ls:
            dsr = ["dissimilatory_sulfite_reductase_DsrD", "dissimilatory_sulfite_reductase-sulfur_oxidation_dsrA",
                   "dissimilatory_sulfite_reductase-sulfur_oxidation_dsrB", "dsrK", "dsrM", "dsrE", "dsrF", "dsrH",
                   "dsrC"]
            if unique(ls, ["dsrE", "dsrF", "dsrH"]):
                for j in clusterDict[i]["line"]:
                    if j[2] not in dsr:
                        out.write(
                            j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                            j[
                                7] + "," + j[8] + "," + j[9] + "\n")
                    else:
                        out.write(
                            "sulfur-oxidation" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," +
                            j[6] + "," + j[
                                7] + "," + j[8] + "," + j[9] + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    if j[2] not in dsr:
                        out.write(
                            j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                            j[
                                7] + "," + j[8] + "," + j[9] + "\n")
                    else:
                        out.write(
                            "sulfite-reduction" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," +
                            j[6] + "," + j[
                                7] + "," + j[8] + "," + j[9] + "\n")

        ############################################################################################################
        else:
            linels = (clusterDict[i]["line"])
            for j in linels:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                    7] + "," + j[8] + "," + j[9] + "\n")

            out.write(
                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    out.close()

    os.system("rm %s/%s.csv" % (args.outdir, args.out))
    os.system("rm %s/%s-2.csv" % (args.outdir, args.out))
    os.system("mv %s/%s-3.csv %s/%s-summary.csv" % (args.outdir, args.out, args.outdir, args.out))

    # ****************************** FILTERING OUT FALSE POSITIVES *************************************
    out = open("%s/%s-summary-fixed.csv" % (args.outdir, args.out), "w")
    clusterDict = defaultdict(list)
    memoryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    summary = open("%s/%s-summary.csv" % (args.outdir, args.out))
    for i in summary:
        if not re.match(r'#', i):
            try:
                ls = i.rstrip().split(",")
                clu = ls[9]
                cat = ls[1]
                dataset = ls[3]
                orf = ls[4]
                gene = ls[2]
                element = ls[0]
                e = ls[5]
                bit = ls[6]
                cutoff = ls[7]
                seq = ls[8]
                clusterDict[clu].append(gene + "|" + dataset + "|" + orf)
                memoryDict[dataset][orf]["cat"] = cat
                memoryDict[dataset][orf]["gene"] = gene
                memoryDict[dataset][orf]["bit"] = bit
                memoryDict[dataset][orf]["cutoff"] = cutoff
                memoryDict[dataset][orf]["clu"] = clu
                memoryDict[dataset][orf]["seq"] = seq
                memoryDict[dataset][orf]["e"] = e
                memoryDict[dataset][orf]["element"] = element
            except IndexError:
                pass

    FLEET = ["EetA", "EetB", "FmnA", "DmkA", "FmnB", "PplA", "Ndh2", "DmkB"]
    for i in clusterDict.keys():
        out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        for j in clusterDict[i]:
            gene = j.split("|")[0]
            dataset = j.split("|")[1]
            orf = j.split("|")[2]
            cat = memoryDict[dataset][orf]["cat"]
            if memoryDict[dataset][orf]["gene"] in FLEET:
                if checkFleet(clusterDict[i]) < 6:
                    pass
                else:
                    out.write(memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["cat"] + "," +
                              memoryDict[dataset][orf]["gene"] + "," + dataset + "," +
                              orf + "," + memoryDict[dataset][orf]["e"] + "," + memoryDict[dataset][orf]["bit"] + "," +
                              memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["seq"] + "," +
                              memoryDict[dataset][orf]["clu"] + "\n")

            else:
                out.write(memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["cat"] + "," +
                          memoryDict[dataset][orf]["gene"] + "," + dataset + "," +
                          orf + "," + memoryDict[dataset][orf]["e"] + "," + memoryDict[dataset][orf]["bit"] + "," +
                          memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["seq"] + "," +
                          memoryDict[dataset][orf][
                              "clu"] + "\n")
    out.close()

    os.system("mv %s/%s-summary-fixed.csv %s/%s-summary.csv" % (args.outdir, args.out, args.outdir, args.out))

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
                depth = open("%s/%s.depth" % (args.outdir, cell))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]
                        total += LS[2]
                normDict[cell] = total/1000000

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (args.outdir, cell, string))
                print("processing... " + cell)
                depth = open("%s/%s.depth" % (args.outdir, cell))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]
                        total += float(LS[2])
                normDict[cell] = total / 1000000

        os.system("mkdir %s/contigDepths" % args.outdir)
        os.system("mv %s/*depth %s/contigDepths/" % (args.outdir, args.outdir))

        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            elements = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        contig = allButTheLast(orf, "_")
                        gene = ls[2]
                        Dict[cell][process].append(float(depthDict[cell][contig]))

            outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write(str(SUM(Dict[j][i])/normDict[j]) + ",")
                outHeat.write("\n")

            outHeat.close()

        else:
            elementDict = defaultdict(lambda: defaultdict(list))
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            for i in bits:
                ls = (i.rstrip().split("\t"))
                if ls[0] != "hmm":
                    if ls[3] == args.element:
                        elementDict[ls[3]][ls[2]].append(ls[0])

            catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in elementDict.keys():
                for j in elementDict[i]:
                    for k in elementDict[i][j]:
                        catDict[k] = j

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        element = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        contig = allButTheLast(orf, "_")
                        gene = ls[2]
                        reaction = ls[1]
                        Dict[cell][gene].append(float(depthDict[cell][contig]))

            cats = []
            CatOrgDict = defaultdict(list)
            for i in Dict.keys():
                for j in Dict[i]:
                    if len(catDict[j]) > 0:
                        if j not in cats:
                            cats.append(j)
                            CatOrgDict[catDict[j]].append(j)

            outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in CatOrgDict.keys():
                for j in CatOrgDict[i]:
                    outHeat.write(i + "--" + j + ",")
                    for k in sorted(Dict.keys()):
                        outHeat.write(str(SUM(Dict[k][j])) + ",")
                    outHeat.write("\n")

            outHeat.close()

    # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE
    elif args.bam != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

        try:
            total = 0
            depth = open("%s.depth" % (args.bam))
            for k in depth:
                LS = k.rstrip().split("\t")
                if LS[0] != "contigName":
                    depthDict[LS[0]] = LS[2]
                    total += float(LS[2])

        except FileNotFoundError:
            os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.bam, args.bam))
            depth = open("%s.depth" % (args.bam))
            total = 0
            for k in depth:
                LS = k.rstrip().split("\t")
                if LS[0] != "contigName":
                    depthDict[LS[0]] = LS[2]
                    total += float(LS[2])

        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            elements = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        contig = allButTheLast(orf, "_")
                        gene = ls[2]
                        Dict[cell][process].append(float(depthDict[contig]))

            outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
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

        else:

            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)

            if args.element in elements:
                elementDict = defaultdict(lambda: defaultdict(list))
                bits = open(HMMdir + "/hmm-meta.txt", "r")
                for i in bits:
                    ls = (i.rstrip().split("\t"))
                    if ls[0] != "hmm":
                        if ls[3] == args.element:
                            elementDict[ls[3]][ls[2]].append(ls[0])

                catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                for i in elementDict.keys():
                    for j in elementDict[i]:
                        for k in elementDict[i][j]:
                            catDict[k] = j

                Dict = defaultdict(lambda: defaultdict(list))
                final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
                for i in final:
                    ls = (i.rstrip().split(","))
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                        if not re.match(r'#', i):
                            element = ls[0]
                            cell = ls[3]
                            orf = ls[4]
                            contig = allButTheLast(orf, "_")
                            gene = ls[2]
                            reaction = ls[1]
                            Dict[cell][gene].append(float(depthDict[contig]))

                cats = []
                CatOrgDict = defaultdict(list)
                for i in Dict.keys():
                    for j in Dict[i]:
                        if len(catDict[j]) > 0:
                            if j not in cats:
                                cats.append(j)
                                CatOrgDict[catDict[j]].append(j)

                outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
                outHeat.write("X" + ',')
                for i in sorted(Dict.keys()):
                    outHeat.write(i + ",")
                outHeat.write("\n")

                for i in CatOrgDict.keys():
                    for j in CatOrgDict[i]:
                        outHeat.write(i + "--" + j + ",")
                        for k in sorted(Dict.keys()):
                            outHeat.write(str(SUM(Dict[k][j])) + ",")
                        outHeat.write("\n")

                outHeat.close()
                print('......')
                print(".......")
                print("Finished!")

            else:
                print("Looks like the element you have chosen is not one that is recognized by LithoGenie. Please"
                      "try again by choosing another element, or checking your spelling.")

    # GENE COUNTS-BASED ABUNDANCE
    else:
        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            elements = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        gene = ls[2]
                        Dict[cell][process].append(gene)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in os.listdir(args.bin_dir):
                if lastItem(i.split(".")) == args.bin_ext:
                    file = open("%s/%s-proteins.faa" % (args.bin_dir, i), "r")
                    file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/%s.%s.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
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
            elementDict = defaultdict(lambda: defaultdict(list))
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            for i in bits:
                ls = (i.rstrip().split("\t"))
                if ls[0] != "hmm":
                    if ls[3] == args.element:
                        elementDict[ls[3]][ls[2]].append(ls[0])

            catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in elementDict.keys():
                for j in elementDict[i]:
                    for k in elementDict[i][j]:
                        catDict[k] = j

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        element = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        gene = ls[2]
                        reaction = ls[1]
                        Dict[cell][gene].append(orf)

            cats = []
            CatOrgDict = defaultdict(list)
            for i in Dict.keys():
                for j in Dict[i]:
                    if len(catDict[j]) > 0:
                        if j not in cats:
                            cats.append(j)
                            CatOrgDict[catDict[j]].append(j)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in binDirLS:
                if lastItem(i.split(".")) == args.bin_ext:
                    file = open("%s/%s-proteins.faa" % (binDir, i), "r")
                    file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/%s.%s.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in CatOrgDict.keys():
                for j in CatOrgDict[i]:
                    outHeat.write(i + "--" + j + ",")
                    for k in sorted(Dict.keys()):
                        if not re.match(r'#', j):
                            if args.norm == "y":
                                outHeat.write(str((len(Dict[k][j]) / int(normDict[k])) * float(100)) + ",")
                            else:
                                outHeat.write(str((len(Dict[k][j]))) + ",")
                    outHeat.write("\n")

            outHeat.close()
            print('......')
            print(".......")
            print("Finished!")


else:
    print("....")
    print(".....")
    # COVERAGE-BASED ABUNDANCE
    if args.bams != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
        normDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
        BAMmap = open(args.bams)
        for i in BAMmap:
            string = ''
            ls = i.rstrip().split("\t")
            cell = ls[0]
            for j in ls[1:]:
                string += " "
                string += j

            try:
                print("processing... " + cell)
                depth = open("%s/contigDepths/%s.depth" % (args.outdir, cell))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]
                        total += float(LS[2])
                normDict[cell] = total

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth%s" % (cell, string))
                print("processing... " + cell)
                depth = open("%s/contigDepths/%s.depth" % (args.outdir, cell))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]
                        total += float(LS[2])
                normDict[cell] = total

        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        contig = allButTheLast(orf, "_")
                        gene = ls[2]
                        Dict[cell][process].append(float(depthDict[cell][contig]))

            outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write(str(SUM(Dict[j][i])/normDict[j]) + ",")
                outHeat.write("\n")

            outHeat.close()
            print('......')
            print(".......")
            print("Finished!")

        else:
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)

            if args.element in elements:
                elementDict = defaultdict(lambda: defaultdict(list))
                bits = open(HMMdir + "/hmm-meta.txt", "r")
                for i in bits:
                    ls = (i.rstrip().split("\t"))
                    if ls[0] != "hmm":
                        if ls[3] == args.element:
                            elementDict[ls[3]][ls[2]].append(ls[0])

                catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                for i in elementDict.keys():
                    for j in elementDict[i]:
                        for k in elementDict[i][j]:
                            catDict[k] = j

                Dict = defaultdict(lambda: defaultdict(list))
                final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
                for i in final:
                    ls = (i.rstrip().split(","))
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                        if not re.match(r'#', i):
                            element = ls[0]
                            cell = ls[3]
                            orf = ls[4]
                            contig = allButTheLast(orf, "_")
                            gene = ls[2]
                            reaction = ls[1]
                            Dict[cell][gene].append(float(depthDict[cell][contig]))

                cats = []
                CatOrgDict = defaultdict(list)
                for i in Dict.keys():
                    for j in Dict[i]:
                        if len(catDict[j]) > 0:
                            if j not in cats:
                                cats.append(j)
                                CatOrgDict[catDict[j]].append(j)

                outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
                outHeat.write("X" + ',')
                for i in sorted(Dict.keys()):
                    outHeat.write(i + ",")
                outHeat.write("\n")

                for i in CatOrgDict.keys():
                    for j in CatOrgDict[i]:
                        outHeat.write(i + "--" + j + ",")
                        for k in sorted(Dict.keys()):
                            outHeat.write(str(SUM(Dict[k][j])) + ",")
                        outHeat.write("\n")

                outHeat.close()
                print('......')
                print(".......")
                print("Finished!")

            else:
                print("Looks like the element you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling.")

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

        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(cats)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        contig = allButTheLast(orf, "_")
                        gene = ls[2]
                        Dict[cell][process].append(float(depthDict[contig]))

            outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
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

        else:
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)

            if args.element in elements:
                elementDict = defaultdict(lambda: defaultdict(list))
                bits = open(HMMdir + "/hmm-meta.txt", "r")
                for i in bits:
                    ls = (i.rstrip().split("\t"))
                    if ls[0] != "hmm":
                        if ls[3] == args.element:
                            elementDict[ls[3]][ls[2]].append(ls[0])

                catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                for i in elementDict.keys():
                    for j in elementDict[i]:
                        for k in elementDict[i][j]:
                            catDict[k] = j

                Dict = defaultdict(lambda: defaultdict(list))
                final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
                for i in final:
                    ls = (i.rstrip().split(","))
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                        if not re.match(r'#', i):
                            element = ls[0]
                            cell = ls[3]
                            orf = ls[4]
                            contig = allButTheLast(orf, "_")
                            gene = ls[2]
                            reaction = ls[1]
                            Dict[cell][gene].append(float(depthDict[contig]))

                cats = []
                CatOrgDict = defaultdict(list)
                for i in Dict.keys():
                    for j in Dict[i]:
                        if len(catDict[j]) > 0:
                            if j not in cats:
                                cats.append(j)
                                CatOrgDict[catDict[j]].append(j)

                outHeat = open("%s/%s.%s.readDepth.heatmap.csv" % (args.outdir, args.out, args.element), "w")
                outHeat.write("X" + ',')
                for i in sorted(Dict.keys()):
                    outHeat.write(i + ",")
                outHeat.write("\n")

                for i in CatOrgDict.keys():
                    for j in CatOrgDict[i]:
                        outHeat.write(i + "--" + j + ",")
                        for k in sorted(Dict.keys()):
                            outHeat.write(str(SUM(Dict[k][j])) + ",")
                        outHeat.write("\n")

                outHeat.close()
                print('......')
                print(".......")
                print("Finished!")

            else:
                print("Looks like the element you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling.")

    # GENE COUNTS-BASED ABUNDANCE
    else:
        if args.element == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                    if not re.match(r'#', i):
                        process = ls[0]
                        cell = ls[3]
                        orf = ls[4]
                        gene = ls[2]
                        Dict[cell][process].append(gene)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in os.listdir(args.bin_dir):
                if lastItem(i.split(".")) == args.bin_ext:
                    file = open("%s/%s-proteins.faa" % (args.bin_dir, i), "r")
                    file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/%s.%s.heatmap.csv" % (args.outdir, args.out, args.element), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
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
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)
            # elements = ["sulfur", "hydrogen", "methane", "nitrogen", "oxygen",
            #             "carbon-monoxide", "C1compounds", "carbon",
            #             "urea", "halogenetated-compounds", "arsenic", "selenium",
            #             "nitriles", "iron"]

            if args.element in elements:
                elementDict = defaultdict(lambda: defaultdict(list))
                bits = open(HMMdir + "/hmm-meta.txt", "r")
                for i in bits:
                    ls = (i.rstrip().split("\t"))
                    if ls[0] != "hmm":
                        if ls[3] == args.element:
                            elementDict[ls[3]][ls[2]].append(ls[0])

                catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                for i in elementDict.keys():
                    for j in elementDict[i]:
                        for k in elementDict[i][j]:
                            catDict[k] = j

                Dict = defaultdict(lambda: defaultdict(list))
                final = open("%s/%s-summary.csv" % (args.outdir, args.out), "r")
                for i in final:
                    ls = (i.rstrip().split(","))
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
                        if not re.match(r'#', i):
                            element = ls[0]
                            cell = ls[3]
                            orf = ls[4]
                            gene = ls[2]
                            reaction = ls[1]
                            Dict[cell][gene].append(orf)

                cats = []
                CatOrgDict = defaultdict(list)
                for i in Dict.keys():
                    for j in Dict[i]:
                        if len(catDict[j]) > 0:
                            if j not in cats:
                                cats.append(j)
                                CatOrgDict[catDict[j]].append(j)

                normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                for i in binDirLS:
                    if lastItem(i.split(".")) == args.bin_ext:
                        file = open("%s/%s-proteins.faa" % (binDir, i), "r")
                        file = fasta(file)
                        normDict[i] = len(file.keys())

                outHeat = open("%s/%s.%s.heatmap.csv" % (args.outdir, args.out, args.element), "w")
                outHeat.write("X" + ',')
                for i in sorted(Dict.keys()):
                    outHeat.write(i + ",")
                outHeat.write("\n")

                for i in CatOrgDict.keys():
                    for j in CatOrgDict[i]:
                        outHeat.write(i + "--" + j + ",")
                        for k in sorted(Dict.keys()):
                            if not re.match(r'#', j):
                                if args.norm == "y":
                                    outHeat.write(str((len(Dict[k][j]) / int(normDict[k])) * float(100)) + ",")
                                else:
                                    outHeat.write(str((len(Dict[k][j]))) + ",")
                        outHeat.write("\n")

                outHeat.close()
                print('......')
                print(".......")
                print("Finished!")
            else:
                print("Looks like the element you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling.")

# ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
if args.makeplots == "y":
    if conda == 0:
        Rdir = args.R
    else:
        os.system("echo ${rscripts} > r.txt")
        file = open("r.txt")
        for i in file:
            Rdir = (i.rstrip())
        os.system("rm r.txt")

    if args.bams == "NA":
        os.system("Rscript --vanilla %s/DotPlot.R %s/%s.%s.heatmap.csv %s" % (
        Rdir, args.outdir, args.out, args.element, args.outdir))
        os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%s.%s.heatmap.csv %s" % (
        Rdir, args.outdir, args.out, args.element, args.outdir))
    else:
        os.system("Rscript --vanilla %s/DotPlot.R %s/%s.%s.readDepth.heatmap.csv %s" % (
        Rdir, args.outdir, args.out, args.element, args.outdir))
        os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%s.%s.readDepth.heatmap.csv %s" % (
        Rdir, args.outdir, args.out, args.element, args.outdir))


