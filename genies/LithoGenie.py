#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def main():
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

    def checkFe(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_reduction", "iron_oxidation"]:
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

    def check1(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_transport", "iron_aquisition-heme_transport"]:
                    count += 1
        return count

    def check1_2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count

    def check2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"]:
                    count += 1
        return count

    def check3(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count

    def checkReg(ls):
        count = 0
        for i in ls:
            hmm = i.split("|")[0]
            if re.findall(r'aquisition', geneToCatDict[hmm]):
                count += 1
        return count

    def checkMam(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] == "magnetosome_formation":
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
                    seq = ''
                else:
                    header = i[1:]
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
        prog="MagicLamp.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************

        Developed by Arkadiy Garber and Nancy Merino;
        University of Southern California, Earth Sciences
        Please send comments and inquiries to arkadiyg@usc.edu

     .')      .'(  .-,.-.,-.      .'(     .-./(     )\.-.    )\.---.   )\  )\  .'(   )\.---.  
    ( /       \  ) ) ,, ,. (  ,') \  )  ,'     )  ,' ,-,_)  (   ,-._( (  \, /  \  ) (   ,-._( 
     ))       ) (  \( |(  )/ (  '-' (  (  .-, (  (  .   __   \  '-,    ) \ (   ) (   \  '-,   
     )'._.-.  \  )    ) \     ) .-.  )  ) '._\ )  ) '._\ _)   ) ,-`   ( ( \ \  \  )   ) ,-`   
    (       )  ) \    \ (    (  ,  ) \ (  ,   (  (  ,   (    (  ``-.   `.)/  )  ) \  (  ``-.  
     )/,__.'    )/     )/     )/    )/  )/ ._.'   )/'._.'     )..-.(      '.(    )/   )..-.(    
 
                                         .....,/%@@&%#*....              
                                 .........,#%%###%&@@@&&%(.....           
                              ........,,,*@@@% ,,,.  ,(%&%#/,......     
                             .....,,****%@@%..**,.    /&@@@&%*%@@/     
                            ..,,********&@&, .**...     ,(&@@@((##(     
                          ....,,,,*&&%//(%       ...    ((%@&***/#(      .,*..       
                        ...,,,,(%%###%%&&,  .,**///,.    /%#/*,,**/(&%,,%/*....,,,
                  (//(#/*,.*/(***(%%%####&%  ..,**///,.   %%(   .,,,,,(#%&%/,.,**/(((#%%%%,  
                 ,.  .*((((#%@@&&%%##(##/  ........,,,/(((%#/.    *   ***/(####((/*,,***///,.  
                .    .,,. .,,.   ,/(/,      .****/((/,..,,.   ,.  .#. ..       ..    .*///(#%%(.   
            .     .  ..,,*///(/,.   .               .*/*. ,  ..         ..,*////*,..  ,(#%&@@@. 
               .   ..,,**///(/*.         ... .    .*/,.. *.     *@(    ..,,*/((((((*,... *(&@@@& 
               .,    ...,,**///,,.     ,      ..     ,*..   .,,  ,. .*/(##(*.,*//((/,....  (&@@#*
         ,..,*/*     ...,,,***,..      ,*  .,    ,*. .,..   .,. .        .*#&@@#.//*,.....  (&%(,
         *#%&%#.         .....       .  /* ..   .,/,       .,, #%      .,,***,  /&@%.,.  ....  #(**
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
  /##(***                      .                     .,,,,.     .  ,         ,**/**(##(#%@@&&/..,,,
   /##*,*#%&&%(//*,......,*(((####&%%                  ,,,   . ,,            ...  .,,/(##(*,.......
   ,*/*,....,*/(((/***/**,,,,,*(%@@&&%/    .......       ,  .. ..      .,.        ..,.........
   .,,*,,.             ..../*,,*&@@@(.   .,**,,..                   ,*,         .........      
     .,,,,..             ..**/////#%#((.   .,**/**,,,.... ..          .*, //       ......          
    .....                 ..,#@@@@@@&( .,  .,,*******,,,,..     .    ..,*(&&      ....       
      ...,.              .. ,&&&&&@@&%#,  .(%%(/*,,....,,*(##(,        ,,,*//      ...         
        .,,..             .,(%&%#((((,.   ,#%#/*,...     .,/((*     ,@@&&&(.     ...      
        ..,,..     ,,.   ./((%&%#/,....   /##(*,. ..*(#*.    ,/*,. &@@@@@@@@&   ....   
          ..,.,..,(/*    /%&@@@&%/,,,/(*  /,.  /&@@&@@@@@#/,   .*,**,*#@@@@&&*  ...   
          /((#%#(*..   ,*#&@@@&(,/((.      .%@@@&&&@@@@@&&&(/. .,   #@@@/. ..    ,(/  
            @@(,,,......,*/%@@@@#**##/*. //#@@/*//#&@@@@@&.   .(#@@&%##/**/*    ,((/,     
             (/,,,,,,,,...,**(&@@@@&%(/*//%&&@@%#/*,,,,(&@@@@%//**%@@@%/*,,,*##.    
             ...,,,,//,....****,,*/%&@@@@@@@%(///*,,,,**/%%%(/%&&@&%%%#/* ./%@#.  /*.       
              ....,,(#,,,,.......,*((#&&&(//////*****((##//(((##/****,  *%@@%..*#/,.     
               ....../%,,,,,,,,,,..,*/(#(((///////(###((((/*,,.        #%@@@%/,,,.       
                 ....,/&@&/,,,,,,,*#/**,,/#%#(/**/(((//((%#/*((....... ,*#(//****,       
                 ..,,**/#&@@@@&&%#//,,,,*%@@@//((/////(###%&/,,,.*@///*,.,,,,.       
                      .,,*/((((((/**,*,,,,*(%@@@@@@@@#(#%&&&&&%##//%%%///,...*,     
                                  %(?/////////&//%                                                
                                     (%((&@@@#/*.                              
                                       @(((/&@@@#///**                         
                                         #&((///////////////*/@                                
                                                             #*@.                             
                                      ()                   * )//*
                                      <^^>             *     (/*   .
                                     .-""-.                  *)
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

        Image design: Kazuki Takahashi (1996);
        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/
        https://ascii.co.uk/text
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=5)", default=5)

    parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format",
                        default="NA")

    parser.add_argument('-out', type=str, help="name output directory (default=litho_out)",
                        default="lithogenie_out")

    parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                        default=1000)

    parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                             "(default=1, max=16)", default=1)

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
    parser.add_argument('-cat', type=str,
                        help="which element do you want to highlight in the heatmap output? [sulfur, hydrogen, "
                             "methane, nitrogen, oxygen, carbon-monoxide, C1compounds, carbon, urea,"
                             "halogenetated-compounds, arsenic, selenium, nitriles, iron, ALL] (default = ALL)",
                        default="ALL")

    parser.add_argument('--skip', type=str, help="skip the main algorithm and just redo the heatmap with different parameters", const=True,
                        nargs="?")
    
    parser.add_argument('--all_results', type=str,
                        help="report all results, regardless of clustering patterns and operon structure", const=True, nargs="?")

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--orfs', type=str,
                        help="include this flag if you are providing bins as open-reading frames or genes in FASTA amino-acid format",
                        const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each iron gene category to be normalized to "
                             "the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, LithoGenie will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, LithoGenie will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--makeplots', type=str,
                        help="include this flag if you would like LithoGenie to make some figures from your data?. "
                             "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Please be sure to install all the required R packages as instrcuted in the LithoGenie Wiki: "
                             "https://github.com/Arkadiy-Garber/LithoGenie/wiki/Installation. "
                             "If you see error or warning messages associated with Rscript, you can still expect to "
                             "see the main output (CSV files) from LithoGenie.", const=True, nargs="?")

    # CHECKING FOR CONDA INSTALL
    os.system("echo ${litho_hmms} > HMMlib.txt")
    os.system("echo ${rscripts} > rscripts.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = i.rstrip()
    bits = HMMdir + "/" + "hmm-meta.txt"

    file = open("rscripts.txt")
    for i in file:
        rscriptDir = i.rstrip()

    try:
        test = open(bits)

    except FileNotFoundError:
        os.system("which MagicLamp.py > mainDir.txt")

        file = open("mainDir.txt")
        for i in file:
            location = i.rstrip()
        location = allButTheLast(location, "/")

        HMMdir = location + "/hmms/litho/"
        bits = HMMdir + "/" + "hmm-meta.txt"
        rscriptDir = location + "/rscripts/"

        try:
            test = open(bits)
        except FileNotFoundError:
            print("MagicLamp script could not locate the required directories. Please run the setup.sh script if \n"
                  "you have Conda installed. Otherwise, please run the setupe-noconda.sh script and put MagicLamp.py \n"
                  "into your $PATH")
            raise SystemExit

    os.system("rm -f HMMlib.txt rscripts.txt mainDir.txt")

    args = parser.parse_known_args()[0]

    # ************** Checking for the required arguments ******************* #
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

    if args.bam != "NA" and args.bams != "NA":
        print("Please provide only one of the following flags: \'-bam\' or \'-bams\'.")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print(
            'Looks like you did not provide an extension for your genomes/bins or assemblies, so LithoGenie does not know'
            ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)

        answer = input("Would you like LithoGenie to proceed and potentially overwrite files in this directory? (y/n): ")
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
        outDirectory = "%s" % args.out[0:len(args.out)-1]
        outDirectoryLS = os.listdir(outDirectory)
    else:
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)

    print("All required arguments provided!")
    print("")

    # *************** MAKE NR A DIAMOND DB AND READ THE FILE INTO HASH MEMORY ************************ #
    if args.ref != "NA":
        try:
            testFile = open(args.ref + ".dmnd")

        except FileNotFoundError:
            print("Making diamond database out of provided reference file")
            os.system("diamond makedb --in %s -d %s" % (args.ref, args.ref))


    # *************** CALL ORFS FROM BINS AND READ THE ORFS INTO HASH MEMORY ************************ #
    if not args.skip:
        BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in binDirLS:
            if lastItem(i.split(".")) == args.bin_ext:
                cell = i
                if not args.gbk:

                    if args.orfs:
                        testFile = open("%s/%s" % (binDir, i), "r")
                        for line in testFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for LithoGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                    else:
                        try:
                            testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                            print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                            for line in testFile:
                                if re.match(r'>', line):
                                    if re.findall(r'\|]', line):
                                        print(
                                            "Looks like one of your fasta files has a header containing the character: \|")
                                        print(
                                            "Unfortunately, this is a problem for LithoGenie because it uses that character as delimiter to store important information.")
                                        print("Please rename your FASTA file headers")
                                        raise SystemExit

                        except FileNotFoundError:
                            binFile = open("%s/%s" % (binDir, i), "r")
                            for line in binFile:
                                if re.match(r'>', line):
                                    if re.findall(r'\|]', line):
                                        print("Looks like one of your fasta files has a header containing the character: \|")
                                        print(
                                            "Unfortunately, this is a problem for LithoGenie because it uses that character as delimiter to store important information.")
                                        print("Please rename your FASTA file headers")
                                        raise SystemExit

                            print("Finding ORFs for " + cell)
                            if args.meta:
                                os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -p meta -q" % (
                                    binDir, i, outDirectory, i, outDirectory, i))
                            else:
                                os.system(
                                    "prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -q" % (
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
                                locus = (ls)
                                locus = (locus.split("       ")[1])
                                locus = locus.split(" ")[0]
                            if re.findall(r'gene   ', ls):
                                gene = (ls)
                                gene = (gene.split("            ")[1])
                                start = (gene.split("..")[0])
                                end = (gene.split("..")[1])
                                start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                                end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                                altContigName = (locus + "_" + start + "_" + end)

                            if re.findall(r'/locus_tag', ls):
                                locusTag = (ls)
                                locusTag = (locusTag.split("=")[1])
                                locusTag = remove(locusTag, ["\""])
                                counter += 1

                            if counter > 0:
                                gbkDict[locus].append(locusTag)
                                counter = 0
                    else:
                        # print(i)
                        gbk = open("%s/%s" % (binDir, i))
                        for gbkline in gbk:
                            ls = gbkline.rstrip()
                            if re.findall(r'LOCUS', ls):
                                locus = (ls)
                                locus = (locus.split("       ")[1])
                                locus = locus.split(" ")[0]
                            if re.findall(r'gene   ', ls):
                                gene = (ls)
                                gene = (gene.split("            ")[1])
                                start = (gene.split("..")[0])
                                end = (gene.split("..")[1])
                                start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                                end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                                altContigName = (locus + "_" + start + "_" + end)
                                counter += 1

                            if re.findall(r'/locus_tag', ls):
                                locusTag = (ls)
                                locusTag = (locusTag.split("=")[1])
                                locusTag = remove(locusTag, ["\""])

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

                if args.orfs:
                    file = open("%s/%s" % (binDir, i))
                else:
                    file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
                file = fasta(file)
                for j in file.keys():
                    orf = j.split(" # ")[0]
                    BinDict[cell][orf] = file[j]

        # ******************** READ BITSCORE CUT-OFFS INTO HASH MEMORY ****************************** #

        meta = open(bits, "r")
        metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        print("\nreading in HMM bitscore cut-offs...")
        for i in meta:
            ls = i.rstrip().split("\t")
            metaDict[ls[0]]["bit"] = ls[1]
            metaDict[ls[0]]["process"] = ls[2]
            metaDict[ls[0]]["element"] = ls[3]
        print("...")

        # ******************* BEGINNING MAIN ALGORITHM **********************************))))
        print("starting main pipeline...")
        HMMdirLS = os.listdir(HMMdir)
        HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
        for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
            if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                os.system(
                    "mkdir -p " + outDirectory + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

                count = 0
                for hmm in HMMdirLS:  # ITERATING THROUGH ALL THE HMM FILES IN THE HMM DIRECTORY
                    if lastItem(hmm.split(".")) == "hmm":
                        count += 1
                        perc = (count / len(HMMdirLS)) * 100
                        sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc+2))
                        sys.stdout.flush()
                        if len(metaDict[hmm.split(".")[0]]["bit"]) == 0:
                            bit = 0
                        else:
                            bit = metaDict[hmm.split(".")[0]]["bit"]
                            if args.orfs:
                                os.system(
                                    "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/%s"
                                    % (int(args.t), float(bit), outDirectory, i, hmm, outDirectory, i, hmm, HMMdir, hmm, binDir, i)
                                )
                            else:
                                os.system(
                                    "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                                    % (int(args.t), float(bit), outDirectory, i, hmm, outDirectory, i, hmm, HMMdir, hmm, outDirectory, i)
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
                                    hmmName = hmm.split(".")[0]
                                    substrate = metaDict[hmmName]["element"]
                                    reaction = metaDict[hmmName]["process"]
                                    bitcut = metaDict[hmmName]["bit"]
                                    ls = delim(line)
                                    evalue = float(ls[4])
                                    bit = float(ls[5])
                                    orf = ls[0]
                                    if evalue < float(1E-1):  # FILTERING OUT BACKGROUND NOISE
                                        # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS

                                        if orf not in HMMdict[i]:
                                            HMMdict[i][orf]["hmm"] = hmm
                                            HMMdict[i][orf]["evalue"] = evalue
                                            HMMdict[i][orf]["bit"] = bit
                                            HMMdict[i][orf]["substrate"] = substrate
                                            HMMdict[i][orf]["reaction"] = reaction
                                            HMMdict[i][orf]["bitcut"] = bitcut
                                        else:
                                            # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                                            if bit > HMMdict[i][orf]["bit"]:
                                                HMMdict[i][orf]["hmm"] = hmm
                                                HMMdict[i][orf]["evalue"] = evalue
                                                HMMdict[i][orf]["bit"] = bit
                                                HMMdict[i][orf]["substrate"] = substrate
                                                HMMdict[i][orf]["reaction"] = reaction
                                                HMMdict[i][orf]["bitcut"] = bitcut

                print("")

        out = open("%s/summary.csv" % (outDirectory), "w")
        out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "reaction" + "," + "substrate" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "\n")
        for key in HMMdict.keys():
            for j in HMMdict[key]:
                out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," + HMMdict[key][j]["reaction"] + "," + HMMdict[key][j]["substrate"] + "," +
                          str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "," + str(HMMdict[key][j]["bitcut"]) + "\n")

        out.close()
        # ****************************************** DEREPLICATION *********************************************************
        summary = open(outDirectory + "/summary.csv", "r")
        SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[0] != "category" and ls[0] != "LithoGenie":
                if len(ls) > 0:
                    category = ls[0]
                    cell = ls[0]
                    orf = ls[1]
                    hmm = ls[2]
                    reaction = ls[3]
                    substrate = ls[4]
                    evalue = ls[5]
                    hmmBit = ls[6]
                    bitcut = ls[7]
                    seq = BinDict[cell][orf]

                    if cell not in SummaryDict.keys():
                        SummaryDict[cell][orf]["hmm"] = hmm
                        SummaryDict[cell][orf]["hmmBit"] = hmmBit
                        SummaryDict[cell][orf]["category"] = category
                        SummaryDict[cell][orf]["e"] = evalue
                        SummaryDict[cell][orf]["bitcut"] = bitcut
                        SummaryDict[cell][orf]["seq"] = seq
                        SummaryDict[cell][orf]["reaction"] = reaction
                        SummaryDict[cell][orf]["substrate"] = substrate

                    else:
                        if orf not in SummaryDict[cell]:
                            SummaryDict[cell][orf]["hmm"] = hmm
                            SummaryDict[cell][orf]["hmmBit"] = hmmBit
                            SummaryDict[cell][orf]["category"] = category
                            SummaryDict[cell][orf]["e"] = evalue
                            SummaryDict[cell][orf]["bitcut"] = bitcut
                            SummaryDict[cell][orf]["seq"] = seq
                            SummaryDict[cell][orf]["reaction"] = reaction
                            SummaryDict[cell][orf]["substrate"] = substrate

                        else:
                            if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                                SummaryDict[cell][orf]["hmm"] = hmm
                                SummaryDict[cell][orf]["hmmBit"] = hmmBit
                                SummaryDict[cell][orf]["category"] = category
                                SummaryDict[cell][orf]["e"] = evalue
                                SummaryDict[cell][orf]["bitcut"] = bitcut
                                SummaryDict[cell][orf]["seq"] = seq
                                SummaryDict[cell][orf]["reaction"] = reaction
                                SummaryDict[cell][orf]["substrate"] = substrate

        # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
        print("Identifying genomic proximities and putative operons")
        CoordDict = defaultdict(lambda: defaultdict(list))
        for i in SummaryDict.keys():
            if i != "cell":
                for j in SummaryDict[i]:
                    contig = allButTheLast(j, "_")
                    numOrf = lastItem(j.split("_"))
                    CoordDict[i][contig].append(int(numOrf))

        counter = 0
        print("Clustering ORFs...")
        print("")
        out = open(outDirectory + "/summary-2.csv", "w")
        for i in CoordDict.keys():
            print(".")
            for j in CoordDict[i]:
                LS = (CoordDict[i][j])
                clusters = (cluster(LS, args.d))
                for k in clusters:
                    if len(RemoveDuplicates(k)) == 1:
                        orf = j + "_" + str(k[0])

                        out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["reaction"] + "," + SummaryDict[i][orf]["substrate"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                        out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1

                    else:
                        for l in RemoveDuplicates(k):
                            orf = j + "_" + str(l)

                            out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["reaction"] + "," + SummaryDict[i][orf]["substrate"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                        out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1
        out.close()

        # ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************
        clusterDict = defaultdict(lambda: defaultdict(list))
        summary = open("%s/summary-2.csv" % (args.out), "r")
        for i in summary:
            if not re.match(r'#', i):
                ls = i.rstrip().split(",")
                clu = int(ls[8])  ###############################
                clusterDict[clu]["line"].append(ls)
                clusterDict[clu]["gene"].append(ls[2].split(".hm")[0])
                clusterDict[clu]["category"].append(ls[4])

        print("..")
        print("...")
        out = open("%s/summary-3.csv" % (args.out), "w")

        if not args.all_results:
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
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
                                    j[0] + "," + j[1] + "," + j[2] + "," + "iron-oxidation" + "," + j[4] + "," + j[
                                        5] + "," + j[6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            else:
                                out.write(
                                    j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                        6] + "," +
                                    j[7] + "," + j[8] + "," + j[9] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    if "MtrA" in ls and "MtrB_TIGR03509" in ls:
                        for j in clusterDict[i]["line"]:
                            if j[2] in ["MtrA", "MtrB_TIGR03509"]:
                                out.write(
                                    j[0] + "," + j[1] + "," + j[2] + "," + "iron-reduction" + "," + j[4] + "," + j[
                                        5] + "," + j[
                                        6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            else:
                                out.write(
                                    j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                        6] + "," +
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
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + "iron-reduction" + "," + j[4] + "," + j[5] + "\n")
                        else:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + "iron-oxidation" + "," + j[4] + "," + j[5] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                ############################################################################################################
                elif "formate_dehydrogenase_alpha" in ls or "formate_dehydrogenase_beta" in ls or "formate_dehydrogenase_gamma" in ls:
                    formatedehyd = ["formate_dehydrogenase_alpha", "formate_dehydrogenase_beta",
                                    "formate_dehydrogenase_gamma"]

                    if unique(ls, formatedehyd) < 2:
                        if len(remove2(ls, formatedehyd)) < 1:
                            pass
                        else:
                            for j in clusterDict[i]["line"]:
                                if j[2] not in formatedehyd:
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[7] + "," + j[8] + "," + j[9] + "\n")

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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                    out.write(
                                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                            6] + "," + j[7] + "," + j[8] + "," + j[9] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                                j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," +
                                j[
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
                    if unique(ls, ["dsrE", "dsrF", "dsrH"]) > 1:
                        for j in clusterDict[i]["line"]:
                            if j[2] not in dsr:
                                out.write(
                                    j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                        6] + "," +
                                    j[
                                        7] + "," + j[8] + "," + j[9] + "\n")
                            else:
                                out.write(
                                    "sulfur-oxidation" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[
                                        5] + "," +
                                    j[6] + "," + j[
                                        7] + "," + j[8] + "," + j[9] + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            if j[2] not in dsr:
                                out.write(
                                    j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[
                                        6] + "," +
                                    j[
                                        7] + "," + j[8] + "," + j[9] + "\n")
                            else:
                                out.write(
                                    "sulfite-reduction" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[
                                        5] + "," +
                                    j[6] + "," + j[
                                        7] + "," + j[8] + "," + j[9] + "\n")

                ############################################################################################################
                else:
                    linels = (clusterDict[i]["line"])
                    for j in linels:
                        out.write(
                            j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                                7] + "," + j[8] + "," + j[9] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        else:
            for i in sorted(clusterDict.keys()):
                linels = (clusterDict[i]["line"])
                for j in linels:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[
                            7] + "," + j[8] + "," + j[9] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        out.close()

        os.system("rm %s/summary.csv" % (args.out))
        os.system("rm %s/summary-2.csv" % (args.out))
        os.system("mv %s/summary-3.csv %s/lithogenie-summary.csv" % (args.out, args.out))

        # ****************************** FILTERING OUT FALSE POSITIVES *************************************
        out = open("%s/lithogenie-summary-fixed.csv" % (args.out), "w")
        out.write(
            "cell" + "," + "ORF" + "," + "HMM" + "," + "reaction" + "," + "substrate" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "\n")
        clusterDict = defaultdict(list)
        memoryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
        summary = open("%s/lithogenie-summary.csv" % (args.out))
        for i in summary:
            if not re.match(r'#', i):
                try:
                    ls = i.rstrip().split(",")
                    dataset = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    cat = ls[3]
                    element = ls[4]
                    e = ls[5]
                    bit = ls[6]
                    cutoff = ls[7]
                    clu = ls[8]
                    seq = ls[9]
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

        if not args.all_results:
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
                            out.write(dataset + "," + orf + "," + memoryDict[dataset][orf]["gene"] + "," + memoryDict[dataset][orf]["cat"] + "," +
                                      memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["e"] + "," + memoryDict[dataset][orf]["bit"] + "," +
                                      memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["seq"] + "\n")

                    elif memoryDict[dataset][orf]["gene"] == "t4ap":
                        aromatics = ["F", "Y", "W", "H"]
                        seq = memoryDict[dataset][orf]["seq"]
                        aromaticAAs = 0
                        aromaticFreeGap = 0
                        gapThreshold = 0
                        for aa in seq:
                            if aa in aromatics:
                                if aromaticFreeGap > 35:
                                    gapThreshold += 1
                                aromaticAAs += 1
                                aromaticFreeGap = 0
                            else:
                                aromaticFreeGap += 1

                        percAromatic = aromaticAAs / len(seq)
                        if percAromatic > 0.097 and gapThreshold == 0 and \
                                re.findall(r'[FYWH](......................)[FYWH](..)[FYWH](....)[FYWH](.................)[FYWH][FYWH](.....)[FYWH]', seq):
                            out.write(
                                dataset + "," + orf + "," + memoryDict[dataset][orf]["gene"] + "," +
                                memoryDict[dataset][orf]["cat"] + "," +
                                memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["e"] + "," +
                                memoryDict[dataset][orf]["bit"] + "," +
                                memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                memoryDict[dataset][orf]["seq"] + "\n")

                    elif memoryDict[dataset][orf]["gene"] in ["GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"]:
                        if checkGACE(clusterDict[i]) < 3:
                            pass
                        else:
                            out.write(
                                dataset + "," + orf + "," + memoryDict[dataset][orf]["gene"] + "," +
                                memoryDict[dataset][orf]["cat"] + "," +
                                memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["e"] + "," +
                                memoryDict[dataset][orf]["bit"] + "," +
                                memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                memoryDict[dataset][orf]["seq"] + "\n")

                    elif memoryDict[dataset][orf]["gene"] in ["DFE_0465", "DFE_0464", "DFE_0463", "DFE_0462", "DFE_0461"]:
                        if checkDFE1(clusterDict[i]) < 3:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf][
                                        "clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                    memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf][
                                        "clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")

                    elif memoryDict[dataset][orf]["gene"] in ["DFE_0451", "DFE_0450", "DFE_0449", "DFE_0448"]:
                        if checkDFE2(clusterDict[i]) < 3:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf][
                                        "clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                    memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf][
                                        "clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")

                    else:
                        out.write(
                            dataset + "," + orf + "," + memoryDict[dataset][orf]["gene"] + "," + memoryDict[dataset][orf]["cat"] + "," +
                            memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["e"] + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
        else:
            for i in clusterDict.keys():
                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                for j in clusterDict[i]:
                    gene = j.split("|")[0]
                    dataset = j.split("|")[1]
                    orf = j.split("|")[2]
                    cat = memoryDict[dataset][orf]["cat"]
                    out.write(dataset + "," + orf + "," + memoryDict[dataset][orf]["gene"] + "," +
                              memoryDict[dataset][orf]["cat"] + "," +
                        memoryDict[dataset][orf]["element"] + "," + memoryDict[dataset][orf]["e"] + "," +
                        memoryDict[dataset][orf]["bit"] + "," +
                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                        memoryDict[dataset][orf]["seq"] + "\n")

        out.close()

        os.system("mv %s/lithogenie-summary-fixed.csv %s/lithogenie-summary.csv" % (args.out, args.out))

    try:
        hmmout = os.listdir("%s/HMM_results" % outDirectory)
        os.system("rm -rf %s/HMM_results/*" % outDirectory)
        os.system("mv %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))
    except FileNotFoundError:
        if not args.skip:
            os.system("mkdir %s/HMM_results" % outDirectory)
            os.system("mv %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

    # if prodigal == 1:
    #     os.system("rm %s/ORF_calls/*-prodigal.out" % outDirectory)


    # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
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
                depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
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
                depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]
                        total += float(LS[2])
                normDict[cell] = total

        if args.cat == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
                        contig = allButTheLast(orf, "_")
                        Dict[cell][process].append(float(depthDict[cell][contig]))

            outHeat = open("%s/lithogenie.%s.readDepth.heatmap.csv" % (args.out, args.cat), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write(str(SUM(Dict[j][i]) / normDict[j]) + ",")
                outHeat.write("\n")

            outHeat.close()
            print('......')
            print(".......")
            print("Finished!")

        else:
            cats = []
            catList = args.cat
            catList = catList.split(",")
            for CAT in catList:
                cats.append(CAT)

            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)

            cats = []
            catList = args.cat
            catList = catList.split(",")
            for CAT in catList:
                if CAT in elements:
                    cats.append(CAT)
                else:
                    print("Looks like an element (%s) that you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling." % CAT)

            elementDict = defaultdict(lambda: defaultdict(list))
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            for i in bits:
                ls = (i.rstrip().split("\t"))
                if ls[0] != "hmm":
                    if ls[3] in cats:
                        elementDict[ls[3]][ls[2]].append(ls[0])

            catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in elementDict.keys():
                for j in elementDict[i]:
                    for k in elementDict[i][j]:
                        catDict[k] = j

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
                        contig = allButTheLast(orf, "_")
                        Dict[cell][gene].append(float(depthDict[cell][contig]))

            cats = []
            CatOrgDict = defaultdict(list)
            for i in Dict.keys():
                for j in Dict[i]:
                    if len(catDict[j]) > 0:
                        if j not in cats:
                            cats.append(j)
                            CatOrgDict[catDict[j]].append(j)

            outHeat = open("%s/lithogenie.%s.readDepth.heatmap.csv" % (args.out, args.cat), "w")
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

        if args.cat == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(cats)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
                        contig = allButTheLast(orf, "_")
                        Dict[cell][process].append(float(depthDict[contig]))

            outHeat = open("%s/lithogenie.%s.readDepth.heatmap.csv" % (args.out, args.cat), "w")
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

            cats = []
            catList = args.cat
            catList = catList.split(",")
            for CAT in catList:
                if CAT in elements:
                    cats.append(CAT)
                else:
                    print("Looks like an element (%s) that you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling." % CAT)

            elementDict = defaultdict(lambda: defaultdict(list))
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            for i in bits:
                ls = (i.rstrip().split("\t"))
                if ls[0] != "hmm":
                    if ls[3] in cats:
                        elementDict[ls[3]][ls[2]].append(ls[0])

            catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in elementDict.keys():
                for j in elementDict[i]:
                    for k in elementDict[i][j]:
                        catDict[k] = j

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
                        contig = allButTheLast(orf, "_")
                        Dict[cell][gene].append(float(depthDict[contig]))

            cats = []
            CatOrgDict = defaultdict(list)
            for i in Dict.keys():
                for j in Dict[i]:
                    if len(catDict[j]) > 0:
                        if j not in cats:
                            cats.append(j)
                            CatOrgDict[catDict[j]].append(j)

            outHeat = open("%s/lithogenie.%s.readDepth.heatmap.csv" % (args.out, args.cat), "w")
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


    # GENE COUNTS-BASED ABUNDANCE
    else:
        if args.cat == "ALL":
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            cats = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in cats:
                    cats.append(element)
            cats = sorted(cats)

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
                        Dict[cell][substrate].append(gene)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in os.listdir(args.bin_dir):
                if lastItem(i.split(".")) == args.bin_ext:
                    if args.orfs:
                        file = open("%s/%s" % (args.bin_dir, i), "r")
                        file = fasta(file)
                    else:
                        file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                        file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/lithogenie.%s.heatmap.csv" % (args.out, args.cat), "w")
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

            print('......')
            print(".......")
            print("Finished!")

        else:
            cats = []
            catList = args.cat
            catList = catList.split(",")
            for CAT in catList:
                cats.append(CAT)

            bits = open(HMMdir + "/hmm-meta.txt", "r")
            elements = []
            for i in bits:
                ls = (i.rstrip().split("\t"))
                element = ls[3]
                if element not in elements:
                    elements.append(element)
            elements = sorted(elements)

            cats = []
            catList = args.cat
            catList = catList.split(",")
            for CAT in catList:
                if CAT in elements:
                    cats.append(CAT)
                else:
                    print("Looks like an element (%s) that you have chosen is not one that is recognized by LithoGenie. Please "
                      "try again by choosing another element, or checking your spelling." % CAT)

            elementDict = defaultdict(lambda: defaultdict(list))
            bits = open(HMMdir + "/hmm-meta.txt", "r")
            for i in bits:
                ls = (i.rstrip().split("\t"))
                if ls[0] != "hmm":
                    if ls[3] in cats:
                        elementDict[ls[3]][ls[2]].append(ls[0])

            catDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in elementDict.keys():
                for j in elementDict[i]:
                    for k in elementDict[i][j]:
                        catDict[k] = j

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/lithogenie-summary.csv" % (args.out), "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if ls[0] != "" and ls[0] != "cell":
                    if not re.match(r'#', i):
                        cell = ls[0]
                        orf = ls[1]
                        gene = ls[2].split(".hm")[0]
                        process = ls[3]
                        substrate = ls[4]
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
                    if args.orfs:
                        file = open("%s/%s" % (args.bin_dir, i), "r")
                        file = fasta(file)
                    else:
                        file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                        file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/lithogenie.%s.heatmap.csv" % (args.out, args.cat), "w")
            outHeat.write("X" + ',')
            for i in sorted(Dict.keys()):
                outHeat.write(i + ",")
            outHeat.write("\n")

            for i in CatOrgDict.keys():
                for j in CatOrgDict[i]:
                    outHeat.write(i + "--" + j + ",")
                    for k in sorted(Dict.keys()):
                        if not re.match(r'#', j):
                            if args.norm:
                                outHeat.write(str((len(Dict[k][j]) / int(normDict[k])) * float(100)) + ",")
                            else:
                                outHeat.write(str((len(Dict[k][j]))) + ",")
                    outHeat.write("\n")

            outHeat.close()
            print('......')
            print(".......")
            print("Finished!")

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
            os.system("Rscript --vanilla %s/DotPlot.R %s/lithogenie.%s.heatmap.csv %s" % (
                Rdir, args.out, args.cat, args.out))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lithogenie.%s.heatmap.csv %s" % (
                Rdir, args.out, args.cat, args.out))
        else:
            os.system("Rscript --vanilla %s/DotPlot.R %s/lithogenie.%s.readDepth.heatmap.csv %s" % (
                Rdir, args.out, args.cat, args.out))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lithogenie.%s.readDepth.heatmap.csv %s" % (
                Rdir, args.out, args.cat, args.out))


if __name__ == '__main__':
    main()




















