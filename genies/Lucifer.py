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

           _-_-                |/\/\             
         /,                 '  ||              
         ||     \\ \\  _-_ \\ =||=  _-_  ,._-_            @
        ~||     || || ||   ||  ||  || \\  ||              /@*   @
         ||     || || ||   ||  ||  ||/    ||        @      &   @@
        (  -__, \\/\\ \\,/ \\  \\, \\,/   \\,       @@     %   (%/
                                                    &%,   /#   %&                   
                                                   /#    .&   *.                   
                                                    ,   ..    ,                    
                                                   ,*   *,,   *.                   
                           ....                    ,,   ***/  ./                   
                        ...                        .,* ..... /.,(                  
                         ..             .            ,..../*,,                     
                          ...             ..           (.....                      
                           ....             .          /./..,                      
                          ., ...,*,     ......          ,*,.                       
                           ,.*%,,.,.........            ./*.                       
                            .,*.,,.....                  ,*                        
                            ..,. ,                        *.,                      
                           ,,.,.,,                      , ...                      
                          .....*,.*.  ,*,               , ,                        
                   ,.**,,....,,....,,**** ..         /(#* /.                       
                   ,**.,,,,,,.,.,,,*,.....,,,.     #####/./.                       
                 .. **,.*************..,,...,**,***..///*,/                        
               ..,***, ..,****.***,. .., % ,,,,.......   ./.                       
               ,,**.. ..,,*...,,.**,,..                  ./.                       
             ,*. .*   /.,..**,***.,*..                   ./.                       
            .,.,*,,     ,..***.**.,*.                     /.                       
              ,.,.#     .,.,,,.,*.,,                     ./.                       
               ,//(/,.  ,,*,***,.*,.                     ./.                       
               .*//##** *.,,,,,, ,,                      ./.                       
                  /.,,    ,,*,*,,.,                      ./.                       
                   ,,,. . ...... .                       ./.                       
                 ,.., ..  .......... .                   ./.                       
                      . .  ...... ... .                  ./.                       
                    , . . ..  ... ... ..                  /.                       
                      ..... ..  ..........               ./.                       
                   .. .................. . .             ./.                       
                   ........... ........ .... .           ./                        
                   .. .......... ....... ..... .         ./                        
                 ,.. . .......... ... ... .. ... ..       .                        
                ,,.. ................. ... ...........*, .* ,                      
              ,, . .. . ............ ...... ... ... ,.... /                        
             *,    .. .. ................... .... , .... ./                        
             ,#  . ...... .......  ............ , ......,./ ...                    
             *   . ... ... ....... ...... ....,.......... .,,* ..                  
             ,*  . .... ........... . .......... .. ... . / .. ., .                
              .* .......... ..............*....... .. ... / ..... , .              
                .*...... ... ......... *,.............(  */    . ,  , ,            
                 . .,,**  ... .. ,*,,.. ..   .... .....  //.     . ,**,*           
                   ..  ..         /   /  . ,,   ... ....../*        .*,(           
                 .**,.*                  ..,*.,,,. ...               ..            
                 ,.**.@                    ( , .             

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

    parser.add_argument('-out', type=str, help="name output directory (default=lucifer_out)",
                        default="lucifer_out")

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
                             "normalization, Lucifer will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, Lucifer will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--makeplots', type=str,
                        help="include this flag if you would like Lucifer to make some figures from your data?. "
                             "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Please be sure to install all the required R packages as instrcuted in the Lucifer Wiki: "
                             "https://github.com/Arkadiy-Garber/Lucifer/wiki/Installation. "
                             "If you see error or warning messages associated with Rscript, you can still expect to "
                             "see the main output (CSV files) from Lucifer.", const=True, nargs="?")

    # CHECKING FOR CONDA INSTALL
    os.system("echo ${lux_hmms} > HMMlib.txt")
    os.system("echo ${rscripts} > rscripts.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = i.rstrip()
    bits = HMMdir + "/" + "bitscores.txt"

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

        HMMdir = location + "/hmms/lux/"
        bits = HMMdir + "/" + "bitscores.txt"
        rscriptDir = location + "/rscripts/"

        try:
            test = open(bits)
        except FileNotFoundError:
            print("MagicLamp script could not locate the required directories. Please run the setup.sh script if \n"
                  "you have Conda installed. Otherwise, please run the setupe-noconda.sh script and put MagicLamp.py \n"
                  "into your $PATH")
            raise SystemExit

    os.system("rm HMMlib.txt rscripts.txt mainDir.txt")

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
            'Looks like you did not provide an extension for your genomes/bins or assemblies, so Lucifer does not know'
            ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)

        answer = input("Would you like Lucifer to proceed and potentially overwrite files in this directory? (y/n): ")
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
                                    "Unfortunately, this is a problem for Lucifer because it uses that character as delimiter to store important information.")
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
                                        "Unfortunately, this is a problem for Lucifer because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                    except FileNotFoundError:
                        binFile = open("%s/%s" % (binDir, i), "r")
                        for line in binFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for Lucifer because it uses that character as delimiter to store important information.")
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
    for i in meta:
        ls = i.rstrip().split("\t")
        metaDict[ls[0]] = ls[1]

    print("")
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
                count += 1
                perc = (count / len(HMMdirLS)) * 100
                sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                sys.stdout.flush()
                if len(metaDict[hmm.split(".")[0]]) == 0:
                    bit = 0
                else:
                    bit = metaDict[hmm.split(".")[0]]

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
                                else:
                                    # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                                    if bit > HMMdict[i][orf]["bit"]:
                                        HMMdict[i][orf]["hmm"] = hmm
                                        HMMdict[i][orf]["evalue"] = evalue
                                        HMMdict[i][orf]["bit"] = bit

            print("")

    out = open("%s/summary.csv" % (outDirectory), "w")
    out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "\n")
    for key in HMMdict.keys():
        for j in HMMdict[key]:
            out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                      str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "\n")

    out.close()
    # ****************************************** DEREPLICATION *********************************************************
    summary = open(outDirectory + "/summary.csv", "r")
    SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] != "category" and ls[0] != "Lucifer":
            if len(ls) > 0:
                category = ls[0]
                cell = ls[0]
                orf = ls[1]
                hmm = ls[2]
                evalue = ls[3]
                hmmBit = ls[4]
                bitcut = metaDict[hmm.split(".")[0]]
                seq = BinDict[cell][orf]

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

                    out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1

                else:
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)

                        out.write(i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + SummaryDict[i][orf]["e"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(SummaryDict[i][orf]["bitcut"]) + "," + str(counter) + "," + str(SummaryDict[i][orf]["seq"]) + "\n")
                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1
    out.close()

# ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************
    clusterDict = defaultdict(lambda: defaultdict(list))
    summary = open("%s/summary-2.csv" % (args.out), "r")
    for i in summary:
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            clu = int(ls[6])
            clusterDict[clu]["line"].append(ls)
            clusterDict[clu]["gene"].append(ls[2].split(".")[0])

    print("..")
    print("...")
    out = open("%s/summary-3.csv" % (args.out), "w")
    out.write("file" + "," + "ORF" + "," "gene" + "," "evalue" + "," "bit_score" + "," "bit_score_cutoff" + "," "cluster_id" + "," "seq" + "\n")
    for i in sorted(clusterDict.keys()):
        ls = (clusterDict[i]["gene"])

        if "LuxAB" in ls or "LuxC" in ls or "LuxD" in ls or "LuxE" in ls:

            lux = ["LuxAB", "LuxC", "LuxD", "LuxE"]

            if unique(ls, lux) < 3:
                pass

            else:
                for j in clusterDict[i]["line"]:
                    out.write(
                        j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[7] + "\n")
                out.write("####################################################" + "\n")

        else:
            for j in clusterDict[i]["line"]:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "," + j[6] + "," + j[7] + "\n")
            out.write("####################################################" + "\n")

    out.close()

    os.system("rm %s/summary.csv" % (args.out))
    os.system("rm %s/summary-2.csv" % (args.out))
    os.system("mv %s/summary-3.csv %s/lucifer-summary.csv" % (args.out, args.out))

    os.system("mkdir -p %s/HMM_results" % outDirectory)
    os.system("rm -f %s/ORF_calls/*-prodigal.out" % outDirectory)
    os.system("rm -rf %s/HMM_results/*-HMM" % outDirectory)
    os.system("mv -f %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

# ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    print("....")
    print(".....")
    # COVERAGE-BASED ABUNDANCE
    if args.bams != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
        BAMmap = open(args.bams)
        for i in BAMmap:
            string = ''
            ls = i.rstrip().split("\t")
            cell = ls[0]
            for j in ls[1:]:
                string += " "
                string += j
            try:
                depth = open("%s/%s.depth" % (args.out, cell))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (args.out, cell, string))
                print("processing... " + cell)
                depth = open("%s/%s.depth" % (args.out, cell))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[cell][LS[0]] = LS[2]

        os.system("mkdir %s/contigDepths" % args.out)
        os.system("mv %s/*depth %s/contigDepths/" % (args.out, args.out))

        cats = ["Proteorhodopsin", "Xantharhodopsin", "Heliorhodopsin", "Heliorhodopsin_Pfam", "Bac_rhodopsin",
                "GpcrRhopsn4", "Htr2", "7tm_1", "ASRT", "Photo_RC", "PsaA_PsaB", "PSII", "PCP", "Chloroa_b-bind",
                "Bac_chlorC", "BChl_A", "PsbH", "LuxAB", "LuxC", "LuxD", "LuxE"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/lucifer-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[cell][contig]))

        outHeat = open("%s/lucifer.readDepth.heatmap.csv" % (args.out), "w")
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

        cats = ["Proteorhodopsin", "Xantharhodopsin", "Heliorhodopsin", "Heliorhodopsin_Pfam", "Bac_rhodopsin",
                "GpcrRhopsn4", "Htr2", "7tm_1", "ASRT", "Photo_RC", "PsaA_PsaB", "PSII", "PCP", "Chloroa_b-bind",
                "Bac_chlorC", "BChl_A", "PsbH", "LuxAB", "LuxC", "LuxD", "LuxE"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/lucifer-summary.csv" % (args.out), "r")
        for i in final:
            ls = (i.rstrip().split(","))
            if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                if not re.match(r'#', i):
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2]
                    contig = allButTheLast(orf, "_")
                    Dict[cell][gene].append(float(depthDict[contig]))

        outHeat = open("%s/lucifer.readDepth.heatmap.csv" % (args.out), "w")
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

    # GENE COUNTS-BASED ABUNDANCE
    else:
        cats = ["Proteorhodopsin", "Xantharhodopsin", "Heliorhodopsin", "Heliorhodopsin_Pfam", "Bac_rhodopsin",
                "GpcrRhopsn4", "Htr2", "7tm_1", "ASRT", "Photo_RC", "PsaA_PsaB", "PSII", "PCP", "Chloroa_b-bind",
                "Bac_chlorC", "BChl_A", "PsbH", "LuxAB", "LuxC", "LuxD", "LuxE"]

        Dict = defaultdict(lambda: defaultdict(list))
        final = open("%s/lucifer-summary.csv" % (args.out), "r")
        for i in final:
            if not re.match(r'#', i):
                ls = (i.rstrip().split(","))
                if ls[0] != "bin" and ls[1] != "assembly" and ls[1] != "genome" and ls[0] != "file":
                    cell = ls[0]
                    orf = ls[1]
                    gene = ls[2].split(".")[0]
                    Dict[cell][gene].append(gene)

        normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in os.listdir(args.bin_dir):
            if lastItem(i.split(".")) == args.bin_ext:
                if args.orfs:
                    file = open("%s/%s" % (binDir, i), "r")
                    file = fasta(file)
                else:
                    file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    file = fasta(file)
                normDict[i] = len(file.keys())

        outHeat = open("%s/lucifer.heatmap.csv" % (outDirectory), "w")
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


    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots:
        print("Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
              "This will not affect any of the output data that was already created. If you see plots generated, great! "
              "If not, you can plot the data as you wish on your own, or start an issue on Lucifer's GitHub repository\n")

        if args.norm:
            os.system("Rscript --vanilla %s/DotPlot.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
        else:
            os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))

        print("\n\n\n")
        print("...")

    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots:
        print("Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
              "This will not affect any of the output data that was already created. If you see plots generated, great! "
              "If not, you can plot the data as you wish on your own, or start an issue on lucifer's GitHub repository\n")

        if args.bam == "NA" and args.bams == "NA":
            if args.norm:
                os.system("Rscript --vanilla %s/DotPlot.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
                os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            else:
                os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
                os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
        else:
            os.system("Rscript --vanilla %s/DotPlot.R %s/lucifer.readDepth.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/lucifer.readDepth.heatmap.csv %s/" % (rscriptDir, outDirectory, outDirectory))

    print("")
    print("Results are written to %s/lucifer-summary.csv and %s/lucifer.heatmap.csv" % (args.out, args.out))
    if args.makeplots:
        print("Pipeline finished without crashing!!! You might have just seen some errors/warnings from Rscript\n"
          "if you uaed the \'--makeplots\' feature. Don't be alarmed. If you are reading this message, \n"
          "then that means that the main pipeline finished without any issues. Thanks for using :)")
    else:
        print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()











