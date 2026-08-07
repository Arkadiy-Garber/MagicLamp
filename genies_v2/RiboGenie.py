#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# RiboGenie: HMM-based identification and subcategorization of ribosomal & translation machinery profiler
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
        prog="MagicLamp.py RiboGenie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r"""
        *******************************************************

        RiboGenie - part of the MagicLamp suite.
        Developed by Arkadiy Garber, with curated gene categories from
        Ileana Perez-Rodriguez, Abhijit Mukherjee, Maya Mehta, and
        Isabella Viney (Penn Geomicrobiology Laboratory).
        Please send comments and inquiries to ark@midauthorbio.com
                      ____  _ _           ____            _
                     |  _ \(_) |__   ___ / ___| ___ _ __ (_) ___
                     | |_) | | '_ \ / _ \ |  _ / _ \ '_ \| |/ _ \
                     |  _ <| | |_) | (_) | |_| |  __/ | | | |  __/
                     |_| \_\_|_.__/ \___/\____|\___|_| |_|_|\___|
                  Ribosomal & translation machinery profiler

        Single-pass analysis:
          Profiles genomes/proteomes with TIGRFAM HMMs for ribosomal
          proteins, ribosome biogenesis factors, translation factors,
          tRNA biogenesis, and ribosome-associated proteins.

        Each hit is assigned to a functional subcategory defined in
        the curated gene survey, and counts are aggregated per
        subcategory in the heatmap output.

        *******************************************************
        """))

    parser.add_argument('-bin_dir', type=str, help="directory of bins/genomes/assemblies", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic 'cluster'."
                                             " This number should be an integer and should reflect the maximum number "
                                             "of genes in between putative ribosomal & translation machinery profiler identified by the HMM "
                                             "database (default=5)", default=5)

    parser.add_argument('-out', type=str, help="name output directory (default=ribogenie_out)",
                        default="ribogenie_out")

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
    #   Layout: $$ribo_hmms points at a directory containing the .hmm
    #   files directly (flat layout, matching PortGenie convention).
    # =================================================================

    os.system("echo ${ribo_hmms} > HMMlib.txt")
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
            HMMrootdir = location + "/hmms/ribo_hmms"
        except FileNotFoundError:
            pass

    if not have_lib(HMMrootdir):
        print("RiboGenie could not locate its HMM library. Please run the setup.sh script if\n"
              "you have Conda installed. Otherwise, please run setup-noconda.sh and put MagicLamp.py\n"
              "into your $PATH. Set the $ribo_hmms environment variable to point at a directory\n"
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
        "TIGR00001": ('large ribosomal subunit', 'rpmI', 'rpmI'),
        "TIGR00002": ('small ribosomal subunit', 'rpsP', 'rpsP'),
        "TIGR00005": ('ribosome biogenesis', 'rluD', 'rluD'),
        "TIGR00006": ('ribosome biogenesis', 'mraW', 'mraW'),
        "TIGR00008": ('translation initiation factor', 'infA', 'infA'),
        "TIGR00009": ('large ribosomal subunit', 'rpmB', 'rpmB'),
        "TIGR00011": ('tRNA biogenesis', 'ybaK', 'ybaK'),
        "TIGR00012": ('large ribosomal subunit', 'rpmC', 'rpmC'),
        "TIGR00029": ('small ribosomal subunit', 'rpsT', 'rpsT'),
        "TIGR00030": ('small ribosomal subunit', 'rpsU', 'rpsU'),
        "TIGR00035": ('miscellaneous', 'racD', 'racD'),
        "TIGR00038": ('translation elongation factor', 'efp', 'efp'),
        "TIGR00043": ('ribosome biogenesis', 'ybeY', 'ybeY'),
        "TIGR00046": ('ribosome biogenesis', 'rsmE', 'rsmE'),
        "TIGR00048": ('ribosome biogenesis', 'rlmN', 'rlmN'),
        "TIGR00049": ('tRNA biogenesis', 'iscA', 'iscA'),
        "TIGR00050": ('tRNA biogenesis', 'lasT', 'lasT'),
        "TIGR00057": ('ribosome biogenesis', 'sua', 'sua'),
        "TIGR00059": ('large ribosomal subunit', 'rplQ', 'rplQ'),
        "TIGR00060": ('large ribosomal subunit', 'rpIR', 'rpIR'),
        "TIGR00061": ('large ribosomal subunit', 'rplU', 'rplU'),
        "TIGR00062": ('large ribosomal subunit', 'rpmA', 'rpmA'),
        "TIGR00064": ('translation initiation factor', 'ftsY', 'ftsY'),
        "TIGR00071": ('tRNA biogenesis', 'truA', 'truA'),
        "TIGR00082": ('ribosome biogenesis', 'rbfA', 'rbfA'),
        "TIGR00086": ('miscellaneous', 'smpB', 'smpB'),
        "TIGR00088": ('tRNA biogenesis', 'trmD', 'trmD'),
        "TIGR00090": ('ribosome biogenesis', 'rsfS', 'rsfS'),
        "TIGR00091": ('tRNA biogenesis', 'trmB', 'trmB'),
        "TIGR00092": ('ribosome biogenesis', 'ychF', 'ychF'),
        "TIGR00093": ('ribosome biogenesis', 'rluB', 'rluB'),
        "TIGR00094": ('tRNA biogenesis', 'truD', 'truD'),
        "TIGR00095": ('ribosome biogenesis', 'rsmD', 'rsmD'),
        "TIGR00096": ('ribosome biogenesis', 'rsmI', 'rsmI'),
        "TIGR00105": ('large ribosomal subunit', 'rmpE', 'rmpE'),
        "TIGR00113": ('tRNA biogenesis', 'queA', 'queA'),
        "TIGR00116": ('ribosome biogenesis', 'tsf', 'tsf'),
        "TIGR00136": ('tRNA biogenesis', 'gid', 'gid'),
        "TIGR00137": ('tRNA biogenesis', 'trmFO', 'trmFO'),
        "TIGR00138": ('ribosome biogenesis', 'rsmG', 'rsmG'),
        "TIGR00150": ('tRNA biogenesis', 'tsaE', 'tsaE'),
        "TIGR00157": ('ribosome biogenesis', 'rsgA', 'rsgA'),
        "TIGR00158": ('large ribosomal subunit', 'rplI', 'rplI'),
        "TIGR00165": ('small ribosomal subunit', 'rpsR', 'rpsR'),
        "TIGR00166": ('small ribosomal subunit', 'rpsF', 'rpsF'),
        "TIGR00168": ('translation initiation factor', 'infC', 'infC'),
        "TIGR00174": ('tRNA biogenesis', 'miaA', 'miaA'),
        "TIGR00185": ('tRNA biogenesis', 'trmL', 'trmL'),
        "TIGR00186": ('ribosome biogenesis', 'yacO', 'yacO'),
        "TIGR00188": ('tRNA biogenesis', 'rnpA', 'rnpA'),
        "TIGR00202": ('ribosome biogenesis', 'csrA', 'csrA'),
        "TIGR00211": ('tRNA biogenesis', 'glyS', 'glyS'),
        "TIGR00216": ('small ribosomal subunit', 'ispH', 'ispH'),
        "TIGR00233": ('tRNA biogenesis', 'trpS', 'trpS'),
        "TIGR00234": ('tRNA biogenesis', 'tyrS', 'tyrS'),
        "TIGR00246": ('ribosome biogenesis', 'rlmH', 'rlmH'),
        "TIGR00250": ('16S rRNA', 'ruvX', 'ruvX'),
        "TIGR00253": ('ribosome biogenesis', 'yhbY', 'yhbY'),
        "TIGR00269": ('tRNA biogenesis', 'ttuA', 'ttuA'),
        "TIGR00276": ('tRNA biogenesis', 'queG', 'queG'),
        "TIGR00279": ('large ribosomal subunit', 'rpl10e', 'rpl10e'),
        "TIGR00280": ('large ribosomal subunit', 'rpl37ae', 'rpl37ae'),
        "TIGR00285": ('miscellaneous', 'ssh10b', 'ssh10b'),
        "TIGR00291": ('ribosome biogenesis', 'SDO1', 'SDO1'),
        "TIGR00307": ('small ribosomal subunit', 'rps8e', 'rps8e'),
        "TIGR00308": ('tRNA biogenesis', 'trm1', 'trm1'),
        "TIGR00323": ('ribosome biogenesis', 'eif6', 'eif6'),
        "TIGR00324": ('tRNA biogenesis', 'endA', 'endA'),
        "TIGR00329": ('tRNA biogenesis', 'KAE1', 'KAE1'),
        "TIGR00342": ('tRNA biogenesis', 'thiI', 'thiI'),
        "TIGR00344": ('tRNA biogenesis', 'alaS', 'alaS'),
        "TIGR00364": ('tRNA biogenesis', 'queC', 'queC'),
        "TIGR00388": ('tRNA biogenesis', 'glyQ', 'glyQ'),
        "TIGR00392": ('tRNA biogenesis', 'ileS', 'ileS'),
        "TIGR00396": ('tRNA biogenesis', 'leuS', 'leuS'),
        "TIGR00398": ('tRNA biogenesis', 'metG', 'metG'),
        "TIGR00405": ('ribosome biogenesis', 'spt5', 'spt5'),
        "TIGR00406": ('ribosome biogenesis', 'prmA', 'prmA'),
        "TIGR00408": ('tRNA biogenesis', 'proS', 'proS'),
        "TIGR00414": ('tRNA biogenesis', 'serS', 'serS'),
        "TIGR00418": ('tRNA biogenesis', 'thrS', 'thrS'),
        "TIGR00420": ('tRNA biogenesis', 'mnmA', 'mnmA'),
        "TIGR00422": ('tRNA biogenesis', 'valS', 'valS'),
        "TIGR00430": ('tRNA biogenesis', 'tgt', 'tgt'),
        "TIGR00431": ('ribosome biogenesis', 'truB', 'truB'),
        "TIGR00435": ('tRNA biogenesis', 'cysS', 'cysS'),
        "TIGR00436": ('ribosome biogenesis', 'era', 'era'),
        "TIGR00438": ('ribosome biogenesis', 'rrmJ', 'rrmJ'),
        "TIGR00442": ('tRNA biogenesis', 'hisS', 'hisS'),
        "TIGR00447": ('tRNA biogenesis', 'pth', 'pth'),
        "TIGR00449": ('tRNA biogenesis', 'tgtA', 'tgtA'),
        "TIGR00450": ('large ribosomal subunit', 'mnmE', 'mnmE'),
        "TIGR00456": ('tRNA biogenesis', 'argS', 'argS'),
        "TIGR00458": ('tRNA biogenesis', 'aspS', 'aspS'),
        "TIGR00464": ('tRNA biogenesis', 'gltX', 'gltX'),
        "TIGR00468": ('tRNA biogenesis', 'pheS', 'pheS'),
        "TIGR00470": ('tRNA biogenesis', 'sepS', 'sepS'),
        "TIGR00472": ('tRNA biogenesis', 'pheT', 'pheT'),
        "TIGR00478": ('ribosome biogenesis', 'tlyA', 'tlyA'),
        "TIGR00479": ('ribosome biogenesis', 'rumA', 'rumA'),
        "TIGR00485": ('ribosome biogenesis', 'tuf', 'tuf'),
        "TIGR00487": ('translation initiation factor', 'infB', 'infB'),
        "TIGR00490": ('translation elongation factor', 'fusA', 'fusA'),
        "TIGR00496": ('ribosome recycling factor', 'frr', 'frr'),
        "TIGR00499": ('tRNA biogenesis', 'lysS', 'lysS'),
        "TIGR00503": ('translation elongation factor', 'prfC', 'prfC'),
        "TIGR00543": ('miscellaneous', 'menF', 'menF'),
        "TIGR00563": ('ribosome biogenesis', 'rsmB', 'rsmB'),
        "TIGR00581": ('miscellaneous', 'moaC', 'moaC'),
        "TIGR00691": ('ribosome biogenesis', 'relA', 'relA'),
        "TIGR00717": ('small ribosomal subunit', 'rpsA', 'rpsA'),
        "TIGR00731": ('large ribosomal subunit', 'Ctc', 'Ctc'),
        "TIGR00737": ('tRNA biogenesis', 'dusB', 'dusB'),
        "TIGR00741": ('ribosome biogenesis', 'hpf', 'hpf'),
        "TIGR00755": ('ribosome biogenesis', 'ksgA', 'ksgA'),
        "TIGR00757": ('ribosome biogenesis', 'rng', 'rng'),
        "TIGR00855": ('large ribosomal subunit', 'rpIL', 'rpIL'),
        "TIGR00922": ('ribosome biogenesis', 'nusG', 'nusG'),
        "TIGR00952": ('small ribosomal subunit', 'rpsO', 'rpsO'),
        "TIGR00959": ('miscellaneous', 'ffh', 'ffh'),
        "TIGR00981": ('small ribosomal subunit', 'rpsL', 'rpsL'),
        "TIGR00982": ('small ribosomal subunit', 'rps12', 'rps12'),
        "TIGR01009": ('small ribosomal subunit', 'rpsC', 'rpsC'),
        "TIGR01011": ('small ribosomal subunit', 'rpsB', 'rpsB'),
        "TIGR01017": ('small ribosomal subunit', 'rpsD', 'rpsD'),
        "TIGR01020": ('small ribosomal subunit', 'rps5', 'rps5'),
        "TIGR01021": ('small ribosomal subunit', 'rpsE', 'rpsE'),
        "TIGR01022": ('large ribosomal subunit', 'rpmJ', 'rpmJ'),
        "TIGR01023": ('large ribosomal subunit', 'rpmG', 'rpmG'),
        "TIGR01024": ('large ribosomal subunit', 'rplS', 'rplS'),
        "TIGR01025": ('small ribosomal subunit', 'rps19p', 'rps19p'),
        "TIGR01028": ('small ribosomal subunit', 'rps7', 'rps7'),
        "TIGR01029": ('small ribosomal subunit', 'rpsG', 'rpsG'),
        "TIGR01030": ('large ribosomal subunit', 'rpmH', 'rpmH'),
        "TIGR01031": ('large ribosomal subunit', 'rpmF', 'rpmF'),
        "TIGR01032": ('large ribosomal subunit', 'prlT', 'prlT'),
        "TIGR01044": ('large ribosomal subunit', 'rplV', 'rplV'),
        "TIGR01046": ('small ribosomal subunit', 'rps10', 'rps10'),
        "TIGR01049": ('small ribosomal subunit', 'rpsJ', 'rpsJ'),
        "TIGR01050": ('small ribosomal subunit', 'rpsS', 'rpsS'),
        "TIGR01066": ('large ribosomal subunit', 'rplM', 'rplM'),
        "TIGR01067": ('large ribosomal subunit', 'rplN', 'rplN'),
        "TIGR01071": ('large ribosomal subunit', 'rplO', 'rplO'),
        "TIGR01077": ('large ribosomal subunit', 'rpl13', 'rpl13'),
        "TIGR01079": ('large ribosomal subunit', 'rplX', 'rplX'),
        "TIGR01125": ('ribosome biogenesis', 'rimO', 'rimO'),
        "TIGR01164": ('large ribosomal subunit', 'rplP', 'rplP'),
        "TIGR01169": ('large ribosomal subunit', 'rplA', 'rplA'),
        "TIGR01171": ('large ribosomal subunit', 'rplB', 'rplB'),
        "TIGR01177": ('tRNA biogenesis', 'trm-G10', 'trm-G10'),
        "TIGR01211": ('tRNA biogenesis', 'ELP3', 'ELP3'),
        "TIGR01213": ('tRNA biogenesis', 'pus10', 'pus10'),
        "TIGR01308": ('large ribosomal subunit', 'rpmD', 'rpmD'),
        "TIGR01309": ('large ribosomal subunit', 'rpl30', 'rpl30'),
        "TIGR01331": ('tRNA biogenesis', 'cysQ', 'cysQ'),
        "TIGR01354": ('ribosome biogenesis', 'cdd', 'cdd'),
        "TIGR01380": ('miscellaneous', 'gshB', 'gshB'),
        "TIGR01388": ('tRNA biogenesis', 'rnd', 'rnd'),
        "TIGR01393": ('translation elongation factor', 'lepA', 'lepA'),
        "TIGR01394": ('translation elongation factor', 'typA', 'typA'),
        "TIGR01574": ('tRNA biogenesis', 'miaB', 'miaB'),
        "TIGR01575": ('ribosome biogenesis', 'rimI', 'rimI'),
        "TIGR01579": ('tRNA biogenesis', 'mtaB', 'mtaB'),
        "TIGR01632": ('large ribosomal subunit', 'rplK', 'rplK'),
        "TIGR01683": ('tRNA biogenesis', 'thiS', 'thiS'),
        "TIGR01902": ('tRNA biogenesis', 'lysK', 'lysK'),
        "TIGR01951": ('ribosome biogenesis', 'nusB', 'nusB'),
        "TIGR01953": ('ribosome biogenesis', 'nusA', 'nusA'),
        "TIGR01966": ('tRNA biogenesis', 'rph', 'rph'),
        "TIGR02006": ('tRNA biogenesis', 'iscS', 'iscS'),
        "TIGR02013": ('miscellaneous', 'rpoB', 'rpoB'),
        "TIGR02021": ('miscellaneous', 'bchM', 'bchM'),
        "TIGR02027": ('large ribosomal subunit', 'rpoA', 'rpoA'),
        "TIGR02034": ('translation elongation factor', 'cysNC', 'cysNC'),
        "TIGR02063": ('tRNA biogenesis', 'vacB', 'vacB'),
        "TIGR02191": ('ribosome biogenesis', 'rnc', 'rnc'),
        "TIGR02258": ('tRNA biogenesis', 'thpR', 'thpR'),
        "TIGR02273": ('small ribosomal subunit', 'rimM', 'rimM'),
        "TIGR02386": ('miscellaneous', 'rpoC', 'rpoC'),
        "TIGR02393": ('translation initiation factor', 'rpoD', 'rpoD'),
        "TIGR02395": ('translation initiation factor', 'rpoN', 'rpoN'),
        "TIGR02420": ('ribosome biogenesis', 'dksA', 'dksA'),
        "TIGR02433": ('tRNA biogenesis', 'tilS', 'tilS'),
        "TIGR02469": ('miscellaneous', 'cbiT', 'cbiT'),
        "TIGR02539": ('tRNA biogenesis', 'sepcysS', 'sepcysS'),
        "TIGR02651": ('tRNA biogenesis', 'rnz', 'rnz'),
        "TIGR02692": ('tRNA biogenesis', 'cca', 'cca'),
        "TIGR02729": ('ribosome biogenesis', 'obg', 'obg'),
        "TIGR03073": ('tRNA biogenesis', 'rtcB', 'rtcB'),
        "TIGR03139": ('tRNA biogenesis', 'queF', 'queF'),
        "TIGR03156": ('ribosome biogenesis', 'hflX', 'hflX'),
        "TIGR03317": ('tRNA biogenesis', 'ygfZ', 'ygfZ'),
        "TIGR03342": ('tRNA biogenesis', 'tusE', 'tusE'),
        "TIGR03365": ('tRNA biogenesis', 'queE', 'queE'),
        "TIGR03367": ('tRNA biogenesis', 'queD', 'queD'),
        "TIGR03591": ('tRNA biogenesis', 'pnp', 'pnp'),
        "TIGR03594": ('ribosome biogenesis', 'engA', 'engA'),
        "TIGR03596": ('ribosome biogenesis', 'ylqF', 'ylqF'),
        "TIGR03598": ('large ribosomal subunit', 'engB', 'engB'),
        "TIGR03625": ('large ribosomal subunit', 'rplC', 'rplC'),
        "TIGR03626": ('large ribosomal subunit', 'rpl3', 'rpl3'),
        "TIGR03627": ('small ribosomal subunit', 'rps9', 'rps9'),
        "TIGR03630": ('small ribosomal subunit', 'rps17', 'rps17'),
        "TIGR03631": ('small ribosomal subunit', 'rpsM', 'rpsM'),
        "TIGR03632": ('small ribosomal subunit', 'rpsK', 'rpsK'),
        "TIGR03635": ('small ribosomal subunit', 'rpsQ', 'rpsQ'),
        "TIGR03636": ('large ribosomal subunit', 'rpl23', 'rpl23'),
        "TIGR03654": ('large ribosomal subunit', 'rplF', 'rplF'),
        "TIGR03677": ('ribosome biogenesis', 'rpl7ae', 'rpl7ae'),
        "TIGR03685": ('large ribosomal subunit', 'rpl12', 'rpl12'),
        "TIGR03722": ('tRNA biogenesis', 'kae1', 'kae1'),
        "TIGR03723": ('tRNA biogenesis', 'tsaD', 'tsaD'),
        "TIGR03724": ('tRNA biogenesis', 'TP53RK', 'TP53RK'),
        "TIGR03838": ('tRNA biogenesis', 'gluQ', 'gluQ'),
        "TIGR03953": ('large ribosomal subunit', 'rplD', 'rplD'),
        "TIGR03972": ('tRNA biogenesis', 'taw1', 'taw1'),
        "TIGR04170": ('tRNA biogenesis', 'rnr', 'rnr'),
    }

    SUBCATEGORIES = [
        '16S rRNA',
        'large ribosomal subunit',
        'miscellaneous',
        'ribosome biogenesis',
        'ribosome recycling factor',
        'small ribosomal subunit',
        'tRNA biogenesis',
        'translation elongation factor',
        'translation initiation factor',
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
        print("Looks like you did not provide an extension for your genomes/bins or assemblies, so RiboGenie does "
              "not know which files in the provided directory are FASTA files that you would like analyzed.")
        print("Exiting")
        raise SystemExit

    try:
        os.listdir(args.out)
        print("Looks like you already have a directory with the name: " + args.out)
        # answer = input("Would you like RiboGenie to proceed and potentially overwrite files in this directory? (y/n): ")
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
                                print("Unfortunately, this is a problem for RiboGenie because it uses that character "
                                      "as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                except FileNotFoundError:
                    binFile = open("%s/%s" % (binDir, i), "r")
                    for line in binFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \\|")
                                print("Unfortunately, this is a problem for RiboGenie because it uses that character "
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
    print("Profiling ribosomal & translation machinery profiler with TIGRFAM HMMs (hmmsearch --cut_tc)")

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
        if ls[0] in ("cell", "RiboGenie"):
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
    out = open("%s/ribogenie-summary.csv" % args.out, "w")
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
    final = open("%s/ribogenie-summary.csv" % args.out, "r")
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

    outHeat = open("%s/ribogenie.heatmap.csv" % outDirectory, "w")
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
    print("Results are written to %s/ribogenie-summary.csv and %s/ribogenie.heatmap.csv" % (args.out, args.out))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()
