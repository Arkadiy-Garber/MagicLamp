#!/usr/bin/env python3

from sys import argv, stderr
from lamp import FeGenie, LithoGenie, RosGenie, MagnetoGenie, WspGenie, Lucifer, Genie

"""
MagicLamp.py: A script for querying HMMs against provided datasets and processing output.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.
Installation requirements:
    *python3
    *libraries from the Python standard library: see FeGenie.py, LithoGenie.py, Lucifer, RosGenie, WspGenie, and Genie
 """
__author__ = "Arkadiy Garber"
__version__ = "1"
__maintainer__ = "Arkadiy Garber"
__email__ = "rkdgarber@gmail.com"

errorMessage = "Options: MagicLamp.py [ FeGenie | LithoGenie | RosGenie.py | MagnetoGenie.py | WspGenie | Lucifer | Genie | help ]\n"

try:
    argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

if argv[1] == "FeGenie":
    FeGenie.main()
elif argv[1] == "LithoGenie":
    LithoGenie.main()
elif argv[1] == "Lucifer":
    Lucifer.main()
elif argv[1] == "RosGenie":
    RosGenie.main()
elif argv[1] == "WspGenie":
    WspGenie.main()
elif argv[1] == "MagnetoGenie":
    WspGenie.main()
elif argv[1] == "Genie":
    WspGenie.main()
elif argv[1] == "help":
    stderr.write("\tMagicLamp.py FeGenie: HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenomes.\n"
                 
                 "\tMagicLamp.py LithoGenie: HMM-based identification and categorization of genes and operons relevant to chemolithoautotrophic metabolisms.\n"
                 
                 "\tMagicLamp.py RosGenie: HMM-based identifications of all genes responsible for neutralization of reactive-oxygen species.\n"
                 
                 "\tMagicLamp.py MagnetoGenie: HMM-based identification of genes responsible for magnetosome formation.\n"
                 
                 "\tMagicLamp.py WspGenie: HMM-based identification of the Wsp operon.\n"
                 
                 "\tMagicLamp.py Lucifer: HMM-based identification of light-sensing and light-producing genes\n"
                 
                 "\tMagicLamp.py Genie: Identification of a user-provided set of HMMs\n")

    exit()
else:
    stderr.write(errorMessage)
    exit()

