#!/usr/bin/env python3

from sys import argv, stderr
from genies import FeGenie, LithoGenie, RosGenie, MagnetoGenie, WspGenie, Lucifer, HmmGenie, GasGenie, MnGenie

"""
MagicLamp.py: A script for querying HMMs against provided datasets and processing output.
Tested mostly on .gbk files annotated by Prokka with the --compliant flag.
Installation requirements:
    *python3
    *libraries from the Python standard library: see FeGenie.py, LithoGenie.py, Lucifer, RosGenie, WspGenie, GasGenie, 
    and HmmGenie
 """
__author__ = "Arkadiy Garber"
__version__ = "1"
__maintainer__ = "Arkadiy Garber"
__email__ = "rkdgarber@gmail.com"

errorMessage = "Options: MagicLamp.py [ FeGenie | LithoGenie | RosGenie.py | MagnetoGenie.py | WspGenie | Lucifer | GasGenie | MnGenie | HmmGenie | help ]\n"

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
    MagnetoGenie.main()
elif argv[1] == "GasGenie":
    GasGenie.main()
elif argv[1] == "MnGenie":
    MnGenie.main()
elif argv[1] == "HmmGenie":
    HmmGenie.main()
elif argv[1] == "help":
    stderr.write("\tMagicLamp.py FeGenie: HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenomes.\n"
                 
                 "\tMagicLamp.py LithoGenie: HMM-based identification and categorization of genes and operons relevant to chemolithoautotrophic metabolisms.\n"
                 
                 "\tMagicLamp.py RosGenie: HMM-based identifications of all genes responsible for neutralization of reactive-oxygen species.\n"
                 
                 "\tMagicLamp.py MagnetoGenie: HMM-based identification of genes responsible for magnetosome formation.\n"
                 
                 "\tMagicLamp.py WspGenie: HMM-based identification of the Wsp operon.\n"
                 
                 "\tMagicLamp.py Lucifer: HMM-based identification of light-sensing and light-producing genes\n"
                 
                 "\tMagicLamp.py GasGenie: HMM-based identification of genes responsible for gas vesicle formation.\n"

                 "\tMagicLamp.py MnGenie: HMM-based identification of genes related to manganese transport and oxidation. Also includes genes that are known to bind manganese.\n"
    
                 "\tMagicLamp.py HmmGenie: Identification of a user-provided set of HMMs.\n")

    exit()
else:
    stderr.write(errorMessage)
    exit()

