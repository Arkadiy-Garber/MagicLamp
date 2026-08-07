#!/usr/bin/env python3

from sys import argv, stderr
from genies_v3 import (
    FeGenie, LithoGenie, HmmGenie, OmniGenie, Lucifer,
    ATPGenie, PortGenie, RnfGenie,
    AbxGenie, MotiliGenie, RiboGenie, ResistiGenie, SporeGenie,
)

"""
MagicLamp.py: A script for querying HMMs against provided datasets and processing output.
Installation requirements:
    *python3
    *libraries from the Python standard library: see FeGenie.py and HmmGenie
 """
__author__ = "Arkadiy Garber"
__version__ = "3"
__maintainer__ = "Arkadiy Garber"
__email__ = "ark@midauthorbio.com"

errorMessage = (
    "Options: MagicLamp.py [ FeGenie | LithoGenie | PortGenie | RnfGenie | Lucifer | ATPGenie |\n"
    "                       AbxGenie | MotiliGenie | RiboGenie | ResistiGenie | SporeGenie |\n"
    "                       OmniGenie | HmmGenie | help ]\n"
)

try:
    argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

if argv[1] == "FeGenie":
    FeGenie.main()
elif argv[1] == "LithoGenie":
    LithoGenie.main()
elif argv[1] == "OmniGenie":
    OmniGenie.main()
elif argv[1] == "HmmGenie":
    HmmGenie.main()
elif argv[1] == "Lucifer":
    Lucifer.main()
elif argv[1] == "ATPGenie":
    ATPGenie.main()
elif argv[1] == "PortGenie":
    PortGenie.main()
elif argv[1] == "RnfGenie":
    RnfGenie.main()
elif argv[1] == "AbxGenie":
    AbxGenie.main()
elif argv[1] == "MotiliGenie":
    MotiliGenie.main()
elif argv[1] == "RiboGenie":
    RiboGenie.main()
elif argv[1] == "ResistiGenie":
    ResistiGenie.main()
elif argv[1] == "SporeGenie":
    SporeGenie.main()
elif argv[1] == "help":
    stderr.write(
        "\tMagicLamp.py FeGenie: HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenomes.\n"
        "\tMagicLamp.py LithoGenie: HMM-based identification and categorization of genes and operons relevant to chemolithoautotrophic metabolisms.\n"
        "\tMagicLamp.py Lucifer: HMM-based identification and categorization of genes and operons relevant to light-sensing and light-producing reactions.\n"
        "\tMagicLamp.py ATPGenie: HMM-based identification and categorization of genes and operons relevant to the ATP synthases.\n"
        "\tMagicLamp.py PortGenie: HMM-based identification and categorization of genes and operons relevant to sodium antiporters and symporters.\n"
        "\tMagicLamp.py RnfGenie: HMM-based identification and categorization of genes and operons relevant to Rnf complex.\n"
        "\tMagicLamp.py AbxGenie: HMM-based identification and subcategorization of antibiotic-biosynthesis genes.\n"
        "\tMagicLamp.py MotiliGenie: HMM-based identification and subcategorization of motility, chemotaxis, thermotaxis, and pili genes.\n"
        "\tMagicLamp.py RiboGenie: HMM-based identification and subcategorization of ribosomal-protein and translation-machinery genes.\n"
        "\tMagicLamp.py ResistiGenie: HMM-based identification and subcategorization of resistance genes (heavy metals, UV, ROS, radioactivity, antibiotics).\n"
        "\tMagicLamp.py SporeGenie: HMM-based identification and subcategorization of sporulation genes by sporulation stage.\n"
        "\tMagicLamp.py OmniGenie: HMM-based identification for a given genie.\n"
        "\tMagicLamp.py HmmGenie: Identification of a user-provided set of HMMs.\n"
    )
    exit()
else:
    stderr.write(errorMessage)
    exit()
