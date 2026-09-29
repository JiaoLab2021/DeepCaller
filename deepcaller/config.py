"""
Central configuration for DeepCaller.

Every tunable constant that is shared by more than one module lives here, so the
encoder, the network and the command line interface can never drift apart.
"""

import os

import numpy as np

# --------------------------------------------------------------------------- #
# Feature tensor layout
# --------------------------------------------------------------------------- #
# Each candidate locus is represented by a (ploidy, WINDOW_LENGTH, FEATURE_DIM)
# tensor:
#   channel 0        local read depth, normalised by the chromosome mean depth
#   channels 1-7     fraction of reads supporting each pileup token class
#   channels 8-14    mean mapping quality of each token class, scaled by MAPQ_SCALE
WINDOW_SIZE = 10                      # flanking positions kept on each side
WINDOW_LENGTH = WINDOW_SIZE * 2 + 1   # total window length (centre included)
TOKEN_CLASSES = 7                     # A, C, G, T, insertion, deletion, gap/none
FEATURE_DIM = 1 + 2 * TOKEN_CLASSES   # depth + counts + mapping qualities
MAPQ_SCALE = 60.0                     # mapping qualities are normalised by this

# float16 halves the resident size of the encoding tensor. All features are
# ratios or scaled means bounded by ~1, so the reduced mantissa is harmless;
# the network casts back to float32 before the forward pass.
ENCODING_DTYPE = np.float16

# Number of loci held per allocation block while a region is being encoded.
# Larger blocks mean fewer allocations, smaller blocks mean a lower peak.
ENCODING_BLOCK_SIZE = 4096

# Reference context fetched around a region, in base pairs. Must comfortably
# exceed WINDOW_SIZE and the longest deletion that can be reported.
REFERENCE_PADDING = 100

# Allele slots that are always materialised, even when no read supports them.
# Keeping alt1 mandatory guarantees that every locus has at least one active
# head for the network to score.
MANDATORY_GROUPS = frozenset({"alt1"})

# --------------------------------------------------------------------------- #
# Pretrained models
# --------------------------------------------------------------------------- #
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

MODELS_ROOT = os.environ.get(
    "DEEPCALLER_MODELS_ROOT",
    os.path.abspath(os.path.join(_PACKAGE_DIR, "..", "models")),
)

WEIGHTS_FILENAME = "deepcaller.weights.h5"

# ploidy -> {default species alias, alias -> directory under MODELS_ROOT}
# The first alias of each ploidy is the default one.
MODEL_REGISTRY = {
    4: {
        "default": "C88_Potato",
        "choices": {
            "C88_Potato": "01_C88_Potato",
            "Bolivia_Alfalfa": "02_Bolivia_Alfalfa",
            "Samantha_Rose": "03_Samantha_Rose",
        },
    },
    6: {
        "default": "SyntheticPotato_Potato",
        "choices": {
            "SyntheticPotato_Potato": "04_SyntheticPotato_Potato",
            "Tanzania_Sweetpotato": "05_Tanzania_Sweetpotato",
        },
    },
}

SUPPORTED_PLOIDY = tuple(sorted(MODEL_REGISTRY))

# Target mean depth used by --downsample, per ploidy.
DOWNSAMPLE_TARGET_DEPTH = {4: 50, 6: 80}

# --------------------------------------------------------------------------- #
# Pileup filters
# --------------------------------------------------------------------------- #
# Discard unmapped, secondary, QC-failed, duplicate and supplementary records.
PILEUP_FLAG_FILTER = 0x704 | 0x800


def species_choices(ploidy):
    """Return the species aliases available for a given ploidy."""
    return list(MODEL_REGISTRY[ploidy]["choices"])


def default_species(ploidy):
    """Return the species alias used when --species is omitted."""
    return MODEL_REGISTRY[ploidy]["default"]


def describe_species_options():
    """Render the --species help text directly from the registry."""
    parts = ["Species model to use; the default depends on --ploidy."]
    for ploidy in SUPPORTED_PLOIDY:
        entry = MODEL_REGISTRY[ploidy]
        parts.append(
            "ploidy={p}: default={d}, choices={{{c}}}.".format(
                p=ploidy, d=entry["default"], c=", ".join(entry["choices"])
            )
        )
    return " ".join(parts)


def resolve_weights(species, ploidy):
    """Map (ploidy, species alias) to the absolute path of a weights file."""
    directory = MODEL_REGISTRY[ploidy]["choices"][species]
    return os.path.join(MODELS_ROOT, directory, WEIGHTS_FILENAME)
