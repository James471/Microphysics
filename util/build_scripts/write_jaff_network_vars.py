#!/usr/bin/env python3

"""Derive SPECIES_ENUM / SPEC_NAMES / NSPEC / NUM_CHEM_BANDS / CHEM_BANDS for a
jaff-generated photochemistry network, so problem CMakeLists.txt files don't
have to hand-duplicate species names/order or radiation band edges that jaff
already emits.

Species names and order come from the network's own _parameters file
(species_N_name entries, written by jaff in its own species order). jaff's
raw names (e.g. "H+", "e-", "c-C3H2", "H2O_DUST") are not valid or not
unambiguous C++ identifiers, so each name is sanitized into one here: the
formula is kept as-is and every structural marker becomes an "_"-separated
suffix ("H+" -> "H_p", "c-C3H2+" -> "C3H2_c_p", "H2O_DUST" -> "H2O_dust"). See
sanitize_species_name() for the full set of rules and for why "_" has to be a
reserved separator rather than the markers simply running onto the formula.

The same sanitized name is used for both SPECIES_ENUM (as the enum tag) and
SPEC_NAMES (as the lookup string), so network_spec_index() and the enum agree
on the species vocabulary by construction instead of by two people typing the
same list twice. SHORT_SPEC_NAMES keeps jaff's raw spelling, so the original
name is always recoverable from the generated header.

aion[] comes from the masses jaff writes alongside those names. zion[] is not
in the _parameters file, so it is derived by parsing each species name into its
constituent atoms and summing their proton numbers (see
species_proton_number()); PROTON_NUMBERS covers every element in jaff's own
mass table, so any network jaff can generate can be configured here.

Radiation band edges come from the network's jaff.toml [radiation].bands
list and are passed through in eV, jaff's native unit for this -- Quokka's
own C++ side (RadSystem::GetChemBandQuanta) is responsible for converting
eV to erg/Hz/whatever it needs, so this script never has to know about or
convert to Quokka's internal units. An open upper edge ("inf" in the TOML,
matching jaff's own sentinel) is passed through as the C++ double infinity
literal so the band array stays a plain fixed-size GpuArray.

[radiation].power_law_index (jaff's photon-number spectrum index alpha, in
n(E) ~ E^(alpha-2)) is passed through unchanged so GetChemBandQuanta can
reproduce jaff's own band-average-energy weighting instead of assuming a
flat spectrum.

Output is written as a small CMake include file that sets NSPEC,
SPECIES_ENUM, SPEC_NAMES, NUM_CHEM_BANDS, CHEM_BANDS, and
POWER_LAW_INDEX, meant to be include()-d by a problem's CMakeLists.txt
before configure_file() of network_header.template.
"""

import argparse
import re
from pathlib import Path

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None

VALID_CXX_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

# jaff's structural name decorations, stripped/rewritten by
# sanitize_species_name() and species_proton_number(). "_" is reserved as the
# separator between a formula and these markers, so a marker can never be
# confused with formula content (see sanitize_species_name).
#
# Isomer prefixes: the cyclic/linear/triplet forms of the same formula.
ISOMER_PREFIXES = ("c-", "l-", "t-")

# Ice-mantle phases of a molecule frozen onto a grain (uclchem).
PHASE_SUFFIXES = ("_BULK", "_DUST")

# jaff's reaction-partner sentinels: not species at all, but emitted as
# reactants/products by the parser. Each has mass 0 in jaff's own atom_mass.csv,
# so a proton number of zero is their definition rather than a fallback.
PSEUDO_SPECIES = frozenset({"_CR", "_CRP", "_CRPHOT", "_PHOTON", "_GRAIN", "_DUMMY"})

# Dust grains are macroscopic particles ("GRAIN0", "GRAIN-"), not molecules, so
# they have no meaningful proton count to contribute to y_e.
GRAIN_PREFIX = "GRAIN"

# KIDA spells a species adsorbed on a grain surface with a leading "X": "XH" is
# H sitting on a surface site (H -> XH, then XH + XH -> H2, the standard surface
# H2 formation route). The following character must be upper-case so real
# elements starting with X -- xenon, "Xe" -- are not mistaken for a prefix.
SURFACE_ADSORBED_RE = re.compile(r"^X[A-Z]")

# 1 atomic mass unit in grams; jaff writes species_N_mass in grams, while
# aion[] is conventionally in amu.
AMU_IN_GRAMS = 1.66053906660e-24

# Proton number Z of every element jaff's mass table knows, so any network it
# can generate can be built here. Kept as a literal table rather than read from
# jaff: this script runs under CMake execute_process() with a bare python3 and
# must work with only the standard library, so it cannot import jaff.
#
# "D" is deuterium, an isotope rather than a distinct element, so it shares
# Z = 1 with H; it is listed because primordial-chemistry networks use it.
PROTON_NUMBERS = {
    "H": 1, "D": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7,
    "O": 8, "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14,
    "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21,
    "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
    "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35,
    "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93,
}

# Element symbols longest-first, so multi-character symbols match before the
# single-character ones they start with ("He" before "H", "Cl" before "C").
# This mirrors how jaff itself decomposes a species name.
_ELEMENTS_LONGEST_FIRST = sorted(PROTON_NUMBERS, key=len, reverse=True)


def parse_species_names(parameters_file):
    pattern = re.compile(r'^species_(\d+)_name\s+string\s+"(.*)"\s*$')
    entries = {}
    for line in parameters_file.read_text().splitlines():
        m = pattern.match(line.strip())
        if m:
            entries[int(m.group(1))] = m.group(2)
    if not entries:
        raise RuntimeError(f"no species_N_name entries found in {parameters_file}")
    return [entries[i] for i in sorted(entries)]


def parse_species_masses_g(parameters_file):
    """Species masses in grams, in jaff's own species order."""
    pattern = re.compile(r"^species_(\d+)_mass\s+real\s+(\S+)\s*$")
    entries = {}
    for line in parameters_file.read_text().splitlines():
        m = pattern.match(line.strip())
        if m:
            entries[int(m.group(1))] = float(m.group(2))
    if not entries:
        raise RuntimeError(f"no species_N_mass entries found in {parameters_file}")
    return [entries[i] for i in sorted(entries)]


def split_species_name(raw_name):
    """Split a jaff species name into its formula and structural decorations.

    Returns (formula, phase, adsorbed, isomer, charge_run), where the formula is
    the bare chemical formula and the rest are the markers stripped off it.
    ``phase`` is "" or one of PHASE_SUFFIXES; ``adsorbed`` is True for a
    surface-adsorbed species; ``isomer`` is "" or one of "c"/"l"/"t";
    ``charge_run`` is the trailing run of "+"/"-" ("" if neutral).

    Both sanitize_species_name() and species_proton_number() decompose names
    through here, so the identifier and the proton number can never disagree
    about what part of a name is formula and what part is decoration.
    """
    name = raw_name

    adsorbed = False
    if SURFACE_ADSORBED_RE.match(name):
        adsorbed = True
        name = name[1:]

    isomer = ""
    for prefix in ISOMER_PREFIXES:
        if name.startswith(prefix):
            isomer = prefix[0]
            name = name[len(prefix):]
            break

    # Taken off before the charge run so the charge parses against the formula
    # itself rather than against the phase suffix's trailing letters.
    phase = ""
    for suffix in PHASE_SUFFIXES:
        if name.endswith(suffix):
            phase = suffix
            name = name[: -len(suffix)]
            break

    m = re.match(r"^(.*?)([+-]*)$", name)
    return m.group(1), phase, adsorbed, isomer, m.group(2)


def species_proton_number(raw_name):
    """Total proton number Z of a species from its jaff name.

    Examples: "H"/"H+" -> 1, "e-" -> 0, "He++" -> 2, "H2O" -> 10, "CH4" -> 10.

    Microphysics' zion[] is a proton count, not a signed charge: composition()
    forms y_e = sum(xn*zion*aion_inv), the electron fraction, which must stay
    non-negative. Using signed charge here would make y_e (and hence mu_e = 1/y_e)
    negative, since aion_inv is ~1823x larger for the electron than for anything
    else. The trailing charge run is therefore parsed off and discarded, and Z is
    summed over the constituent atoms -- for a molecule that sum is what zion[]
    wants, not the Z of any single element.

    The name grammar follows jaff's own species decomposition: a sequence of
    element symbols, each optionally followed by a repeat count, then a trailing
    run of "+"/"-" giving the charge. Symbols are matched longest-first so "He"
    is not read as H followed by unknown "e". A repeat count applies only to the
    symbol immediately before it, so "H2O" is H,H,O rather than (HO),(HO).

    Structural decorations are stripped before decomposition, so a decorated
    species inherits the proton number of its underlying formula: "H2O_DUST"
    and "XH" give Z of "H2O" (10) and "H" (1). Three families carry no protons
    of their own: the electron, jaff's zero-mass reaction-partner sentinels
    (PSEUDO_SPECIES), and dust grains, which are macroscopic particles rather
    than molecules.
    """
    # Sentinels, not species: zero mass in jaff's own table.
    if raw_name in PSEUDO_SPECIES:
        return 0

    base, _phase, _adsorbed, _isomer, _charge = split_species_name(raw_name)

    # The electron is not an element and carries no protons.
    if base == "e":
        return 0

    if base.startswith(GRAIN_PREFIX):
        return 0

    # A trailing "x" marks a lumped species with an unspecified hydrogen count
    # ("CHx", "OHx" in GOW), not an element. jaff itself strips it to recover
    # the base formula (Reaction.serialize), so "CHx" is counted as "CH". The
    # count is by definition approximate for a lumped species; this is the same
    # approximation jaff makes internally.
    if base.endswith("x"):
        base = base[:-1]

    total_z = 0
    pos = 0
    while pos < len(base):
        for symbol in _ELEMENTS_LONGEST_FIRST:
            if base.startswith(symbol, pos):
                break
        else:
            raise RuntimeError(
                f"cannot determine proton number for species {raw_name!r}: "
                f"unrecognised element symbol at {base[pos:]!r}. If this is a "
                "real element, add it to PROTON_NUMBERS."
            )
        pos += len(symbol)

        # An optional repeat count for the symbol just matched.
        count_match = re.match(r"\d+", base[pos:])
        count = 1
        if count_match:
            count = int(count_match.group())
            pos += count_match.end()

        total_z += PROTON_NUMBERS[symbol] * count

    if total_z == 0:
        raise RuntimeError(
            f"cannot determine proton number for species {raw_name!r}: "
            "no element symbols found"
        )
    return total_z


def load_radiation_table(jaff_toml_file):
    """Return the [network.radiation] table from a jaff TOML config.

    Uses tomllib/tomli when available (Python 3.11+, or the tomli backport).
    Neither exists on older interpreters, and Microphysics' other build scripts
    only require Python 3.6, so fall back to a minimal reader that extracts just
    the keys this script needs: the bands array and power_law_index.
    """
    if tomllib is not None:
        with jaff_toml_file.open("rb") as f:
            return tomllib.load(f)["network"]["radiation"]

    table = {}
    in_section = False
    for raw_line in jaff_toml_file.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            in_section = line.strip("[]").strip() in ("network.radiation",)
            continue
        if not in_section or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if value.startswith("["):
            items = value.strip("[]").split(",")
            table[key] = [_parse_scalar(i.strip()) for i in items if i.strip()]
        else:
            table[key] = _parse_scalar(value)

    if not table:
        raise RuntimeError(
            f"{jaff_toml_file}: could not read [network.radiation]; install tomli "
            "(pip install tomli) or use Python 3.11+"
        )
    return table


def _parse_scalar(token):
    token = token.strip()
    if token.startswith(('"', "'")):
        return token[1:-1]
    if token in ("true", "false"):
        return token == "true"
    try:
        return int(token)
    except ValueError:
        return float(token)


def sanitize_species_name(raw_name):
    """Turn a jaff species name into a valid, unambiguous C++ identifier.

    The formula comes first and every structural marker jaff spells with a
    prefix or punctuation becomes an "_"-separated suffix, in the fixed order
    formula + phase + adsorbed + isomer + charge:

      * charge -- a trailing run of "+"/"-" becomes "_p"/"_m", or "_<n>p"/
        "_<n>m" for multiply-charged ions:  "H+" -> "H_p", "He++" -> "He_2p";
      * the electron is spelled "e";
      * isomer prefixes become a trailing "_c"/"_l"/"_t":
        "c-C3H2+" -> "C3H2_c_p";
      * ice phases are lower-cased in place:  "H2O_DUST" -> "H2O_dust";
      * a surface-adsorbed species trades its "X" prefix for "_ads":
        "XH" -> "H_ads";
      * jaff's reaction-partner sentinels drop their leading "_":
        "_CRPHOT" -> "CRPHOT".

    Reserving "_" as the separator is what makes the encoding injective. An
    unseparated suffix runs together with a formula that itself ends in letters
    and digits, so "He++" and "He2+" both collapse onto "He2p", and "N+" gives
    "Np", also the formula for neptunium. Separated, the two readings of the
    digit stay distinct by construction -- "He++" -> "He_2p" against "He2+" ->
    "He2_p" -- rather than by luck. That distinction is not hypothetical: 154 of
    the ~940 species jaff ships end in a digit immediately before their charge
    run ("C2+" -> "C2_p"), and they avoid collision today only because no
    shipped network happens to also contain the doubly-charged twin.

    Markers are stripped by split_species_name(), which species_proton_number()
    also uses, so the identifier and the proton number always agree on which
    part of a name is formula. The formula itself is passed through unchanged
    and must already be identifier-safe; main() validates that, and rejects
    duplicate sanitized names -- register_microphysics_network() turns this
    script's nonzero exit into a CMake FATAL_ERROR.
    """
    if raw_name == "e-":
        return "e"

    if raw_name in PSEUDO_SPECIES:
        return raw_name.lstrip("_")

    base, phase, adsorbed, isomer, charge_run = split_species_name(raw_name)

    # Lower-cased so the marker reads as a suffix rather than as part of the
    # formula, matching the casing of every other marker here.
    name = base + phase.lower()

    if adsorbed:
        name += "_ads"

    if isomer:
        name += "_" + isomer

    if charge_run:
        sign = "p" if charge_run[0] == "+" else "m"
        n = len(charge_run)
        name += "_" + (sign if n == 1 else f"{n}{sign}")

    return name


def parse_bands_ev(jaff_toml_file):
    bands_ev = load_radiation_table(jaff_toml_file)["bands"]

    parsed = []
    for i, edge in enumerate(bands_ev):
        if isinstance(edge, str):
            if edge == "inf" and i == len(bands_ev) - 1:
                parsed.append(float("inf"))
                continue
            raise ValueError(
                f"{jaff_toml_file}: non-numeric band edge {edge!r}; only "
                '"inf" as the last edge is supported'
            )
        parsed.append(float(edge))
    return parsed


def parse_power_law_index(jaff_toml_file):
    return float(load_radiation_table(jaff_toml_file)["power_law_index"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network-parameters", required=True, type=Path,
                         help="Path to the network's own _parameters file (species_N_name entries)")
    parser.add_argument("--jaff-toml", required=True, type=Path,
                         help="Path to the network's jaff.toml (for [radiation].bands)")
    parser.add_argument("--output", required=True, type=Path,
                         help="Path to write the generated .cmake include file")
    args = parser.parse_args()

    raw_names = parse_species_names(args.network_parameters)
    sanitized_names = [sanitize_species_name(n) for n in raw_names]

    for raw_name, sanitized_name in zip(raw_names, sanitized_names):
        if not VALID_CXX_IDENTIFIER.match(sanitized_name):
            raise RuntimeError(
                f"sanitized species name {sanitized_name!r} (from raw name "
                f"{raw_name!r}) is not a valid C++ identifier"
            )

    if len(set(sanitized_names)) != len(sanitized_names):
        raise RuntimeError(
            f"sanitized species names collide: {raw_names} -> {sanitized_names}"
        )

    species_enum_lines = [f"{sanitized_names[0]} = 0"] + sanitized_names[1:]
    species_enum = ",\n  ".join(species_enum_lines)

    spec_names = ",\n  ".join(f'\\"{n}\\"' for n in sanitized_names)
    short_spec_names = ",\n  ".join(f'\\"{n}\\"' for n in raw_names)

    # aion/zion are unused by EOS/photoionization itself, which works from
    # spmasses/gammas -- but eos() in interfaces/eos.H unconditionally calls
    # composition(), which forms sum(xn[n] * aion_inv[n]) and inverts it. Emitting
    # these keeps that inversion well-defined; leaving them unset makes
    # configure_file() substitute empty (i.e. all-zero) arrays and divide by zero.
    masses_g = parse_species_masses_g(args.network_parameters)
    if len(masses_g) != len(raw_names):
        raise RuntimeError(
            f"{args.network_parameters}: {len(raw_names)} species_N_name entries "
            f"but {len(masses_g)} species_N_mass entries"
        )
    masses_amu = [m / AMU_IN_GRAMS for m in masses_g]
    charges = [species_proton_number(n) for n in raw_names]

    aion = ",\n  ".join(f"{a:.6e}" for a in masses_amu)
    aion_inv = ",\n  ".join(f"{1.0 / a:.6e}" for a in masses_amu)
    zion = ",\n  ".join(f"{float(z):.6e}" for z in charges)
    aion_constexpr = "\n  ".join(
        f"case {n}: a = {a:.6e}; break;" for n, a in zip(sanitized_names, masses_amu)
    )
    zion_constexpr = "\n  ".join(
        f"case {n}: z = {float(z):.6e}; break;" for n, z in zip(sanitized_names, charges)
    )

    bands_ev = parse_bands_ev(args.jaff_toml)
    num_chem_bands = len(bands_ev) - 1
    chem_bands = ", ".join(
        "std::numeric_limits<double>::infinity()" if edge == float("inf") else f"{edge:.6e}"
        for edge in bands_ev
    )
    power_law_index = parse_power_law_index(args.jaff_toml)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "# Auto-generated by write_jaff_network_vars.py -- do not edit\n"
        f"set(NSPEC {len(sanitized_names)})\n"
        f'set(SPECIES_ENUM\n"{species_enum}"\n)\n'
        f'set(SPEC_NAMES\n"{spec_names}"\n)\n'
        f'set(SHORT_SPEC_NAMES\n"{short_spec_names}"\n)\n'
        f'set(AION\n"{aion}"\n) # amu\n'
        f'set(AION_INV\n"{aion_inv}"\n)\n'
        f'set(ZION\n"{zion}"\n)\n'
        f'set(AION_CONSTEXPR\n"{aion_constexpr}"\n)\n'
        f'set(ZION_CONSTEXPR\n"{zion_constexpr}"\n)\n'
        f"set(NUM_CHEM_BANDS {num_chem_bands})\n"
        f'set(CHEM_BANDS "{chem_bands}") # eV\n'
        f"set(POWER_LAW_INDEX {power_law_index:.6e}) # jaff network.radiation.power_law_index\n"
    )


if __name__ == "__main__":
    main()
