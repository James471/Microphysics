#!/usr/bin/env python3

"""Derive SPECIES_ENUM / SPEC_NAMES / NSPEC / NUM_CHEM_BANDS / CHEM_BANDS for a
jaff-generated photochemistry network, so problem CMakeLists.txt files don't
have to hand-duplicate species names/order or radiation band edges that jaff
already emits.

Species names and order come from the network's own _parameters file
(species_N_name entries, written by jaff in its own species order). jaff's
raw names (e.g. "H+", "e-") are not valid C++ identifiers, so each name is
sanitized into one here: a trailing run of "+"/"-" characters is collapsed
into a count+sign suffix (single charge: "p"/"m"; multi-charge: "<n>p"/"<n>m"),
and "e-" is special-cased to "e". The same sanitized name is used for both
SPECIES_ENUM (as the enum tag) and SPEC_NAMES (as the lookup string), so
network_spec_index() and the enum agree on the species vocabulary by
construction instead of by two people typing the same list twice.

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

# 1 atomic mass unit in grams; jaff writes species_N_mass in grams, while
# aion[] is conventionally in amu.
AMU_IN_GRAMS = 1.66053906660e-24


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


def species_proton_number(raw_name):
    """Proton number Z of a species from its jaff name, e.g. "H"/"H+" -> 1, "e-" -> 0.

    Microphysics' zion[] is a proton count, not a signed charge: composition()
    forms y_e = sum(xn*zion*aion_inv), the electron fraction, which must stay
    non-negative. Using signed charge here would make y_e (and hence mu_e = 1/y_e)
    negative, since aion_inv is ~1823x larger for the electron than for anything
    else.

    jaff names are element symbols with a trailing charge run ("H", "H+", "e-").
    Only the hydrogen/electron species the photoionization networks use are
    recognised; anything else is rejected rather than guessed at, since a wrong Z
    would silently corrupt y_e.
    """
    base = re.match(r"^(.*?)([+-]*)$", raw_name).group(1)
    if base == "e":
        return 0
    if base == "H":
        return 1
    raise RuntimeError(
        f"cannot determine proton number for species {raw_name!r}: "
        "extend species_proton_number() with this element"
    )


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
    if raw_name == "e-":
        return "e"

    m = re.match(r"^(.*?)([+-]*)$", raw_name)
    base, charge_run = m.group(1), m.group(2)

    if not charge_run:
        return base

    sign = "p" if charge_run[0] == "+" else "m"
    n = len(charge_run)
    suffix = sign if n == 1 else f"{n}{sign}"
    return base + suffix


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
