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

from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path

VALID_CXX_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def parse_species_names(parameters_file: Path) -> list[str]:
    pattern = re.compile(r'^species_(\d+)_name\s+string\s+"(.*)"\s*$')
    entries: dict[int, str] = {}
    for line in parameters_file.read_text().splitlines():
        m = pattern.match(line.strip())
        if m:
            entries[int(m.group(1))] = m.group(2)
    if not entries:
        raise RuntimeError(f"no species_N_name entries found in {parameters_file}")
    return [entries[i] for i in sorted(entries)]


def sanitize_species_name(raw_name: str) -> str:
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


def parse_bands_ev(jaff_toml_file: Path) -> list[float]:
    with jaff_toml_file.open("rb") as f:
        config = tomllib.load(f)

    bands_ev = config["network"]["radiation"]["bands"]

    parsed: list[float] = []
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


def parse_power_law_index(jaff_toml_file: Path) -> float:
    with jaff_toml_file.open("rb") as f:
        config = tomllib.load(f)

    return float(config["network"]["radiation"]["power_law_index"])


def main() -> None:
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
        f"set(NUM_CHEM_BANDS {num_chem_bands})\n"
        f'set(CHEM_BANDS "{chem_bands}") # eV\n'
        f"set(POWER_LAW_INDEX {power_law_index:.6e}) # jaff network.radiation.power_law_index\n"
    )


if __name__ == "__main__":
    main()
