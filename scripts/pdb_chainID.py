import csv
import re
import os
import sys
import glob
import requests
from collections import defaultdict

script_dir = os.path.dirname(os.path.abspath(__file__))  # location of pages directory

def download_cif_file(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_code_dl}.cif", 'wb') as f:
            f.write(response.content)
        return f"{pdb_code_dl}.cif"
    else:
        print(f"\033[1mFailed to download CIF file for {pdb_code_dl}\033[0m")
        return None

def parse_resolution(cif_file):
    resolution = None
    method = ""
    symmetry_type_found = False

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith('_reflns.d_resolution_high'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "Xtal"
                break
        elif line.startswith('_em_3d_reconstruction.resolution'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "EM"
                break
        elif line.startswith("_exptl.method 'Solution NMR'"):
            method = "NMR"
            break
        elif line.startswith('_em_3d_reconstruction.symmetry_type'):
            symmetry_type_found = True
        elif symmetry_type_found:
            parts = line.split()
            if len(parts) >= 7:
                resolution = parts[6]
                method = "EM"
            break
    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None

    return resolution, method

def parse_hetatms(cif_file):
    hetatms = set()
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('HETATM'):
            parts = line.split()
            hetatm_name = parts[5]
            hetatms.add(hetatm_name)
    return sorted(hetatms)

def fetch_uniprot_id(uniprot_ac):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get('uniProtkbId', None)
    return None

def parse_uniprot_mapping(cif_file, chains_to_check=None):
    uniprot_mapping = {}
    uniprot_ids = []
    chain_mapping = {}

    fallback_accession = None
    fallback_ID = None
    fallback_chain = None

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            fallback_accession = (line.split()[1] if len(line.split()) > 1
                                  else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref.db_code'):
            fallback_ID = (line.split()[1] if len(line.split()) > 1
                           else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref_seq.pdbx_strand_id'):
            parts = line.split()
            if len(parts) > 1:
                fallback_chain_candidate = parts[1]
                if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                    fallback_chain = fallback_chain_candidate
            else:
                k = i + 1
                while k < len(lines):
                    next_line = lines[k].strip()
                    if not next_line or next_line.startswith('_'):
                        break
                    fallback_chain_candidate = next_line.split()[0]
                    if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                        fallback_chain = fallback_chain_candidate
                        break
                    k += 1

    i = 0
    process_sifts = False
    sifts_mapping = {}
    sifts_chain_mapping = {}

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('_struct_ref.pdbx_db_accession'):
            j = i + 1
            while j < len(lines) and not lines[j].startswith('loop_') and lines[j].strip():
                line_content = lines[j].strip()
                
                # Skip lines that are part of multi-line sequences (start with ; or are continuation)
                if line_content.startswith(';') or line_content.endswith(';'):
                    if line_content.startswith(';') and not line_content.endswith(';'):
                        j += 1
                        while j < len(lines):
                            if lines[j].strip().endswith(';'):
                                break
                            j += 1
                    j += 1
                    continue
                parts = line_content.split()
                if 'UNP' in parts and len(parts) >= 4:  # Need at least entity_id, UNP, db_code, accession
                    entity_id = parts[0]
                    uniprot_id = parts[2]
                    uniprot_accession = parts[3]
                    # Validate that the accession looks like a UniProt accession
                    if len(uniprot_accession) >= 6 and '-' not in uniprot_accession and '_' not in uniprot_accession:
                        uniprot_mapping[entity_id] = uniprot_accession
                        uniprot_ids.append((entity_id, uniprot_id))
                j += 1
            i = j
            continue

        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1
            while i < len(lines):
                next_line = lines[i].strip()
                if not next_line or next_line.startswith(('loop_', '#')):
                    break
                parts = next_line.split()
                if len(parts) > 3:
                    entity_id, chain_id = parts[1], parts[3]
                    chain_mapping.setdefault(entity_id, set()).add(chain_id)
                i += 1
            continue

        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True

        elif process_sifts:
            if line.startswith(('loop_', '#')):
                process_sifts = False
            elif line:
                parts = line.split()
                if len(parts) > 3:
                    entity_id, chain_id, sifts_accession = parts[0], parts[1], parts[2]
                    if len(sifts_accession) >= 6 and '-' not in sifts_accession:
                        sifts_mapping[entity_id] = sifts_accession
                        sifts_chain_mapping.setdefault(entity_id, set()).add(chain_id)
        i += 1

    if not uniprot_mapping:
        if sifts_mapping:
            fallback_entity = sorted(sifts_mapping.keys())[0]
        elif chain_mapping:
            fallback_entity = sorted(chain_mapping.keys())[0]
        else:
            fallback_entity = "1"

        uniprot_mapping[fallback_entity] = fallback_accession
        uniprot_ids.append((fallback_entity, fallback_ID))

        chains_for_entity = chain_mapping.get(fallback_entity)
        if chains_for_entity:
            chain_mapping[fallback_entity] = chains_for_entity
        else:
            chain_mapping[fallback_entity] = set([fallback_chain] if fallback_chain else [])

    final_uniprot_mapping = {}
    final_uniprot_ids = []
    final_chain_mapping = {}

    for entity in chain_mapping:
        original_ac = uniprot_mapping.get(entity, None)
        uniprot_id = next((uid for eid, uid in uniprot_ids if eid == entity), None)
        sifts_ac = sifts_mapping.get(entity, None)

        # Use SIFTS AC if available; else original AC
        final_ac = sifts_ac if sifts_ac else original_ac

        # Use only original chains — NO merging with SIFTS chains
        chains_orig = chain_mapping[entity]
        final_chains = sorted(list(chains_orig))

        final_uniprot_mapping[entity] = final_ac
        final_uniprot_ids.append((entity, uniprot_id))
        final_chain_mapping[entity] = final_chains

    return final_uniprot_mapping, final_uniprot_ids, final_chain_mapping

def check_chimera(table_data):
    chain_uniprot_map = defaultdict(set)
    for row in table_data:
        pdb_code = row['PDB'].upper()
        chain_id = row['Chain ID']
        uniprot_accession = row['UniProt AC']
        chain_uniprot_map[chain_id].add(uniprot_accession)
    warnings = []
    for chain_id, uniprot_accessions in chain_uniprot_map.items():
        if len(uniprot_accessions) > 1:
            warnings.append(f"{pdb_code} chain {chain_id} is associated with mutiple UniProt AC: {', '.join(uniprot_accessions)}.")
    return warnings

def process_cif_file(pdb_code):
    cif_file = download_cif_file(pdb_code)
    if not cif_file:
        print(f"\033[1mFailed to download CIF file for {pdb_code}. Skipping...\033[0m")
        return None

    resolution, method = parse_resolution(cif_file)
    hetatms = parse_hetatms(cif_file)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_file)

    hetatm_names = ', '.join(hetatms)
    table_data = []

    printed_chains = set()
    for entity_id, uniprot_accession in uniprot_mapping.items():
        uniprot_id = next((id for ent_id, id in uniprot_ids if ent_id == entity_id), None)
        chain_ids = chain_mapping.get(entity_id, "A")
        for chain_id in chain_ids:
            if (uniprot_accession, chain_id) in printed_chains:
                continue
            printed_chains.add((uniprot_accession, chain_id))
            table_data.append({
                'PDB': pdb_code.lower(),
                'Entity': entity_id,
                'Chain ID': chain_id,
                'UniProt AC': uniprot_accession,
                'UniProt ID': uniprot_id,
                'HETATM': hetatm_names,
                'Resolution': resolution,
                'Method': method
            })

    warnings = check_chimera(table_data)

    return table_data, warnings

def get_chain_ids(pdb_code, uniprot_ac=None):
    """
    Main function to retrieve chain IDs and UniProt info for a PDB code.
    Optionally filter results for a specific UniProt accession.
    Returns list of dictionaries of chain data and any warnings.
    """
    table_data, warnings = process_cif_file(pdb_code)
    if not table_data:
        return [], []

    if uniprot_ac:
        filtered = [row for row in table_data if row['UniProt AC'] == uniprot_ac]
    else:
        filtered = table_data

    return filtered, warnings

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <PDB_CODE> [UNIPROT_AC]")
        sys.exit(1)
    pdb_code = sys.argv[1]
    uniprot_ac = sys.argv[2] if len(sys.argv) > 2 else None
    results, warnings = get_chain_ids(pdb_code, uniprot_ac)

    for row in results:
        print(f"PDB: {row['PDB']}\tEntity: {row['Entity']}\tChain ID: {row['Chain ID']}\t"
              f"UniProt AC: {row['UniProt AC']}\tUniProt ID: {row['UniProt ID']}\t"
              f"HETATM: {row['HETATM']}\tResolution: {row['Resolution']}\tMethod: {row['Method']}")

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(w)

if __name__ == "__main__":
    main()
