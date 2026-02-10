from pymol import cmd, util
import requests
import re


def fetch_uniprot_tm_helices(uniprot_ac):
    """
    Fetch transmembrane annotations from UniProt.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.txt"
    response = requests.get(url)

    if response.status_code != 200:
        raise RuntimeError(f"Failed to fetch UniProt entry: {uniprot_ac}")

    tm_helices = []
    lines = response.text.splitlines()

    for line in lines:
        match = re.match(r"FT\s+TRANSMEM\s+(\d+)\.\.(\d+)", line)
        if match:
            start, end = map(int, match.groups())
            tm_helices.append((start, end))

    return tm_helices

def generate_pymol_coloring(obj_name, residue_ranges):
    """
    Apply coloring to TM helices usinng UniProt annotations.
    """
    util.cbaw(obj_name)
    color_cmds = []
    for i, (start, end) in enumerate(residue_ranges, 1):
        # Remap similar colors
        color_index = 31 if i == 1 else 30 if i == 7 else i
        color_cmds.append(f"color {color_index}, \"{obj_name}\" and resi {start}-{end}")
    full_command = "; ".join(color_cmds)
    cmd.do(full_command)

def colortm(*args, _self=None):
    if len(args) == 0:
        print("Usage: colortm <object_name>, <uniprot_accession>")
        return
    if len(args) == 1 and "," in args[0]:
        parts = [p.strip() for p in args[0].split(",")]
    else:
        parts = [str(a).strip() for a in args]
    if len(parts) != 2:
        print("Usage: colortm <object_name>, <uniprot_accession>")
        return

    obj_name, uniprot_ac = parts
    uniprot_ac = uniprot_ac.upper()  # Normalize UniProt ID

    # Fetch TM helices
    try:
        tm_helices = fetch_uniprot_tm_helices(uniprot_ac)
    except Exception as e:
        print(f"ERROR: {e}")
        return
    if not tm_helices:
        print(f"ERROR: UniProt entry {uniprot_ac} has no TRANSMEM annotations.")
        return
    # Apply coloring
    generate_pymol_coloring(obj_name, tm_helices)
    # Print TM annotations
    for start, end in tm_helices:
        print(f"  FT TRANSMEM  {start:>3}  {end}")
    print(f"{obj_name} with TRANSMEM annotations for UniProt entry {uniprot_ac}.")

cmd.extend("colortm", colortm)
