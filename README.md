# Visualise transmembrane annotation on a protein structure.

This repository contains a Streamlit app that annotates transmembrane topology on protein structures using DeepTMHMM predictions and UniProt annotations.

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tm-annotations.streamlit.app/)

## Features

- Retrieves transmembrane topology annotations from a UniProt entry.
- Map UniProt annotation onto a protein structure using updated mmCIF files from PDBe..
- Run a DeepTMHMM prediction and map the predicted topology onto the structure.
- Independent Python script to run as a PyMOL plugin.
  
### Input
- **UniProt accession**
- **PDB code**

## Prerequisites

- **Python 3.x**
- Python libraries: `streamlit`, `py3Dmol`, `pybiolib`, `requests`, `matplotlib`.

## PyMOL plugin

In PyMOL 3.0, Plugin -> Plugin Manager -> Install from PyMOLWiki or any URL:
https://github.com/connyyu/TM_annotations/blob/main/scripts/colortm.py
Usage: colortm <object_name>, <uniprot_accession>
    
## Author

- **Conny Yu** – [GitHub Profile](https://github.com/connyyu)  
  Haku series 1.2 _June 2025_
