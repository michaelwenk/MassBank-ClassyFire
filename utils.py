import os
import json
import requests
from rdkit import Chem
from rdkit import RDLogger

# Suppress RDKit warnings globally
RDLogger.DisableLog('rdApp.*')

url = "http://classyfire.wishartlab.com"

def smiles_to_canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    return canonical_smiles


def inchi_to_canonical_smiles(inchi):
    mol = Chem.MolFromInchi(inchi)
    smiles = Chem.MolToSmiles(mol, canonical=True)

    return smiles


def smiles_to_inchi_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    inchi_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))

    return inchi_key


def structure_query(compound, label='IPB-Halle-ClassyFire'):   
    payload = {
        "label": label,
        "query_input": compound,
        "query_type": "STRUCTURE"
    }
    headers = {"Content-Type": "application/json"}
    r = requests.post(f"{url}/queries.json", data=json.dumps(payload), headers=headers)
    r.raise_for_status()
    return r.json()['id']


def get_results(query_id, return_format="json"):  
    r = requests.get(f"{url}/queries/{query_id}.{return_format}", headers={"Content-Type": f"application/{return_format}"})
    r.raise_for_status()
    return r.text


def list_files_recursive(path):
    files = list()
    for entry in os.listdir(path):
        full_path = os.path.join(path, entry)
        if os.path.isdir(full_path):
            files.extend(list_files_recursive(full_path))
        elif 'MSBNK-' in full_path and full_path.endswith('.txt'):
                files.append(full_path)
    return files