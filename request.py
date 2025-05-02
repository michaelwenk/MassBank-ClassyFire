import sys
import pandas as pd
import requests
import json
import time
from rdkit import Chem

from rdkit import RDLogger

# Suppress RDKit warnings globally
RDLogger.DisableLog('rdApp.*')


url = "http://classyfire.wishartlab.com"

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


if __name__ == '__main__':
    # 1. Determine Input File Type and Load InChIs
    df = pd.read_csv(sys.argv[1], sep='|', header=None)
    
    accessions = df.iloc[:, 0].tolist()
    inchis = df.iloc[:, 1].tolist()

    accessions = [accession.strip() for accession in accessions]
    inchis = [inchi.strip() for inchi in inchis]
   
    # 2. Prepare InChI List (deduplicated)
    inchi_input_list = list(set(inchis))
    print(f'Total unique input InChIs: {len(inchi_input_list)}')


    # 3. Convert InChI to SMILES
    smiles_input_list = list()
    for inchi in inchi_input_list:
        try:
            smiles = inchi_to_canonical_smiles(inchi)
            if smiles:
                smiles_input_list.append(smiles)
            else:
                print(f'Invalid InChI: {inchi}')
        except Exception as e:
            print(f'Error converting InChI to SMILES: {e}')
            continue
    # Remove duplicates
    smiles_input_list = list(set(smiles_input_list))
    print(f'Total unique input SMILES: {len(smiles_input_list)}')

    # 4. Split into chunks
    chunk_size = 50
    chunks = [smiles_input_list[i:i + chunk_size] for i in range(0, len(smiles_input_list), chunk_size)]
    print(f'Total chunks: {len(chunks)}')   

    # 5. Submit each chunk to ClassyFire
    for i, chunk in enumerate(chunks):
        print(f'-> Processing chunk {i + 1}/{len(chunks)}')
        compound = '\n'.join(chunk)
        query_id = structure_query(compound)
        print(f'Query ID: {query_id}')

        # Simulate waiting for the results to be ready
        print('Waiting for 10 seconds before fetching results...')        
        time.sleep(10)

        # Fetch results for the current chunk
        print(f'Fetching results for chunk {i + 1}/{len(chunks)}')
        result = get_results(query_id)
        with open(f'results/intermediate_results/chunk_{i + 1}_results.json', 'w') as f:
            f.write(result)
        print(f'Results for chunk {i + 1} saved to results/intermediate_results/chunk_{i + 1}_results.json')

       
    # 6. Process the results
    # Load the results from the intermediate files
    invalid_entities = list()
    results = dict()
    for i in range(len(chunks)):
        with open(f'results/intermediate_results/chunk_{i + 1}_results.json', 'r') as f:
            result = json.load(f)
            print(f'\nProcessing results for chunk {i + 1}')
            if result["invalid_entities"]:
                print(f'Invalid SMILES in chunk {i + 1}: {len(result["invalid_entities"])}')
                invalid_entities.append(result["invalid_entities"])           
            entities = result["entities"]
            if not entities:
                print(f'No entities found in chunk {i + 1}')
            else:
                print(f'Found {len(entities)} entities in chunk {i + 1}')
            
            # Process the entities
            for entity in entities:
                # Extract the canonical SMILES
                smiles = entity.get("smiles")
                try:
                    inchi_key = smiles_to_inchi_key(smiles)
                    if inchi_key:
                        results[inchi_key] = entity
                    else:
                        print(f'Invalid InChIKey for SMILES: {smiles}')                    
                    
                except Exception as e:
                    print(f'Error converting SMILES to canonical SMILES: {e}')
                    continue                
            
    print(f'\n\nInvalid entities: {len(invalid_entities)}')
    print(f'Total results: {len(results)}')
    # Save the results to a JSON file
    with open('results/merged_results.json', 'w') as f:
        json.dump(results, f, indent=4)

    # 7. Check if all SMILES are present in the results collect the data
    not_found_smiles = list()
    mapping = dict()
    for i, inchi in enumerate(inchis):
        smiles = inchi_to_canonical_smiles(inchi)
        inchi_key = smiles_to_inchi_key(smiles)
        if inchi_key not in results:
            not_found_smiles.append(smiles)
            continue
        
        accession = accessions[i]
        mapping[accession] = results[inchi_key]

        # plot the results with plotly


            
    print(f'Total SMILES found in results: {len(smiles_input_list) - len(not_found_smiles)} / {len(smiles_input_list)}')    
    print(f'Total SMILES not found in results: {len(not_found_smiles)} / {len(smiles_input_list)}')

    # 8. Save the mapping to a JSON file
    with open('results/mapping.json', 'w') as f:
        json.dump(mapping, f, indent=4)  
    
    
    
