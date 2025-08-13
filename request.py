import os
import sys
import json
import time
import pandas as pd
from utils import get_results, inchi_to_canonical_smiles, smiles_to_inchi_key, structure_query

if __name__ == '__main__':
    print('\n\n--> Starting request.py ...\n\n')

    if len(sys.argv) < 2:
        print("Usage: python request.py <waiting_time[seconds]>")
        sys.exit(1)

    waiting_time = int(sys.argv[1])

    # 1. Read TSV file
    df = pd.read_csv("results/query_list.tsv", sep='|', header=None)
    
    accessions: list[str] = df.iloc[:, 0].tolist()
    inchis: list[str] = df.iloc[:, 1].tolist()

    print(f'Total accessions: {len(accessions)}')
    print(f'Total InChIs: {len(inchis)}')

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

    os.makedirs("results/intermediate_results", exist_ok=True)

    # 5. Submit each chunk to ClassyFire
    for i, chunk in enumerate(chunks):
        print(f'\n\n-> Processing chunk {i + 1}/{len(chunks)}')
        compound = '\n'.join(chunk)
        query_id = structure_query(compound)
        print(f'Query ID: {query_id}')

        # Simulate waiting for the results to be ready
        print(f'Waiting for {waiting_time} seconds before fetching results...')
        time.sleep(waiting_time)

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

    not_found_smiles = list(set(not_found_smiles))  # Remove duplicates

            
    print(f'Total SMILES found in results: {len(smiles_input_list) - len(not_found_smiles)} / {len(smiles_input_list)}')    
    print(f'Total SMILES not found in results: {len(not_found_smiles)} / {len(smiles_input_list)}')

    # 8. Save the mapping to a JSON file
    with open('results/mapping.json', 'w') as f:
        json.dump(mapping, f, indent=4)  
    
    
    
