import sys
import pandas as pd
from utils import list_files_recursive

if __name__ == '__main__':
    directory_path = sys.argv[1] 
    files = list_files_recursive(directory_path)
    output_file = sys.argv[2]
    
    print(f'Processing files in directory: {directory_path}')
    print(f'Total files: {len(files)}')

    accessions: list[str] = []
    inchis: list[str] = []
    for i, file in enumerate(files):
        lines = list() 
        # print(f'Processing file {i + 1}/{len(files)} -> {file}')
        with open(file, 'r') as f:
            lines = f.readlines()
            chemont_line_index = -1
            last_chlink_line_index = -1
            iupac_line_index = -1
            accession = ""
            for j, line in enumerate(lines):
                if line.startswith('ACCESSION:'):
                    accession = line.split(':')[1].strip()
                elif line.startswith('CH$LINK: ChemOnt'):
                    chemont_line_index = j
                elif line.startswith('CH$LINK:'):
                    last_chlink_line_index = j
                elif line.startswith('CH$IUPAC:'):
                    iupac_line_index = j
                elif line.startswith('CH$SMILES:'):
                    smiles_line_index = j
            if(accession != "" and chemont_line_index == -1 and iupac_line_index != -1):           
                inchi = lines[iupac_line_index].split(":")[1].strip()                
                if inchi.startswith('InChI='):                    
                    accessions.append(accession)
                    inchis.append(inchi)

    print(f'Total accessions without ChemOnt classification: {len(accessions)}')
    print(f'Total InChIs without ChemOnt classification: {len(inchis)}')

    # create dataframe
    df = pd.DataFrame({
        'accession': accessions,
        'inchi': inchis
    })    

    # write the dataframe to a TSV file    
    df.to_csv(output_file, index=False, header=False, sep='|')
    print(f'Query list saved to {output_file}')   