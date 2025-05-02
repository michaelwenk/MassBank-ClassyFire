import os
import sys
import json

def list_files_recursive(path):
    files = list()
    for entry in os.listdir(path):
        full_path = os.path.join(path, entry)
        if os.path.isdir(full_path):
            files.extend(list_files_recursive(full_path))
        elif 'MSBNK-' in full_path and full_path.endswith('.txt'):
                files.append(full_path)
    return files



if __name__ == '__main__':
    directory_path = sys.argv[1] 
    files = list_files_recursive(directory_path)
    
    print(f'Total files: {len(files)}')

    with open(f'{sys.argv[2]}', 'r') as f:
        input_data = json.load(f)
         
        for i, file in enumerate(files):
            lines = list() 
            print(f'Processing file {i + 1}/{len(files)} -> {file}')
            with open(file, 'r') as f:
                lines = f.readlines()
                chemOntLineIndex = -1
                lastChLinkLineIndex = -1
                iupacLineIndex = -1
                accession = ""
                for j, line in enumerate(lines):
                    if line.startswith('ACCESSION:'):
                        accession = line.split(':')[1].strip()
                    elif line.startswith('CH$LINK: ChemOnt'):
                        chemOntLineIndex = j
                    elif line.startswith('CH$LINK:'):
                        lastChLinkLineIndex = j
                    elif line.startswith('CH$IUPAC:'):
                        iupacLineIndex = j
                # print(f'accession: {accession}')
                # print(f'chemOntLineIndex: {chemOntLineIndex}')
                # print(f'lastChLinkLineIndex: {lastChLinkLineIndex}')
                # print(f'iupacLineIndex: {iupacLineIndex}')

                if accession in input_data:
                    item = input_data[accession]
                    kingdom_name = ""
                    superclass_name = ""
                    class_name = ""
                    subclass_name = ""
                    chemont_id = ""
                    if 'kingdom' in item and item['kingdom'] is not None:          
                        kingdom_name = item['kingdom']['name']
                        chemont_id = item['kingdom']['chemont_id']             
                        if 'superclass' in item and item['superclass'] is not None:
                            superclass_name = item['superclass']['name']
                            chemont_id = item['superclass']['chemont_id']
                            if 'class' in item and item['class'] is not None:
                                class_name = item['class']['name']        
                                chemont_id = item['class']['chemont_id']            
                                if 'subclass' in item and item['subclass'] is not None:
                                    subclass_name = item['subclass']['name']
                                    chemont_id = item['subclass']['chemont_id']
                    # print(f'kingdom_name: {kingdom_name}')
                    # print(f'superclass_name: {superclass_name}')
                    # print(f'class_name: {class_name}')
                    # print(f'subclass_name: {subclass_name}')

                    mod_line = "CH$LINK: ChemOnt " + chemont_id + "; " + kingdom_name + "; " + superclass_name
                    if class_name != "":
                        mod_line += "; " + class_name
                        if subclass_name != "":
                            mod_line += "; " + subclass_name
                    mod_line += "\n"

                    # check if the ChemOnt line is present
                    if chemOntLineIndex != -1:                        
                        lines[chemOntLineIndex] = mod_line
                    else:
                        # check if the last CH$LINK line is present
                        if lastChLinkLineIndex != -1:                           
                            lines.insert(lastChLinkLineIndex + 1, mod_line)
                        elif iupacLineIndex != -1:                               
                            lines.insert(iupacLineIndex + 1, mod_line)
                else:
                    if chemOntLineIndex != -1:
                        print(f'Accession {accession} not found in input data -> remove ChemOnt link')
                        lines.pop(chemOntLineIndex)

                # write the lines back to the file
                with open(file, 'w') as f:
                    f.writelines(lines)
                print(f'File {file} updated successfully')     
                                
           


