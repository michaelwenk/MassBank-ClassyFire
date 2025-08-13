import sys
import json
import shutil
from utils import list_files_recursive

if __name__ == '__main__':
    print('\n\n--> Starting modify_massbank_data.py ...\n\n')

    if len(sys.argv) < 3:
        print("Usage: python modify_massbank_data.py <src_dir> <copy_data[true/false]>\n\n")
        sys.exit(1)

    src_dir = sys.argv[1] 
    copy_data_str = sys.argv[2]

    if copy_data_str.lower() == 'true':
        out_dir = "results/" + src_dir        
        print(f'Copying files from {src_dir} to {out_dir} to not modify original data...')
        shutil.copytree(src=src_dir, dst=out_dir, dirs_exist_ok=True)
        print(f'Copied files!')        
    elif copy_data_str.lower() == 'false':
        out_dir = src_dir
        print(f'Not copying files, using original source directory: {out_dir}')
    else:
        print("Invalid value for copy_data. Please use 'true' or 'false'.")
        sys.exit(1)
        
    files = list_files_recursive(out_dir)   
    mappings_json = "results/mapping.json"

    print(f'Total MassBank record files: {len(files)} in directory: {out_dir}')
    print(f'Using mappings from: {mappings_json}')        

    with open(f'{mappings_json}', 'r') as f:
        input_data = json.load(f)
         
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

                    if chemont_id != "" and kingdom_name != "" and superclass_name != "":
                        mod_line = "CH$LINK: ChemOnt " + chemont_id + "; " + kingdom_name + "; " + superclass_name
                        if class_name != "":
                            mod_line += "; " + class_name
                            if subclass_name != "":
                                mod_line += "; " + subclass_name
                        mod_line += "\n"

                        # check if the ChemOnt line is present
                        if chemont_line_index != -1:                        
                            lines[chemont_line_index] = mod_line
                        else:
                            # check if the last CH$LINK line is present
                            if last_chlink_line_index != -1:                           
                                lines.insert(last_chlink_line_index + 1, mod_line)
                            elif iupac_line_index != -1:                               
                                lines.insert(iupac_line_index + 1, mod_line)
                    # else:
                    #     if chemont_line_index != -1:
                    #         print(f'Accession {accession} -> no or insufficient classification data -> remove ChemOnt link')
                    #         lines.pop(chemont_line_index)                    
                # else:
                #     if chemont_line_index != -1:
                #         print(f'Accession {accession} not found in input data -> remove ChemOnt link')
                #         lines.pop(chemont_line_index)

                # write the lines back to the file
                with open(file, 'w') as f:
                    f.writelines(lines)
                # print(f'File {file} updated successfully')     


    print(f'\n\n-> Finished updating {len(files)} files in {out_dir}\n\n')
                                
           


