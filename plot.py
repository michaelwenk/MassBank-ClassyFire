import plotly.graph_objects as go
import sys
import json


if __name__ == '__main__':
    with open(f'{sys.argv[1]}', 'r') as f:
        input_data = json.load(f)
        
        labels = list()
        parents = list()
        values = list()
        for i, accession in enumerate(input_data):
            print(f'Processing item {i + 1}/{len(input_data)} -> {accession}')
            item = input_data[accession]  
            if 'kingdom' in item and item['kingdom'] is not None:          
                kingdom_name = item['kingdom']['name']
                if kingdom_name not in labels:
                    labels.append(kingdom_name)
                    parents.append('')
                    values.append(0)
                else:
                    index = labels.index(kingdom_name)
                    values[index] += 1   

                if 'superclass' in item and item['superclass'] is not None:
                    superclass_name = item['superclass']['name']
                    if superclass_name not in labels:
                        labels.append(superclass_name)
                        parents.append(kingdom_name)
                        values.append(0)
                    else:
                        index = labels.index(superclass_name)
                        values[index] += 1

                    if 'class' in item and item['class'] is not None:
                        class_name = item['class']['name']
                        if class_name not in labels:
                            labels.append(class_name)
                            parents.append(superclass_name)
                            values.append(0)
                        else:
                            index = labels.index(class_name)
                            values[index] += 1

                        if 'subclass' in item and item['subclass'] is not None:
                            subclass_name = item['subclass']['name']
                            if subclass_name not in labels:
                                labels.append(subclass_name)
                                parents.append(class_name)
                                values.append(0)
                            else:
                                index = labels.index(subclass_name)
                                values[index] += 1                           
            
        print(labels)
        print(parents)
        print(values)
        print(len(labels))
        print(len(parents))
        print(len(values))
        
        fig = go.Figure()
        fig.add_trace(go.Sunburst(
            ids=labels,
            labels=labels,
            parents=parents,
            values=values,
            branchvalues="total",
            insidetextorientation='radial',
        ))
        fig.update_layout(
            margin = dict(t=10, l=10, r=10, b=10)
        )
        fig.show()