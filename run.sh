python create_query_list.py data/MassBank-data results/query_list.tsv && \
python request.py results/query_list.tsv && \
python modify_massbank_data.py data/MassBank-data results/mapping.json # && \
# python plot.py results/merged_results.json