#!/bin/bash

python -u create_query_list.py $1 && \
python -u request.py $2 && \
python -u modify_massbank_data.py $1 true