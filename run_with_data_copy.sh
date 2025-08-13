#!/bin/bash

python create_query_list.py $1 && \
python request.py $2 && \
python modify_massbank_data.py $1 true