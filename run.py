#!/usr/bin/env python

import sys
import json
import shutil
import os


sys.path.insert(0, 'src') # add library code to path
from etl import get_data
from process_data import process_data
from conversion import convert_data


DATA_PARAMS = 'config/data-params.json'
TEST_PARAMS = 'config/test-params.json'
CONVERT_PARAMS = 'config/convert-params.json'


def load_params(fp):
    with open(fp) as fh:
        param = json.load(fh)
        
    return param


def main(targets):
    
    # make the clean target
    if 'clean' in targets:
        shutil.rmtree('data/',ignore_errors=True)
        
        
    # make the conversion target
    if 'convert' in targets:
        cfg = load_params(CONVERT_PARAMS)['data']
        convert_data(**cfg)
        

    # make the data target
    if 'data' in targets:
        cfg = load_params(DATA_PARAMS)['data']
        get_data(**cfg)
        

    # make the test target
    if 'data-test' in targets:
        cfg = load_params(TEST_PARAMS)['data']
        get_data(**cfg)


    # make the process target
    if 'process' in targets:
        cfg = load_params(TEST_PARAMS)['process']
        process_data(**cfg)
        
        
    # make the test project target
    if 'test-project' in targets:
        cfg_data = load_params(TEST_PARAMS)['data']
        get_data(**cfg_data)
        
        cfg_process = load_params(TEST_PARAMS)['process']
        process_data(**cfg_process, test=True)

    return


if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)
