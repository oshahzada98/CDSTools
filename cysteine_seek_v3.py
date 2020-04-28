from common_etl.support import *
from datetime import datetime
from Bio.PDB import *
import sys
import yaml
import io
import pandas as pd 
import csv

'''
----------------------------------------------------------------------------------------------
The configuration reader. Parses the YAML configuration into dictionaries
'''
def load_config(yaml_config):
    yaml_dict = None
    config_stream = io.StringIO(yaml_config)
    try:
        yaml_dict = yaml.load(config_stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as ex:
        print(ex)

    if yaml_dict is None:
        return None, None, None

    return  yaml_dict['file_name_filters'],yaml_dict['epitope_positions'],yaml_dict['residue_positions'], yaml_dict['steps']



def main():   
    
    with open(sys.argv[1], mode='r') as yaml_file:
        params, epitope_positions,residue_positions, steps = load_config(yaml_file.read())

    parser = PDBParser() # Parses the PDB file 
    structure = parser.get_structure("S", "6vsb.pdb")
    
    chains = [chain_id for chain_id in structure.get_chains()]
    chain_a, chain_b, chain_c = chains[0], chains[1], chains[2]
    
    if 'first_comparison' in steps:
        print('first_comparison running!')
        results_ab = set(first_comparison(chain_a,chain_b))
        if 'write_to_file' in steps:
            print('write_to_file running!')
            write_to_file(results_ab, params['file_output_name_1'])
    
    if 'second_comparison' in steps:
        print('second_comparison running!')
        results_ac = second_comparison(chain_a,chain_c,epitope_positions)
        if 'write_to_file' in steps:
            write_to_file(results_ac, params['file_output_name_2'])
    
    if 'third_comparison' in steps: 
        print('third_comparison running!')
        results_ab = third_comparison(chain_a,chain_b)
        if 'write_to_file' in steps:
            write_to_file(results_ab, params['file_output_name_3'])

    if 'distance_1' in steps:
        print('distance_1')
        avg_distance_list = []
        for residue_position in residue_positions:
            avg_distance = distance_1(chain_a, epitope_positions,residue_position)
            print(avg_distance)
            avg_distance_list.append(avg_distance)
        
        print(f'min: {min(avg_distance_list)}, max: {max(avg_distance_list)}')


if __name__ == "__main__":
    main()