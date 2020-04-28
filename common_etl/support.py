from datetime import datetime
from Bio.PDB import *
import sys
import re 
import math 
from itertools import islice 

def first_comparison(chain_1, chain_2):

    results = []
    for atom_1 in chain_1.get_atoms():
        for atom_2 in chain_2.get_atoms():
            if atom_1 - atom_2 < 4.5:
                results.append((atom_1.get_parent(), atom_2.get_parent()))

    return results 


def write_to_file(my_list,file_output_name):
    with open(file_output_name, 'w') as your_file:
        for item in my_list:
            your_file.write(f"{item}\n")



def selfcompare(chain1, epitope_difference_list) :
    results = []
    # count = 0
    for i in chain1.get_atoms() :
        for x in chain1.get_atoms() :
           print(i,x)

def second_comparison(chain_1, chain_2):

    results = []
    count = 0

    for atom_1 in chain_1.get_atoms() :
        for atom_2 in chain_2.get_atoms() :
            if (
                    atom_1.get_parent().get_resname() == "SER" or atom_1.get_parent().get_resname() == "THR" or atom_1.get_parent().get_resname() == "ALA" or
                atom_1.get_parent().get_resname() == "CYS") and (
                    atom_1.get_fullname() == " CB ") :
                if (
                        atom_2.get_parent().get_resname() == "SER" or atom_2.get_parent().get_resname() == "THR" or atom_2.get_parent().get_resname() == "ALA"or
                    atom_2.get_parent().get_resname() == "CYS") and (
                        atom_2.get_fullname() == " CB ") :
                    count += 1
                    if (atom_1 - atom_2 <= 7.34):  # 7.34 :  # 5.81 + 0.67*s
                        results.append((atom_1.get_parent(), atom_2.get_parent(), atom_1 - atom_2))
    return results


def third_comparison(chain_1):
    
    results = []
    for atom_1 in chain_1.get_atoms() :
        if atom_1.get_parent().get_resname() == "SER" :
            results.append(atom_1.get_parent())
    
    return results

def distance_2(aa_residue_pos, epitope_positions):
    
    aa_1 = (aa_residue_pos[0][0])
    aa_2 = (aa_residue_pos[1][0])
    residue_pos_1 = int(aa_residue_pos[0][1])
    residue_pos_2 = int(aa_residue_pos[1][1])
    
    epitope_difference_list = []
    average_distance = None
    residue_difference = abs(residue_pos_1 - residue_pos_2)
    if residue_difference == 0:
        pass
    else:
        for epitope_position in epitope_positions:
            distance_1 = abs(residue_difference - epitope_position)
            epitope_difference_list.append(distance_1)
    
    
    if len(epitope_difference_list) == 0:
        pass
    else:
        epitope_difference_sum = sum(epitope_difference_list)
        average_distance = round(epitope_difference_sum / len(epitope_difference_list),3)
                
    return average_distance


    # num_residues = 0
    # chainA_to_residue_site_distances_list = []
    # chainB_to_residue_site_distances_list = []
    # for residue_position in residue_positions:
    #     num_residues += 1
    #     chainA_to_residue_site_distance = abs(residue_pos_1 - residue_pos_2)
    #     chainA_to_residue_site_distances_list.append(chainA_to_residue_site_distance)
    #     chainB_to_residue_site_distance = abs(residue_pos_2 - int(residue_position))
    #     chainB_to_residue_site_distances_list.append(chainB_to_residue_site_distance)
    #     distance_1 = chainA_to_residue_site_distance + chainB_to_residue_site_distance
    
    # distance_output = distance_1 / (num_residues * 2)

    # return chainA_to_residue_site_distances_list, chainB_to_residue_site_distances_list, distance_output, aa_1, aa_2


def distance_1(chain,epitope_list,residue_list):
    '''
        Distance Calculation from specific amino acid to the 
        subsitution site 

        return int
    '''
    results = []
    distance_list = []
    avg_distance_list = []
    for i in chain.get_atoms() :
        for x in chain.get_atoms() :
            if( 
            x.get_parent().id[1] in epitope_list and x.get_fullname() == " CB " 
            and i.get_fullname() == " CB " and i.get_parent().id[1] == residue_list):
                distance_list.append(i-x)
        
    
    
    
    avg_distance = sum(distance_list) / len(epitope_list)
    
    return round(avg_distance,2)
    


def filter_results(line):
    
    aa_res_tuple = ()
    aa_pattern = r"[A-Z]{3}"
    num_pattern = r"\d+\S"
    aa = re.findall(aa_pattern, line)
    res_pos = re.findall(num_pattern, line)
    aa_res_tuple= tuple(zip(aa,res_pos))
   
    return aa_res_tuple