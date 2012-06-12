'''
Created on Apr 27, 2012

@author: marioot
'''
import csv
from utilities.Logger import Logger
import os
import re
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.GeneContainer import GeneContainer
from data_analysis.utilities import generate_structure

def create_suspect_key_dict(csv_path):
    protein_keys = {}
    gene_keys    = {}
    
    with open(csv_path, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[4] == ' Protein':
                species_id = row[6].replace(' ', '')
                ref_prot_id = row[2].replace(' ', '')
                if protein_keys.has_key(species_id):
                    protein_keys[species_id].append(ref_prot_id) 
                else:
                    protein_keys[species_id] = [ref_prot_id]
            else:
                species_id = row[6].replace(' ', '')
                ref_prot_id = row[2].replace(' ', '')
                if gene_keys.has_key(species_id):
                    gene_keys[species_id].append(ref_prot_id) 
                else:
                    gene_keys[species_id] = [ref_prot_id]
    return (protein_keys, gene_keys)
            
def get_csv_path():
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('containers')
    
    #path do logfilea: alignment_logger.handlers[0].baseFilename
    
    containers_log_path = alignment_logger.handlers[0].baseFilename

    if os.path.exists(containers_log_path):
        containers_csv_path = re.sub('.log', '.csv', containers_log_path)
    else:
        containers_csv_path = containers_log_path
    os.rename(containers_log_path, containers_csv_path)
    return containers_csv_path


def main ():
    generate_structure.main()
    
    pc = ProteinContainer.Instance()
    gc = GeneContainer.Instance()
    
    (protein_keys, gene_keys) = create_suspect_key_dict(get_csv_path())

    for key in protein_keys:
        protein_keys[key].append(pc.get(key).ref_protein)
        
    for key in gene_keys:
        gene_keys[key].append(gc.get(key).ref_protein)
   
    suspect_proteins = []
    for key in protein_keys:
        if suspect_proteins.count(protein_keys[key]) == 0:
            suspect_proteins.append(protein_keys[key])
    for key in gene_keys:
        if suspect_proteins.count(gene_keys[key]) == 0:
            suspect_proteins.append(gene_keys[key])
    
    suspect_proteins_sets = []
    for protein_list in suspect_proteins:
        suspect_proteins_sets.append(set(protein_list))
    #print suspect_proteins_sets
    
    final_sets = []
    intersection_hapened = True
    while (intersection_hapened):
        intersection_hapened = False
        #print suspect_proteins_sets
        #print final_sets
        #print "\n"
        for set1 in suspect_proteins_sets:
            suspect_proteins_sets.remove(set1)
            #print "Removing %s"%set1
            for set2 in suspect_proteins_sets:
                if set1.intersection(set2) != set():
                    suspect_proteins_sets.remove(set2)
                    intersection_hapened = True
                    set1 = set1.union(set2)
                    #print "Found intersection: %s"%set2
                    #print "New union: %s"%set1
                    break
            if intersection_hapened:
                suspect_proteins_sets.append(set1)
                break
            else:
                #print "Done: %s"%set1
                final_sets.append(set1)
    #print final_sets           

                
    
    output_file = open('/home/intern/Desktop/wierd_protein_groups1.txt', 'w')
    for protein_list_set in final_sets:
        protein_list = list(protein_list_set)
        output_file.write("{0}\n".format(', '.join(protein_list)))
    output_file.close()
    
    pass
    
if __name__ == '__main__':
    main()