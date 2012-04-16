'''
Created on Apr 15, 2012

@author: intern
'''

import re, sys

from utilities.ConfigurationReader import ConfigurationReader

proteins_known = {}
proteins_abinitio = {}

def parse_descr_file (protein_id):
    
    cr = ConfigurationReader.Instance()
    descr_file_path = cr.get_value('root', 'session_dir') + "/" + protein_id + "/" + protein_id + ".descr"
    
    pattern_known = re.compile("(.*)\t(ENS.*)\t(.*)\t(.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
    pattern_abinitio = re.compile("(.*)\t(GEN.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
    
    
    descr_file = open(descr_file_path, 'r')
    
    for line in descr_file.readlines():
        line = line.strip()
        
        match = re.match(pattern_known, line)
        if match:
            (species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
            proteins_known[species_name] = (spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
            
        match = re.match(pattern_abinitio, line)    
        if match:
            (species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
            proteins_abinitio[species_name] = (spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    
    descr_file.close()
    
    return proteins_known, proteins_abinitio
    
def get_protein_ids (protein_id):
    
    (proteins_known, proteins_abinitio) = parse_descr_file(protein_id)
        
    prot_ids_known = {}
    prot_ids_abinitio = {}
    
    for key, value in proteins_known.items():
        prot_ids_known[key] = list(value)[0]
    for key, value in proteins_abinitio.items():
        prot_ids_abinitio[key] = list(value)[0]
        
    return prot_ids_known, prot_ids_abinitio

def get_gene_regions (protein_id):
    (proteins_known, proteins_abinitio) = parse_descr_file(protein_id)
    
    genes_known = {}
    genes_abinitio = {}
    
    for key, value in proteins_known.items():
        genes_known[key] = list(value)[3:]
    for key, value in proteins_abinitio.items():
        genes_abinitio[key] = list(value)[1:]
        
    return genes_known, genes_abinitio

def get_project_root_dir ():
    ex_path = sys.path[0]
    m = re.match("(.*ExoLocator).*", ex_path)
    proj_root_dir = m.groups()[0]
    return proj_root_dir

def get_species_list ():
    species_file_path = get_project_root_dir() + "/species.txt"
    species_file = open(species_file_path, 'r')
    species_list = []
    
    for line in species_file.readlines():
        species_list.append(line.strip()) 
        
    species_file.close()
    return species_list

def get_protein_list ():
    cr = ConfigurationReader.Instance()
    protein_file_path = cr.get_value('input', 'protein_list')
    protein_file = open(protein_file_path, 'r')
    protein_list = []
    
    for line in protein_file.readlines():
        (prot_id, num_of_exons) = line.strip().split()
        protein_list.append((prot_id, num_of_exons))
    protein_file.close()
    
    return protein_list
    


    
if __name__ == '__main__':
    (prot_ids_known, prot_ids_abinitio) = get_gene_regions("ENSP00000311134")
    print prot_ids_known, prot_ids_abinitio
