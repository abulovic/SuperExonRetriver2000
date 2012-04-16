'''
Created on Apr 15, 2012

@author: intern
'''

import re

from utils.ConfigurationReader import *

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


    
if __name__ == '__main__':
    (prot_ids_known, prot_ids_abinitio) = get_gene_regions("ENSP00000311134")
    print prot_ids_known, prot_ids_abinitio
