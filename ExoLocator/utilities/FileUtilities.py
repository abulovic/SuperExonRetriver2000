'''
Created on Apr 15, 2012

@author: intern
'''

import re, sys

from utilities.ConfigurationReader import ConfigurationReader
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
import os



def parse_descr_file (protein_id):
    
    '''
    Function for parsing the description file associated with the protein_id.
    Description file contains two different types of entries: for known protein and for abinitio.
    Consequently, there are two formats that can be expected. They are both tab delimited.
    known_format:     species protein_id gene_id transcript_id location_type:assembly:location_id:seq_begin:seq_end:strand
    abinitio_format:  species protein_id location_type:assembly:location_id:seq_begin:seq_end:strand
    @param protein_id: protein for which the description file will be parsed. 
    @return: proteins_known_data - dictionary (key is species, value is a tuple of all the available data for that species protein)
             abinitio_known_data (the same), two dictionaries are returned as a tuple
    @raise IOError: in case there is no description file present for the protein_id
    '''
    
    proteins_known_data = {}
    proteins_abinitio_data = {}
    
    cr = ConfigurationReader.Instance()
    descr_file_path = cr.get_value('root', 'session_dir') + "/" + protein_id + "/" + protein_id + ".descr"
    
    pattern_known = re.compile("(.*)\t(ENS.*)\t(.*)\t(.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
    pattern_abinitio = re.compile("(.*)\t(GEN.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
    
    try:
        descr_file = open(descr_file_path, 'r')
    except IOError:
        raise IOError("There is no description file present for protein: %s" % protein_id)
    
    for line in descr_file.readlines():
        line = line.strip()
        
        match = re.match(pattern_known, line)
        if match:
            (species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
            proteins_known_data[species_name] = (spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
            
        match = re.match(pattern_abinitio, line)    
        if match:
            (species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
            proteins_abinitio_data[species_name] = (spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    
    descr_file.close()
    
    return proteins_known_data, proteins_abinitio_data
    
def get_protein_ids (protein_id):
    '''
    Parses the description file for the protein_id, and retrieves only the protein ids for every species
    @param protein_id: protein_id for which protein ids of other species should be retrieved
    @return: (prot_ids_known, prot_ids_abinitio) - two dictionaries (key is species name, value is orthologous protein id)
    '''
    (proteins_known, proteins_abinitio) = parse_descr_file(protein_id)
        
    prot_ids_known = {}
    prot_ids_abinitio = {}
    
    for key, value in proteins_known.items():
        prot_ids_known[key] = list(value)[0]
    for key, value in proteins_abinitio.items():
        prot_ids_abinitio[key] = list(value)[0]
        
    return prot_ids_known, prot_ids_abinitio

def get_gene_regions (protein_id):
    '''
    Parses the description file for the protein_id, and retrieves the information about protein locations for every species.
    Locations are stored as tuples (location_type, assembly, location_id, seq_begin, seq_end, strand)
    @param protein_id: protein_id for which protein ids of other species should be retrieved
    @return: (prot_ids_known, prot_ids_abinitio) - two dictionaries (key is species name, value is gene location data as described)
    '''
    (proteins_known, proteins_abinitio) = parse_descr_file(protein_id)
    
    genes_known = {}
    genes_abinitio = {}
    
    for key, value in proteins_known.items():
        genes_known[key] = list(value)[3:]
    for key, value in proteins_abinitio.items():
        genes_abinitio[key] = list(value)[1:]
        
    return genes_known, genes_abinitio

def get_project_root_dir ():
    '''
    Auxiliary function, retrieves the root directory of the entire project
    '''
    ex_path = sys.path[0]
    m = re.match("(.*ExoLocator).*", ex_path)
    proj_root_dir = m.groups()[0]
    return proj_root_dir

def get_species_list ():
    '''
    @return: speacies_list - list of species as available for the species.txt file in the root project directory.
    '''
    species_file_path = get_project_root_dir() + "/species.txt"
    species_file = open(species_file_path, 'r')
    species_list = []
    
    for line in species_file.readlines():
        species_list.append(line.strip()) 
        
    species_file.close()
    return species_list

def get_protein_list ():
    '''
    Parses the protein input file as declared in the command_line_tools.cfg configuration file. 
    Each line of file must contain the protein ID and the number of exons for that particual protein.
    @return: protein_list list of tuples containing two values (protein id, number of exons). 
    '''
    cr = ConfigurationReader.Instance()
    protein_file_path = cr.get_value('input', 'protein_list')
    protein_file = open(protein_file_path, 'r')
    protein_list = []
    
    for line in protein_file.readlines():
        (prot_id, num_of_exons) = line.strip().split()
        protein_list.append((prot_id, num_of_exons))
    protein_file.close()
    
    return protein_list

def read_status_file (protein_id):
    '''
    @return: status_dict dictionary of mapped status entries to their values
    Status entries may be:
        MUTUAL_BEST: OK/FAILED
    '''
    
    dc = DirectoryCrawler()
    status_file_path = dc.get_mutual_best_status_file_path(protein_id)
    try:
        status_file = open(status_file_path, 'r')
    except IOError:
        raise IOError('No .status file for protein %s' % protein_id)
    
    status_dict = dict(token.split() for token in status_file.read().strip().split('\n'))
    
    return status_dict

def update_entry_in_status_file (protein_id, status_entry, status_entry_value):
    '''
    Updates the status entry to new value. 
    If there is no .status file as to this update, it generates the status file.
    If there exists the status file, it reads it.
    If this status entry is already present, and its value the same as the new value, then nothing is done.
    Otherwise, the value is updated and written in the status file.
    '''
    
    dc = DirectoryCrawler()
    status_file_path = dc.get_mutual_best_status_file_path(protein_id)
    status_dict = {}
    
    if (os.path.isfile(status_file_path)):
        status_dict = read_status_file(protein_id)
    
    if (status_dict.has_key(status_entry)):
        if (status_dict[status_entry] == status_entry_value):
            return
        else:
            status_dict[status_entry] = status_entry_value
            status_file = open(status_file_path, 'w')
            for status_entry, status_entry_value in status_dict:
                status_file.write("%s %s\n")
            status_file.close()
            
    else:
        status_file = open(status_file_path, 'a+')
        status_file.write("%s %s\n" % (status_entry, status_entry_value))
        status_file.close()
        

if __name__ == '__main__':
    (prot_ids_known, prot_ids_abinitio) = get_gene_regions("ENSP00000311134")
    print prot_ids_known, prot_ids_abinitio
