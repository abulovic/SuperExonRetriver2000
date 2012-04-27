'''
Created on Apr 15, 2012

@author: intern
'''

import re, sys, os

from utilities.ConfigurationReader import ConfigurationReader
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

from subprocess import Popen, PIPE, STDOUT
import shutil

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

def get_status_file_keys():
    status_file_keys_path = "{0}/cfg/status_file_keys.txt".format(get_project_root_dir())
    status_file_keys      = open(status_file_keys_path, 'r')
    keys                  = []
    for line in status_file_keys.readlines():
        keys.append(line.strip())
    return keys

def get_reference_species_dictionary():
    reference_species_mapping_path  = "{0}/cfg/referenced_species_mapping.txt".format(get_project_root_dir())
    reference_species_mapping       = open(reference_species_mapping_path, 'r')
    ref_species                     = {}
    for line in reference_species_mapping.readlines():
        key_val = line.split("-")
        ref_species[key_val[0]] = key_val[1].strip()
    return ref_species

def read_status_file (protein_id):
    '''
    @return: status_dict dictionary of mapped status entries to their values
    Status entries may be:
        MUTUAL_BEST:    OK/FAILED
        DATA_RETRIEVAL: OK/PARTIAL/FAILED
    '''
    
    dc = DirectoryCrawler()
    status_file_path = dc.get_mutual_best_status_file_path(protein_id)
    try:
        status_file = open(status_file_path, 'r')
    except IOError:
        raise IOError('No .status file for protein %s' % protein_id)
    
    status_dict = dict(token.split() for token in status_file.read().strip().split('\n'))
    status_file.close()
    
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
            for status_entry, status_entry_value in status_dict.items():
                status_file.write("%s %s\n" % (status_entry, status_entry_value))
            status_file.close()
            
    else:
        status_file = open(status_file_path, 'a+')
        status_file.write("%s %s\n" % (status_entry, status_entry_value))
        status_file.close()
        
def execute_command_and_log (logger, command, arguments = None):
    '''
    Arguments (if provided) should have the format of a list containing (species, protein_id)
    '''
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if output:
        output = output.strip()
        output = " ".join(output.split('\n'))
        if arguments:
            arguments = list(arguments)
            arguments.append(output)
            arguments.append(command)
            logger.error("%s,%s,%s,%s" % tuple(arguments))
        else:
            logger.error(",,%s,%s" % (output, command))

def clear_directory (path):
    for root, dirs, files in os.walk(path):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))
            
def reset_action (protein_id, key):
    update_entry_in_status_file(protein_id, key, 'FAILED')
    crawler = DirectoryCrawler()
    
    if key == 'GENE_RETRIEVAL': 
        clear_directory(crawler.get_gene_path(protein_id))
    elif key == 'EXP_GENE_RETRIEVAL' : 
        clear_directory(crawler.get_expanded_gene_path(protein_id))
    elif key == 'PROTEIN_RETRIEVAL' : 
        clear_directory(crawler.get_protein_path(protein_id))
    elif key == 'ENSEMBL_EXON_RETRIEVAL' : 
        clear_directory(crawler.get_exon_ensembl_path(protein_id))
    elif key == 'GENEWISE_EXON_RETRIEVAL' : 
        clear_directory(crawler.get_exon_genewise_path(protein_id))
        clear_directory(crawler.get_genewise_path(protein_id))
    elif key == 'REF_SP_DB_FORMATTING' : 
        clear_directory(crawler.get_database_path(protein_id))
    elif key == 'BLASTN_ALIGNMENT' : 
        clear_directory(crawler.get_blastn_path(protein_id))
    elif key == 'TBLASTN_ALIGNMENT' : 
        clear_directory(crawler.get_tblastn_path(protein_id))
    elif key == 'SW_GENE_ALIGNMENT' : 
        clear_directory(crawler.get_SW_gene_path(protein_id))
    elif key == 'SW_EXON_ALIGNMENT' : 
        clear_directory(crawler.get_SW_exon_path(protein_id))

def reset_action_global(key):
    protein_list_raw = get_protein_list()
    protein_list = []
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
        
    for protein in protein_list:
        reset_action(protein, key)
        
def main():
    print get_status_file_keys()
    print get_reference_species_dictionary()

if __name__ == '__main__':
    main()
