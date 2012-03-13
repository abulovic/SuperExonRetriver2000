'''
Created on Mar 13, 2012

@author: Ana
'''

import sys
import re
import os
import ConfigParser;

################################
################################
# generate file name
# prot_type abinitio / all
# species species name in Latin (homo_sapiens)
def generate_file_name (prot_type, species):
    file_name = "%s/%s/pep/" % (fasta_db, species.lower())
    generic_file_name=""
    # extrapolate the generic part of the file name
    # (Species.assembly_id.database_version).pep.prot_type.fa
    for f in os.listdir(file_name):
        if (f != "README"):
            generic_file_name = f
            break
    m = re.findall ('(.*).pep.', generic_file_name)   
    file_name = "%s/%s.pep.%s.fa" % (file_name, m[0], prot_type)
    return file_name

##################################
## Load necessary configuration ##

config_file     = "../../config.cfg"
config          = ConfigParser.RawConfigParser()
config.read(config_file)

project_root_folder     = config.get('Project root', 'project_root_folder')
session_folder          = config.get('Project root', 'session_folder')
session_folder          = "%s/%s" % (project_root_folder, session_folder)

protein_list_file       = config.get('Test files', 'protein_list')
input_protein_file      = config.get('Session files', 'input_protein_file')

fasta_db                = config.get('Ensembl cfg', 'ensembldb')

###################################

species         = "Homo_sapiens"
protein_type    = "all"

protein_database = generate_file_name(protein_type, species);

protein_file = open(protein_list_file, 'r')

for protein_id in protein_file.readlines():
    protein_id = protein_id.strip()
    protein_session_dir = "%s/%s" % (session_folder, protein_id)
    if (not os.path.isdir(protein_session_dir)):
        os.makedirs(protein_session_dir)
    input_protein = "%s/%s" % (protein_session_dir, input_protein_file)
    print input_protein
    
    cmd_retrieve_protein = "fastacmd -d %s -s %s -p T -o %s"  %  (protein_database,          # database name
                                                                 protein_id,                # id
                                                                 input_protein)
    
    # extract the protein from the database
    os.system(cmd_retrieve_protein)
    
    cmd_run_mutual_best = "python ../../protein_mutual_best_search/src/ensembl_mutual_best.py %s %s" % (species, protein_session_dir)
    # run the mutual best protein search
    #os.system(cmd_run_mutual_best)
    
    #get the dna data
    cmd_retrieve_dna = "python ../../ensembl_search/src/grab_slices.py %s" % protein_session_dir
    os.system(cmd_retrieve_dna)
    
    
    
    


