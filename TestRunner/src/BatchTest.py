'''
Created on Mar 13, 2012

@author: Ana
'''

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

def generate_directories_tree ():
    for protein_id in protein_file.readlines():
        protein_id = protein_id.strip()
        protein_sessions_dir = "%s/%s" % (session_folder, protein_id);
        
        if (not os.path.isdir(protein_sessions_dir)):
            os.makedirs(protein_sessions_dir)
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, gene_regions_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, gene_regions_folder))
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, expanded_regions_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, expanded_regions_folder))
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, exons_path))):
            os.mkdir("%s/%s" % (protein_sessions_dir, exons_path))
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, mut_best_proteins))):
            os.makedirs("%s/%s" % (protein_sessions_dir, mut_best_proteins))
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, blastout_folder))
        if (not os.path.isdir("%s/%s/dna" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s/dna" % (protein_sessions_dir, blastout_folder))
        if (not os.path.isdir("%s/%s/protein" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s/protein" % (protein_sessions_dir, blastout_folder))
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

gene_regions_folder     = config.get('Gene regions path', 'regions')
expanded_regions_folder = config.get('Gene regions path', 'expanded_regions')
exons_path              = config.get('Exon database path', 'exons')
mut_best_proteins       = config.get('Found proteins path', 'proteins')
blastout_folder         = config.get('Blastout path', 'blastout')

fasta_db                = config.get('Ensembl cfg', 'ensembldb')

###################################

species         = "Homo_sapiens"
protein_type    = "all"

protein_database = generate_file_name(protein_type, species);

protein_file = open(protein_list_file, 'r')

generate_directories_tree()

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
    #os.system(cmd_retrieve_protein)
    
    cmd_run_mutual_best = "python ../../protein_mutual_best_search/src/ensembl_mutual_best.py %s %s" % (species, protein_session_dir)
    # run the mutual best protein search
    #os.system(cmd_run_mutual_best)
    
    #get the dna data
    cmd_retrieve_dna = "python ../../ensembl_search/src/grab_slices.py %s" % protein_session_dir
    #os.system(cmd_retrieve_dna)
    
    #get the base exons
    cmd_generate_exons = "python ../../exon_finder/src/exon_base_generator.py %s %s %s" % (input_protein, 
                                                                                           "%s/gene_regions/%s.fa" % 
                                                                                           (protein_session_dir, species),
                                                                                          protein_session_dir)
    number_of_exons = os.system(cmd_generate_exons)>>8

    #run exon_finder
    cmd_run_exon_finder = "python ../../exon_finder/src/exon_finder.py %s %s" % (number_of_exons, protein_session_dir, protein_id)
    os.system(cmd_run_exon_finder)
    protein_file.close()