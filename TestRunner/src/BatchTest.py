'''
Created on Mar 13, 2012

@author: anana
'''



import re
import os, sys, inspect

sys.path.append("../../ensembl_search/src")
sys.path.append("../../exon_finder/src")

import ConfigParser;
from BioMartSearchEngine import *
from AlignmentGenerator import *
from AlignmentParserBlast import *
from AlignmentParserSW import *
from StatisticsGenerator import StatisticsGenerator



def parseDescriptionFile (descrFileName):
    descrFile = open(descrFileName, 'r')
    abinitioSpecies = []
    knownSpecies = []
    i = 0
    for line in descrFile.readlines():
        line = line.strip()
        
        if (len(line) == 0):
            continue
        if (i % 2 == 0):
            species= line
        else :
            data = line.split()
            if (len(data) == 5):
                knownSpecies.append(species)
            else :
                abinitioSpecies.append(species)
        i = i + 1
        
    
    descrFile.close()
    return (knownSpecies,abinitioSpecies)
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
    protein_file = open(protein_list_file, 'r')
    for line in protein_file.readlines():
        line = line.strip()
        if (len(line) == 0):
            continue
        (proteinId, numberOfExons) = line.split()
        numberOfExons = int(numberOfExons)
        protein_sessions_dir = "%s/%s" % (session_folder, proteinId);
        
        print "%s/%s/dna" % (protein_sessions_dir, swout_folder)
        
        if (not os.path.isdir(protein_sessions_dir)):
            os.makedirs(protein_sessions_dir)
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, gene_regions_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, gene_regions_folder))
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, expanded_regions_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, expanded_regions_folder))
            
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, exons_path))):
            os.mkdir("%s/%s" % (protein_sessions_dir, exons_path))
        if (not os.path.isdir("%s/%s/db" % (protein_sessions_dir, exons_path))):
            os.mkdir("%s/%s/db" % (protein_sessions_dir, exons_path))
        if (not os.path.isdir("%s/%s/wise2" % (protein_sessions_dir, exons_path))):
            os.mkdir("%s/%s/wise2" % (protein_sessions_dir, exons_path))
            
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, mut_best_proteins))):
            os.makedirs("%s/%s" % (protein_sessions_dir, mut_best_proteins))
            
        if (not os.path.isdir("%s/%s" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s" % (protein_sessions_dir, blastout_folder))
        if (not os.path.isdir("%s/%s/dna" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s/dna" % (protein_sessions_dir, blastout_folder))
        if (not os.path.isdir("%s/%s/protein" % (protein_sessions_dir, blastout_folder))):
            os.makedirs("%s/%s/protein" % (protein_sessions_dir, blastout_folder))
            
        if (not os.path.isdir("%s/%s/" % (protein_sessions_dir, swout_folder))):
            os.makedirs("%s/%s/" % (protein_sessions_dir, swout_folder))
        if (not os.path.isdir("%s/%s/dna" % (protein_sessions_dir, swout_folder))):
            os.makedirs("%s/%s/dna/" % (protein_sessions_dir, swout_folder))
        if (not os.path.isdir("%s/%s/exon" % (protein_sessions_dir, swout_folder))):
            os.makedirs("%s/%s/exon/" % (protein_sessions_dir, swout_folder))
    protein_file.close()
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
descr_file              = config.get('Session files', 'descr_output')

gene_regions_folder     = config.get('Gene regions path', 'regions')
expanded_regions_folder = config.get('Gene regions path', 'expanded_regions')
exons_path              = config.get('Exon database path', 'exons_path')
mut_best_proteins       = config.get('Found proteins path', 'proteins')
blastout_folder         = config.get('Blastout path', 'blastout')
swout_folder            = config.get('SWout path', 'swout')

fasta_db                = config.get('Ensembl cfg', 'ensembldb')

statistics_path         = config.get('Statistics', 'exon_finder')

###################################

species         = "Homo_sapiens"
protein_type    = "all"

protein_database = generate_file_name(protein_type, species);

generate_directories_tree()

biomart         = BioMartSearchEngine()
alignmentGen    = AlignmentGenerator()
alParserBlast   = AlignmentParserBlast()
alParserSW      = AlignmentParserSW()
statGen         = StatisticsGenerator()


protein_file = open(protein_list_file, 'r')
for line in protein_file.readlines():
    line = line.strip()
    if (len(line) == 0):
            continue
    (protein_id, numberOfExons) = line.split()
    numberOfExons = int(numberOfExons)
    
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
    
    try:
        alignmentGen.setProteinFolder(protein_id)
        alParserBlast.setProteinFolder(protein_id)
        alParserSW.setProteinFolder(protein_id)
        statGen.setProteinFolder(protein_id)
    except (Exception):
        print "Exception in setting protein folder, protein {0}, exception {1}".format(protein_id, sys.exc_info())
    
    cmd_run_mutual_best = "python ../../protein_mutual_best_search/src/ensembl_mutual_best.py %s %s" % (species, protein_session_dir)
    #run the mutual best protein search
    mutual_best_return_value = os.system(cmd_run_mutual_best)>>8
    if (mutual_best_return_value == -1):
        print ("For protein %s no original protein was retrieved by blastp." % protein_id)
        continue
    
    #get the dna data
    cmd_retrieve_dna = "python ../../ensembl_search/src/LocalEnsemblSearchEngine.py %s" % protein_session_dir
    os.system(cmd_retrieve_dna)
    
    #get the exons

    try:
        biomart.populateExonDatabase(protein_id)
    except (Exception):
        print "Exception in populateExonDatabase, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    try:
        (known_species, abinitioSpecies) = parseDescriptionFile("%s/%s/%s" % (session_folder, protein_id, descr_file))
    except (Exception):
        print "Exception in parse_description_file, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    try:
        alignmentGen.runBatchBlastn(True)
        alignmentGen.runBatchTblastn()
    except (Exception):
        print "Exception in run batch blast, proiten {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    try:
        alignmentGen.runBatchSW()
        alignmentGen.runBatchSW(swType="ex-ex")
    except (Exception):
        print "Exception in run batch SW, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    try:
        alignmentGen.runBatchGenewise(abinitioSpecies)
    except (Exception):
        print "Exception in run batch genewise, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    try:
        exonsBlastn = alParserBlast.batchParseOutput(numberOfExons, "blastn")
    except (Exception):
        print "Exception in batchParseOutput blastn, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    try:
        exonsTblastn = alParserBlast.batchParseOutput(numberOfExons, "tblastn")
    except (Exception):
        print "Exception in batchParseOutput tblastn, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
        
    try:
        exonsSW = alParserSW.batchParseOutputExDna()
    except (Exception):
        print "Exception in SW batchParseOutput, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    try:
        exonsSWEnsembl = alParserSW.batchParseOutputExEx()
    except (Exception):
        print "Exception in SW batchParse output ex-ex, protein {0}, exception {1}".format(protein_id, sys.exc_info())

    try:
        statGen.generate_statistics(exonsTblastn, exonsBlastn, exonsSW, exonsSWEnsembl, known_species+abinitioSpecies, protein_id)
    except (Exception):
        print "Exception in generate_statistics, protein {0}, exception {1}".format(protein_id, sys.exc_info())
        continue
    
    #get the base exons
    '''
    for species in abinitioSpecies:
        cmd_generate_exons = "python ../../exon_finder/src/exon_base_generator.py %s %s %s %s" % ( "%s/%s/%s.fa" % (protein_session_dir, mut_best_proteins, species),
                                                                                                "%s/gene_regions/%s.fa" %  (protein_session_dir, species),
                                                                                                protein_session_dir,
                                                                                                species)
        number_of_exons = os.system(cmd_generate_exons)>>8
        print "Found %d exons for species %s" % (number_of_exons, species)
    generateExonFastaFromDescription(protein_session_dir)
    '''
    
    #print "%s %s %s" % (number_of_exons, protein_session_dir, protein_id)
    #run exon_finder
    #cmd_run_exon_finder = "python ../../exon_finder/src/exon_finder.py %s %s %s" % (number_of_exons, protein_session_dir, protein_id)
    #os.system(cmd_run_exon_finder)
