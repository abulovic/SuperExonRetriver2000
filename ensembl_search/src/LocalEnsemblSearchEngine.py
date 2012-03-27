'''
Created on Mar 6, 2012

@author: Ana
'''

###############################
###############################
# imports
from os import system, path, makedirs, listdir;
import re;
import sys;
import fileinput;
import ConfigParser;

###############################
# check for the command line
# arguments
if (len(sys.argv) < 2 ):
    print ("Too few arguments. Pass the path to resource folder as an argument and try again.")
else:
    session_resource_path = sys.argv[1]


###############################
###############################
#Tools configuration and setting paths
config_file = "../../config.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
#the input file
descr_file_path = config.get('Session files', 'descr_output')

descr_file_path = "%s/%s" % (session_resource_path, descr_file_path);
ensembldb = config.get('Ensembl cfg', 'ensembldb')
masked = config.get('Ensembl cfg', 'masked')
if (masked == 1): 
    masked = True
else:
    masked = False
# how much to expand
ens_expansion = int(config.get('Ensembl cfg', 'expansion'))

#the gene regions folder
gene_regions_path = config.get('Gene regions path', 'regions')
gene_regions_path = "%s/%s" % (session_resource_path, gene_regions_path)
# expanded regions folder
expanded_regions_f = config.get('Gene regions path', 'expanded_regions')
expanded_regions_f = "%s/%s" % (session_resource_path, expanded_regions_f)
# protein path
protein_folder = config.get('Found proteins path', 'proteins')
protein_folder = "%s/%s" % (session_resource_path, protein_folder)

exon_folder = config.get('Exon database path', 'exons_path')
exon_folder = "%s/%s" % (protein_folder, exon_folder)

exon_database = "%s/db/" % (exon_folder)


###############################
# create session dir, if not existing
if (not path.exists(session_resource_path)):
    makedirs(session_resource_path)
    
if (not path.exists(gene_regions_path)):
    makedirs(gene_regions_path)
    
if (not path.exists(expanded_regions_f)):
    makedirs(expanded_regions_f)
    
if (not path.exists(protein_folder)):
    makedirs(protein_folder)
    
if (not path.exists(exon_folder)):
    makedirs(exon_folder)
    
if (not path.exists(exon_database)):
    makedirs(exon_database)
    

    


    

###############################
###############################
# Input files
descr_file = open(descr_file_path, 'r')
print descr_file

###############################
###############################
# dictionaries
ensembl_ids = dict()        #keys - species name, values - all the neccessary information

###############################
###############################
# input file parsing
descr_file_lines = descr_file.readlines()
# all the variables
# mandatory:
species_name = ""
id_type = ""
sid = ""
# optional:
assembly = ""
beginning = 0
ending = 0
strand = 1

###############################
###############################
def generate_file_name (masked, species, id_type, id):
    file_name = "%s/%s/dna/" % (ensembldb, species.lower())
    tmp_file=""
    for file in listdir(file_name):
        if (file != "README"):
            tmp_file = file 
            break
    m = re.findall ('(.*).dna', tmp_file)   
    if (masked == True):
        file_name = "%s/%s.dna_rm." % (file_name, m[0])
    else :
        file_name = "%s/%s.dna." % (file_name, m[0])
    if (id_type == 'chromosome'):
        file_name = "%schromosome.%s.fa" % (file_name, id)
    else :
        file_name = "%stoplevel.fa" % (file_name)
    return file_name
   
################################
# protein_type - all / abinitio
def generate_protein_file_name (species, protein_type):
    file_name = "%s/%s/pep" % (ensembldb, species.lower())
    tmp_file=""
    for file in listdir(file_name):
        if (file != "README"):
            tmp_file = file 
            break
    m = re.findall ('(.*).pep', tmp_file)
    if (protein_type == "all"):
        file_name = "%s/%s.pep.all.fa" % (file_name, m[0])
    else:
        file_name = "%s/%s.pep.abinitio.fa" % (file_name, m[0])
        
    return file_name

def add_sequence_end_tags(output_file_name, seq_begin, seq_end):
    les_stat = "F"
    res_stat = "F"
    for line in fileinput.input(output_file_name, inplace=1):
        if (fileinput.filelineno() == 1):
            m = re.findall(":1:([0-9]+):-*1", line)
            seq_len = int(m[0])
            if (seq_begin - ens_expansion < 1):
                les_stat = "T"
            if (seq_end + ens_expansion > seq_len):
                res_stat = "F"
            print "%s les:%s res:%s" % (line.strip(), les_stat, res_stat)
        else:
            print line,

###############################
###############################
# function for getting the normal gene regions
def get_gene_regions ():

    i = 0
    print "Get gene regions started."    
    for line in descr_file_lines:
        print line
        if (i%3 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s.fa" % (gene_regions_path, species_name)
            
        elif (i%3 == 1):
            
            parameters = line.strip('\n').split(' ')[1].split(':')
            database = generate_file_name(masked, species_name, parameters[0], parameters[2])
            print database
            if (int(parameters[5]) == -1) : # convert into fastacmd friendly format
                parameters[5] = str(2)
            sid = ""
            if (parameters[0] == "chromosome"):
                sid = "chrom%s" % parameters[2]
            else :
                sid = parameters[2]
                
            cmd = "fastacmd -d %s -s %s -S %s -L %s,%s -p F -o %s"  %       (database,          # database name
                                                                             sid,                # id
                                                                             parameters[5],     # strand
                                                                             parameters[3],     # seq beginning
                                                                             parameters[4],     # seq ending
                                                                             output_file_name)
            print (cmd)   
            system(cmd)
            print "Wrote data to %s" % output_file_name
        i = i+1
                
                
###############################################
###############################################
# function for getting the expanded gene regions
def get_expanded_gene_regions ():

    i = 0
    print "Get expanded gene regions started."    
    for line in descr_file_lines:
        print line
        if (i%3 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s.fa" % (expanded_regions_f, species_name)
            
        elif (i%3 == 1):
            
            parameters = line.strip('\n').split(' ')[1].split(':')
            database = generate_file_name(masked, species_name, parameters[0], parameters[2])
            print database
            if (int(parameters[5]) == -1) : # convert into fastacmd friendly format
                parameters[5] = str(2)
            sid = ""
            if (parameters[0] == "chromosome"):
                sid = "chrom%s" % parameters[2]
            else :
                sid = parameters[2]
                
            cmd = "fastacmd -d %s -s %s -S %s -L %d,%d -p F -o %s"  %       (database,          # database name
                                                                             sid,               # id
                                                                             parameters[5],     # strand
                                                                             max(0, int(parameters[3])-ens_expansion),     # seq beginning
                                                                             int(parameters[4])+ens_expansion,     # seq ending
                                                                             output_file_name)   
            system(cmd)
            add_sequence_end_tags(output_file_name, int(parameters[3]), int(parameters[4]))
            print "Wrote data to %s" % output_file_name
        i = i+1
        
        
def get_proteins ():
    i = 0
    print "Get proteins started."    
    for line in descr_file_lines:
        if (i%3 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s.fa" % (protein_folder, species_name)
        
        elif (i%3 == 1):
            prot_id = line.strip().split()[0]
            protein_type = line.strip().split()[2]
            if (protein_type == "known" ):
                database = generate_protein_file_name(species_name, "all")
            else:
                database = generate_protein_file_name(species_name, "abinitio")
                
            cmd = "fastacmd -d %s -s %s -p T -o %s"  %       (database,          # database name
                                                             prot_id,            # id
                                                             output_file_name)
            print(cmd)
            system(cmd)
            
        i = i + 1
            
            
            
            
get_gene_regions()
get_expanded_gene_regions()
get_proteins()
