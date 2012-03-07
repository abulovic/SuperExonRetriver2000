'''
Created on Mar 6, 2012

@author: Ana
'''

###############################
###############################
# imports
from os import remove, system, path;
import re;
import ConfigParser;
from subprocess import *;


###############################
###############################
#Tools configuration and setting paths
config_file = "../../anab.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
session_resource_path = config.get('Session resource', 'session_resource_path')
#the input file
descr_file_path = config.get('Ensembl cfg', 'descr_file')
descr_file_path = "%s/%s" % (session_resource_path, descr_file_path);
# how much to expand
ens_expansion = config.get('Ensembl cfg', 'expansion')

#the gene regions folder
gene_regions_path = config.get('Gene regions path', 'regions')
gene_regions_path = "%s/%s" % (session_resource_path, gene_regions_path)
# the tmp folder
tmp_folder = config.get('Gene regions path', 'temporary_sequences')
tmp_folder = "%s%s" % (session_resource_path, tmp_folder)
# expanded regions folder
expanded_regions_f = config.get('Gene regions path', 'expanded_regions')
expanded_regions_f = "%s%s" % (session_resource_path, expanded_regions_f)
# cached data status file
status_file = "%s.status" % tmp_folder

###############################
###############################
# Input files
descr_file = open(descr_file_path, 'r')

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
# Function for data caching
def cache_data ():
    i = 0
    
    # remove the cached status
    if (is_data_cached()) :
        remove(status_file)
    
    # start caching
    for line in descr_file_lines:
        if (i%4 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s.fa" % (gene_regions_path, species_name)
        elif (i%4 == 2):
            parameters = line.strip('\n').split(' ')[1].split(':')
            if (parameters[0] != 'chromosome'):
                # otherwise, cache the data
                # after caching, write success to .status file
                cache_f = "%s/%s.%s.fa" % (tmp_folder, species_name, parameters[0])
                cmd = "./ensembl_api.sh -s %s -t %s -i %s -a %s -o %s" %(species_name,       # species name 
                                                                                parameters[0],      # id type 
                                                                                parameters[2],      # id 
                                                                                parameters[1],      # assembly 
                                                                                cache_f)
                system(cmd)
        i = i+1  
    
    # cached status
    stat_f = open(status_file, 'w')
    stat_f.close()
     

###############################
###############################
# Function for checking if the data has been cached
def is_data_cached ():
    if path.isfile(status_file):
        return True
    else :
        return False



###############################
###############################
# function for getting the normal gene regions
def get_gene_regions ():

    i = 0
    
    # if data cached, search the cached files
    # if not, search the database
    print "Get gene regions started"
    cached_stat = is_data_cached()
    if (cached_stat) :
        cmd_cached = "-c %s" % (tmp_folder)
    else :
        cmd_cached = ""
    
    for line in descr_file_lines:
        print line
        if (i%4 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s.fa" % (gene_regions_path, species_name)
        elif (i%4 == 2):
            
            parameters = line.strip('\n').split(' ')[1].split(':')
            # chromosome is never cached
            if (parameters[0] == 'chromosome'):
                cmd = "./ensembl_api.sh -s %s -t %s -i %s -a %s -b %s -e %s -r %s -o %s" %(species_name,       # species name 
                                                                                parameters[0],      # id type 
                                                                                parameters[2],      # id 
                                                                                parameters[1],      # assembly 
                                                                                parameters[3],      # sequence beginning
                                                                                parameters[4],      # sequence ending
                                                                                parameters[5],      # strand
                                                                                output_file_name)   
                system(cmd)
                print "got chromosome"
                
            # other sequences might be
            else:
                cmd = "./ensembl_api.sh -s %s -t %s -i %s -a %s -b %s -e %s -r %s -o %s %s" %(species_name,       # species name 
                                                                                parameters[0],      # id type 
                                                                                parameters[2],      # id 
                                                                                parameters[1],      # assembly 
                                                                                parameters[3],      # sequence beginning
                                                                                parameters[4],      # sequence ending
                                                                                parameters[5],      # strand
                                                                                output_file_name,
                                                                                cmd_cached)
                system(cmd)
                print "got other"
                
        i = i+1
                
                
###############################
###############################
# function for getting the expanded gene regions
def get_expanded_gene_regions (exp = int(ens_expansion)):
    
    i = 0
    
    # if data cached, search the cached files
    # if not, search the database
    cached_stat = is_data_cached()
    if (cached_stat) :
        cmd_cached = "-c %s" % (tmp_folder)
    else :
        cmd_cached = ""
    
    for line in descr_file_lines:
        if (i%4 == 0):
            species_name = line.strip('\n')
            output_file_name = "%s/%s_%d.fa" % (expanded_regions_f, species_name, exp)
        elif (i%4 == 2):
            parameters = line.strip('\n').split(' ')[1].split(':')
            # chromosome is never cached
            if (parameters[0] == 'chromosome'):
                cmd = "./ensembl_api.sh -s %s -t %s -i %s -a %s -b %s -e %s -r %s -o %s" %(species_name,       # species name 
                                                                                parameters[0],      # id type 
                                                                                parameters[2],      # id 
                                                                                parameters[1],      # assembly 
                                                                                max(1, int(parameters[3])-exp),      # sequence beginning
                                                                                int(parameters[4])+exp,              # sequence ending
                                                                                parameters[5],      # strand
                                                                                output_file_name)   
                system(cmd)
                
            # other sequences might be
            else:
                cmd = "./ensembl_api.sh -s %s -t %s -i %s -a %s -b %s -e %s -r %s -o %s %s" %(species_name,       # species name 
                                                                                parameters[0],      # id type 
                                                                                parameters[2],      # id 
                                                                                parameters[1],      # assembly 
                                                                                max(1, int(parameters[3])-exp),      # sequence beginning
                                                                                int(parameters[4])+exp,              # sequence endin
                                                                                parameters[5],      # strand
                                                                                output_file_name,
                                                                                cmd_cached)         # if data is cached, add the -c argument
                system(cmd)
                
        i = i+1
                
#cache_data()
get_gene_regions()
#get_expanded_gene_regions()
