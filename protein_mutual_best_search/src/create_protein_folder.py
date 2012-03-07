'''
Created on Mar 5, 2012

@author: marioot
'''
import re;
import sys;
import ConfigParser;
###############################
###############################
#Returns a list of species short names (eg. HOM_SAP)
def parse_descr_file (descr_file) :
    species_name="";
    short_name="";
    all_species = {};
    
    FILE = open(descr_file,'r');
    while FILE :
        species_name = FILE.readline().rstrip('\n');
        if (not species_name):    # end of file
            break
        short_name = FILE.readline().rstrip('\n');
        params_line = FILE.readline().rstrip('\n').split(' ')[1];
        FILE.readline()
        all_species[short_name] = (species_name, params_line);
             
    FILE.close()
    return all_species
###############################
###############################
#For each species that had a mutual best protein, creates a fasta file (eg. HOM_SAP.fasta) with that protein sequence.
#It stores it in ../session_resource/mutual_best_proteins/ folder.
def create_protein_files(all_species, known_proteome, abinitio_proteome, proteins_path): 
    for species in all_species:
        path_in = known_proteome + all_species[species][0] + "/" + all_species[species][0]
        path_out = proteins_path + species + ".fasta"

        INFILE = open(path_in,'r')        
        lines = INFILE.read()
        short_name_line = lines.find(all_species[species][1])
        if short_name_line == -1:
            path_in = abinitio_proteome + all_species[species][0] + "/" + all_species[species][0]
            INFILE = open(path_in,'r')
            lines = INFILE.read()
            short_name_line = lines.find(all_species[species][1])
        data = ""
        while(lines[short_name_line] != ">"):
            data += lines[short_name_line]
            short_name_line += 1
        data_split = data.splitlines()
        data_split[0] = ">" + species + " " + all_species[species][1] + "\n"
        OUTFILE = open(path_out, 'w')
        for d in data_split:
            OUTFILE.write(d)
###############################
###############################
#  Tools configuration and setting paths
config_file = "marioot.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
known_proteome_f = config.get('Database path', 'known_proteome_path')
abinit_proteome_f = config.get('Database path', 'abinit_proteome_path')
#session resource
session_resource_path = config.get('Session resource', 'session_resource_path')
#resulting proteins path
proteins_path = session_resource_path + config.get('Found proteins path', 'proteins')
###############################
#  input negotiation
#
#descr = "test.descr";
if(len(sys.argv) < 1):
    print "Usage: %s <descr file>\n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
descr = sys.argv[1] + ".descr"
descr = session_resource_path + descr
###############################
#  Create protein folder
all_species = parse_descr_file(descr)
create_protein_files(all_species, known_proteome_f, abinit_proteome_f, proteins_path)