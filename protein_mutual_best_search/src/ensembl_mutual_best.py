'''
Created on Oct 9, 2011

@author: Mario, Ana :)
'''
from os import system, popen, remove, listdir;
import re;
import sys;
import ConfigParser;
from find_by_blasting import find_by_blasting, extract_first_seq;

###############################
###############################
## function to do the bidirectional search
def bidirectional_search (protein_type,     # all / abinitio
                          species_name,     # in Latin
                          original_sequence_fa, #fasta file with the original protein sequence
                          ):
    match_found = False
    logger("%s\n" % species_name)
    ## do the forward search
    
    ID_gene = ""
    ID_transcript = ""
    
    proteome_database = generate_file_name (protein_type, species_name)
    
    
    forward_ids = find_by_blasting(proteome_database, original_sequence_fa, working_results_f)
    if (len(forward_ids) == 0):
        logger("%s not found in %s, \"known\" sequences\n\n" % (original_sequence_fa, species_name))
        if (protein_type == "all"):
            return (False, "", "", None, None)
        else:
            return (False, "", "")
    else:
        # if we found something, we still need to check if it is the mutual best match
        #print forward_ids
        logger("best forward hit:  %s  \n\n" % forward_ids[0])
        
        
        if (protein_type == "all"):
            [ID_protein, ID_gene, ID_transcript, gene_location] = forward_ids[0].split()[0:4]
        else :
            [ID_protein, ID_protein2, gene_location] = forward_ids[0].split()[0:3]
        extract_first_seq (working_results_f, ID_protein, "%s.output_fasta" % species_name)
        
        
        #write forward protein to out file
        forward_fasta_f.write(">%s_FWD\n" % species_name)
        forward_fasta_f.write(popen("grep -v \'>\' %s.output_fasta" % species_name).read())
        
        ###back search###
        back_ids = find_by_blasting(original_proteome_f, "%s.output_fasta" % species_name, working_results_f)
        if (len(back_ids) == 0):
            logger("\t reciprocal search in %s using %s as a query produced no hits (?)\n\n" % (original_species_name, species))
            logger("*****************************************\n\n\n")
            return False
        
        # check whether back search retrieves the original query
        # three possibilities here
        # 1) the first hit is the original gene ==> this is the mutual best hit
        # 2) the original gene does not appear on the list -- we take there is no ortohologue in this species
        # 3) the original gene is somewhere down the list  -- for now take that it also means "orthologue not found"
        [ID_protein_b, ID_gene_b, ID_transcript_b, gene_location_b] = back_ids[0].split()[0:4]
        
        
        # if the back gene is the same as the original one, we are done
        # otherwise, go and check in the "ab initio" detected list
        if (ID_gene_b == orig_gene_name):

            species_data = popen("grep -v \'>\' %s.output_fasta" % species_name).read();
            remove("{0}.output_fasta".format(species_name))
            fasta_f.write(">%s\n%s" % (species_name, species_data))
            
            # for chromosomal dna, write to chr.fa file
            dna_type = gene_location.split(":")[0]
            if (dna_type == "chromosome"):
                chromosome_fasta_f.write(">%s pep:%s gene:%s transcript:%s %s\n%s" % (species_name, ID_protein, ID_gene, ID_transcript, gene_location, species_data ))
            
            if (protein_type == "all"):
                logger("\t best mutual hit:  %s  %s\n\n" % (ID_gene, ID_protein))
            else :
                logger ("\t best mutual hit:  %s\n\n" % (ID_protein))
            match_found = True
        else:
            logger("\t back search returns  %s as the best hit\n\n" % ID_gene_b)
            if (protein_type == "all"):
                logger("\t no mutual best found in known -- try ab initio \n\n")
                
    if (protein_type == "all"):
        return (match_found, ID_protein, gene_location, ID_gene, ID_transcript)
    else:
        return (match_found, ID_protein, gene_location)
    
   
########################################
# invoke mafft and do msa on the results
# of the bidirectional search 
def mafft_alignment (input_fasta_file):
    

    afafile = re.sub(".fa", ".mafft.afa", input_fasta_file.name)
    cmd = "%s --quiet %s > %s" % (mafft, input_fasta_file.name, afafile)
    print cmd
    system(cmd);
    
###############################
###############################
#Logs the message
def logger(message):
    print message
    LOG.write(message)
    
################################
################################
# generate file name
# prot_type abinitio / all
# species species name in Latin (homo_sapiens)
def generate_file_name (prot_type, species):
    file_name = "%s/%s/pep/" % (ensembldb, species.lower())
    generic_file_name=""
    # extrapolate the generic part of the file name
    # (Species.assembly_id.database_version).pep.prot_type.fa
    for f in listdir(file_name):
        if (f != "README"):
            generic_file_name = f
            break
    m = re.findall ('(.*).pep.', generic_file_name)   
    file_name = "%s/%s.pep.%s.fa" % (file_name, m[0], prot_type)
    return file_name
        
              
###############################
###############################
#Tools configuration and setting paths
config_file = "../../config.cfg"
config          = ConfigParser.RawConfigParser()
config.read(config_file)

project_root_dir        = config.get('Project root', 'project_root_folder')
species_list            = config.get('Project root', 'species_list')
species_list            = "%s/%s" % (project_root_dir, species_list)

ensembldb               = config.get('Ensembl cfg', 'ensembldb')

output_fasta            = config.get('Session files', 'fasta_output')
output_descr            = config.get('Session files', 'descr_output')
protein_input           = config.get('Session files', 'input_protein_file')

mafft                   = config.get('mafft cfg', 'mafft')

log_file_name           = config.get('LOG', 'mutual_best')

###############################
# temporary results
working_results_f       = "tmp.output_fasta"

species_list            = popen("cat %s" % species_list).read().split("\n")


###############################
#  input negotiation
#
#output_fasta = "test.output_fasta";
#output_descr = "test.output_descr";

if(len(sys.argv) < 2):
    print "ERROR:\tInvalid number of input arguments.\n"
    print "Usage: Recquired command line arguments:"
    print "1) Species name in Latin (underscore delimited)"
    print "2) session resource folder"
    exit
    
[original_species_name, session_resource_path] = [sys.argv[1], sys.argv[2]]

log_path =      "%s/%s" % (session_resource_path, log_file_name)
output_descr =  "%s/%s" % (session_resource_path, output_descr)
output_fasta =  "%s/%s" % (session_resource_path, output_fasta)
protein_input_file = "%s/%s" % (session_resource_path, protein_input)

###############################
#  output files
descr_f = open(output_descr, 'w')
fasta_f = open(output_fasta, 'w')
forward_fasta = re.sub(".fa", ".fwd.fa", output_fasta)
forward_fasta_f = open(forward_fasta, 'w')
chromosome_fasta = re.sub(".fa", ".chr.fa", output_fasta)
chromosome_fasta_f = open(chromosome_fasta, 'w')
LOG = open(log_path, "w")

###############################
#  dictionaries
forward_ids_dict = dict()
ab_init_forward_ids_dict = dict()
seen = dict()

###############################
#  other initialization
output_gene_location = ""
output_protein = ""


###############################
#  find the  query sequence in the "original" species
#  (together with the  gene/protein/transcript entry it belongs to)
original_proteome_f     = generate_file_name("all", original_species_name)
orig_search_ids         = find_by_blasting(original_proteome_f, protein_input_file, working_results_f)
if (len(orig_search_ids) == 0):
    logger("No original species ID retrieved by blasting.\n")
    exit (-1)
logger("the closest match to %s in  %s is %s\n" % (protein_input_file, original_species_name, orig_search_ids[0]))

split_info              = re.split(" ", orig_search_ids[0])
[orig_protein, orig_gene_name, orig_transcript, orig_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
logger("%s id in %s: %s\n" % (protein_input_file, original_species_name, orig_protein))

logger("gene name in %s:   %s\n\n" % (original_species_name, orig_gene_name))
orig_seq_f               = "orig_seq.output_fasta"
extract_first_seq(working_results_f, orig_protein, orig_seq_f)

#print bidirectional_search("all", "Bos_taurus", protein_input_file)

 
    
###############################
#  blast the query sequence against all of the remaining proteomes
#
not_found = []

try:
    for species in species_list:
        (match_found_known, ID_protein_known, gene_location_known, ID_gene, ID_transcript) = bidirectional_search("all", species, protein_input_file)
        print "Found in known: " + str(match_found_known)
        if (match_found_known):
            descr_f.write("%s\n" % (species))
            descr_f.write("%s %s %s %s %s\n\n" % (ID_protein_known, gene_location_known, "known", ID_gene, ID_transcript))
        else:
            (match_found_abinitio, ID_protein_abinitio, gene_location_abinitio) = bidirectional_search("abinitio", species, protein_input_file)
            print "Found in abinitio: " + str(match_found_abinitio)
            if (match_found_abinitio):
                descr_f.write("%s\n" % (species))
                descr_f.write("%s %s %s\n\n" % (ID_protein_abinitio, gene_location_abinitio, "abinitio"))
            else:
                not_found.append(species)
                forward_ids_dict[species] = (ID_protein_known, gene_location_known)
                ab_init_forward_ids_dict[species] = (ID_protein_abinitio, gene_location_abinitio)
                logger("\t\t no mutual best found in ab initio either\n\n\n")
except (Exception):
    print "Exception in ensembl_mutual_best, somewhere out there, exception {0}".format(sys.exc_info())
            

fasta_f.close()
descr_f.close()
chromosome_fasta_f.close()


mafft_alignment(fasta_f)
mafft_alignment(chromosome_fasta_f)
mafft_alignment(forward_fasta_f)
###############################
###############################
#Summary of proteins that were not found
logger("not found:\n")

for item in not_found:
    logger("\t%s\n \tknown\n\n" % item)
    logger("\tprotein ID: %s\ngene location: %s" % forward_ids_dict[item])
    logger("\t%s\n \tabinitio\n\n" % item)
    logger("\tprotein ID: %s\ngene location: %s" % ab_init_forward_ids_dict[item])

    
###############################
###############################
#cleanup
remove(orig_seq_f)
remove(working_results_f)
###############################
###############################
# optional: sort and align


