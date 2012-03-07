'''
Created on Oct 9, 2011

@author: Mario
'''
from os import system, popen, remove;
import re;
import sys;
import ConfigParser;
from find_by_blasting import find_by_blasting, extract_first_seq;

###############################
###############################
#example: HOM_SAP
def short_name(name):
    aux = re.split("_", name)
    aux[0] = aux[0].upper()
    aux[1] = aux[1].upper()
    if(name_tag == ""):
        return aux[0][0:3]+"_"+aux[1][0:3]
    else:
        return aux[0][0:3]+"_"+aux[1][0:3]+"_"+name_tag
    

###############################
###############################
#Logs the message
def logger(message):
    print message
    LOG.write(message)
    


###############################
###############################
#Tools configuration and setting paths
config_file = "../../anab.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
ensembldb = config.get('Ensembl cfg', 'ensembldb')
descr = config.get('Ensembl cfg', 'descr_file')
species_list = '../species.txt' 
working_results_f = "tmp.fasta"
species_list = popen("cat %s" % species_list).read().split("\n")
#session resource
session_resource_path = config.get('Session resource', 'session_resource_path')
log_path = session_resource_path + config.get('LOG', 'mutual_best')
###############################
#  input negotiation
#
#fasta = "test.fasta";
#descr = "test.descr";
if(len(sys.argv) < 3):
    print "Usage: %s <input seq file> <orig_genome>  <fasta file (output)> [<name tag>] \n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
[seqfile, orig_genome, fasta] = [sys.argv[1], sys.argv[2], sys.argv[3]]
descr = session_resource_path + descr
fasta = session_resource_path + fasta
name_tag = ""
if(len(sys.argv) > 4):
    name_tag = sys.argv[4]
###############################
#  output files
descr_f = open(descr, 'w')
fasta_f = open(fasta, 'w')
forward_fasta = re.sub(".fasta", ".fwd.fasta", fasta)
forward_fasta_f = open(forward_fasta, 'w')
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
#  (together with the  gene/protein/trnascript entry it belongs to)
known_proteome_f = ensembldb + 
orig_search_ids = find_by_blasting("%s%s/%s" % (known_proteome_f, orig_genome, orig_genome), seqfile, working_results_f)
logger("the closest match to %s in  %s is %s\n" % (seqfile, orig_genome, orig_search_ids[0]))

split_info = re.split(" ", orig_search_ids[0])
[orig_protein, orig_gene_name, orig_transcript, orig_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
logger("%s id in %s: %s\n" % (seqfile, orig_genome, orig_protein))

logger("gene name in %s:   %s\n\n" % (orig_genome, orig_gene_name))
orig_seq_f = "orig_seq.fasta"
extract_first_seq(working_results_f, orig_protein, orig_seq_f)
###############################
#  blast the query sequence against all of the remaining proteomes
#
not_found = []

for species in species_list:
    ###forward search###
    logger("%s\n" % species)
    found = 0
    found_in_known = 0;
    found_in_abinitio = 0

    forward_ids = find_by_blasting("%s%s/%s" % (known_proteome_f, species, species), orig_seq_f, working_results_f)
    if (len(forward_ids) == 0):
        logger("%s not found in %s, \"known\" sequences\n\n" % (seqfile, species))
        break
    else:
        # if we found something, we still need to check if it is the mutual best match
        logger("best forward hit:  %s  \n\n" % forward_ids[0])
        split_info = re.split(" ", forward_ids[0])
        [protein, gene_name, transcript, gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
        output_protein = protein; output_gene_location = gene_location
        extract_first_seq (working_results_f, protein, "%s.fasta" % species)
        spec_name = short_name(species)
        #write forward protein to out file
        forward_fasta_f.write(">%s_FWD\n" % spec_name)
        forward_fasta_f.write(popen("grep -v \'>\' %s.fasta" % species).read())
        
        ###back search###
        back_ids = find_by_blasting("%s%s/%s" % (known_proteome_f, orig_genome, orig_genome), "%s.fasta" % species, working_results_f)
        if (len(back_ids) == 0):
            logger("\t reciprocal search in %s using %s as a query produced no hits (?)\n\n" % (orig_genome, species))
            logger("*****************************************\n\n\n")
            continue
        
        # check whether back search retrieves the original query
        # three possibilities here
        # 1) the first hit is the original gene ==> this is the mutual best hit
        # 2) the original gene does not appear on the list -- we take there is no ortohologue in this specie
        # 3) the original gene is somewhere down the list  -- for now take that it also means "orthologue not found"
        split_info = re.split(" ", back_ids[0])
        [back_protein, back_gene_name, back_transcript, back_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
        
        # if the back gene is the same as the original one, we are done
        # otherwise, go and check in the "ab initio" detected list
        if (back_gene_name == orig_gene_name):
            logger("\t best mutual hit:  %s  %s\n\n" % (gene_name, protein))
            found = 1;
            found_in_known = 1;
        else:
            logger("\t back search returns  %s as the best hit\n\n" % back_gene_name)
            logger("\t no mutual best found in known -- try ab initio \n\n")
        
    # We didn't find the protein in the "known" set of genes, or it wasn the mutual best
    # Can we find it in the abinitio annotated species?
    if(found ==  0):
        ab_init_forward_ids = find_by_blasting("%s%s/%s" % (abinit_proteome_f, species, species), orig_seq_f, working_results_f)
        if (len(ab_init_forward_ids) == 0):
            not_found.append(species)
            logger("\t\t %s not found in %s, \"ab initio\" sequences\n\n" % (seqfile, species))
        else:
            split_info = re.split(" ", ab_init_forward_ids[0])
            [protein, gene_name, transcript, gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
            output_protein = protein; output_gene_location = gene_location
            extract_first_seq(working_results_f, protein, "%s.fasta" % species)
            
            #if the back gene is the same as original, we found the match in abinitio
            ab_init_back_ids = find_by_blasting("%s%s/%s" % (known_proteome_f, orig_genome, orig_genome), "%s.fasta" % species, working_results_f)
            split_info = re.split(" ", ab_init_back_ids[0])
            [back_protein, back_gene_name, back_transcript, back_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
        
            if (back_gene_name == orig_gene_name):
                found = 1
                found_in_abinitio = 1
                logger("\t\t best mutual hit:  %s  %s\n\n" % (gene_name, protein))
            else:
                not_found.append(species)
                forward_ids_dict[species] = forward_ids
                ab_init_forward_ids_dict[species] = ab_init_forward_ids
                logger("\t\t no mutual best found in ab initio either\n\n\n")
    logger("*****************************************\n\n")
    if (found == 0):
        system("rm %s.fasta" % species)
        continue
    #Done only for found genes
    sequence = popen("cat %s.fasta" % species).read()
    lines = re.split("\n", sequence)
    spec_name = short_name(species)
    
    if(seen.has_key(spec_name)):
        logger("name %s appears twice in %s\n\n" % (spec_name, species))
        system("rm %s.fasta" % species)
        continue
    else:
        seen[spec_name] = 1
    # determine status of the protein
    if (found_in_known == 1):
        stats = 'KNOWN'

    if (found_in_abinitio == 1):
        stats = 'ABINITIO'
    #write to output files
    descr_f.write("%s\n%s\n" % (species, spec_name))
    descr_f.write("%s %s %s\n\n" % (output_protein, output_gene_location, stats))
    lines[0] = ">" + spec_name + "\n"
    text = ""
    fl = 0
    for line in lines:
        if fl == 0:
            fl = 1
            continue
        text += line
    fasta_f.write(lines[0])
    fasta_f.write(re.sub("\n", "", text))
    fasta_f.write("\n")
    
    system("rm %s.fasta" % species)
    #
    
fasta_f.close()
descr_f.close()
###############################
###############################
#Summary of proteins that were not found
logger("not found:\n")

for item in not_found:
    logger("\t%s\n \tknown\n\n" % item)
    for hit in forward_ids_dict[item]:
        hit = re.sub("\n", "", hit)
        logger("\t\t %s\n" % hit)
    logger("\t\t %s\n" % hit)
    for hit in ab_init_forward_ids_dict[item]:
        hit = re.sub("\n", "", hit)
        logger("\t\t %s\n" % hit)
    logger("\n\n")
    
###############################
###############################
#cleanup
remove(orig_seq_f)
remove(working_results_f)
###############################
###############################
# optional: sort and align
mafft =  "mafft --localpair --maxiterate 1000 "

afafile = fasta;
afafile = re.sub(".fasta", ".mafft.afa", afafile)
cmd = "%s --quiet %s > %s" % (mafft, fasta, afafile)
system(cmd);

afafile = forward_fasta
afafile = re.sub(".fasta", ".mafft.afa", afafile)
cmd = "%s --quiet %s > %s" % (mafft, forward_fasta, afafile)
system(cmd)