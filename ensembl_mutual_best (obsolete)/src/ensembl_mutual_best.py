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
def genome2name(genome):
    aux = re.split("_", genome)
    aux[0] = aux[0].upper()
    aux[1] = aux[1].upper()
    if(name_tag == ""):
        return aux[0][0:3]+"_"+aux[1][0:3]
    else:
        return aux[0][0:3]+"_"+aux[1][0:3]+"_"+name_tag

###############################
###############################
#config file example: marioot.cfg
config_file = "marioot.cfg"

config = ConfigParser.RawConfigParser()
config.read(config_file)
genome_path = config.get('Database path', 'genome_path')
abinit_genome_path = config.get('Database path', 'abinit_genome_path')
all_list = "%s/species" % genome_path
bfastacmd = "fastacmd"
blast = "blastall"
blast += " -p blastp -e 1.e-2  -m 8 "; # 8 == output in a tabular format

#orig_genome = "Homo_sapiens"
#seqfile = "parf_human.fasta"
outfile    = "tmp.fasta"
genomes = popen("cat %s" % all_list).read().split("\n")
#genomes = ["Choloepus_hoffmanni"]

################################################
#  input negotiation
#
#fasta = "test.fasta";
#descr = "test.descr";
if(len(sys.argv) < 4):
    print "Usage: %s <input seq file> <orig_genome>  <descr file (output)> <fasta file (output)> [<name tag>] \n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
[seqfile, orig_genome, descr, fasta] = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]

name_tag = ""
if(len(sys.argv) > 5):
    name_tag = sys.argv[5]
################################################
#  output files
#
descr_f = open(descr, 'w')
fasta_f = open(fasta, 'w')
forward_fasta = re.sub(".fasta", ".fwd.fasta", fasta)
forward_fasta_f = open(forward_fasta, 'w')
LOG = open("log", "w")
################################################
#  dictionaries
#
forward_ids_dict = dict()
ab_init_forward_ids_dict = dict()
seen = dict()
################################################
#  other initialization
#
output_gene_location = ""
output_protein = ""
##################################################################
#  find the  query sequence in the "original" genome (together
#  with the  gene/protein/trnascript entry it belongs to)
#
outfile = "tmp.fasta"
orig_search_ids = find_by_blasting("%s/%s/%s" % (genome_path, orig_genome, orig_genome), seqfile, blast, outfile)

msg = "the closest match to %s in  %s is %s\n" % (seqfile, orig_genome, orig_search_ids[0])
print msg
LOG.write(msg)

split_info = re.split(" ", orig_search_ids[0])
[orig_protein, orig_gene_name, orig_transcript, orig_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]

msg = "%s id in %s: %s\n" % (seqfile, orig_genome, orig_protein)
print msg
LOG.write(msg)

msg = "gene name in %s:   %s\n\n" % (orig_genome, orig_gene_name)
print msg
LOG.write(msg)

orig_seq_file = "orig_seq.fasta"
extract_first_seq(outfile, orig_protein, orig_seq_file)
##################################################################
#  blast the query sequence against all of the remaining genomes
#
not_found = []

for genome in genomes:
    #forward search
    msg = "%s\n" % genome
    print msg
    LOG.write(msg)
    
    found = 0
    found_in_known = 0;
    found_in_abinitio = 0
    outfile = "tmp.fasta"
    
    forward_ids = find_by_blasting("%s/%s/%s" % (genome_path, genome, genome), orig_seq_file, blast, outfile)
    if (len(forward_ids) == 0):
        msg = "%s not found in %s, \"known\" sequences\n\n" % (seqfile, genome)
        print msg
        LOG.write(msg)
        break
    else:
        # if we found something, we still need to check if it is mutual best
        msg = "best forward hit:  %s  \n\n" % forward_ids[0]
        print msg
        LOG.write(msg)
        
        split_info = re.split(" ", forward_ids[0])
        [protein, gene_name, transcript, gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
        output_protein = protein; output_gene_location = gene_location
        extract_first_seq (outfile, protein, "%s.fasta" % genome)
        spec_name = genome2name(genome)
        
        forward_fasta_f.write(">%s_FWD\n" % spec_name)
        forward_fasta_f.write(popen("grep -v \'>\' %s.fasta" % genome).read())
        
        #back search
        back_ids = find_by_blasting("%s/%s/%s" % (genome_path, orig_genome, orig_genome), "%s.fasta" % genome, blast, outfile)
        
        if (len(back_ids) == 0):
            msg = "\t reciprocal search in %s using %s as a query produced no hits (?)\n\n" % (orig_genome, genome)
            print msg
            LOG.write(msg)
            msg = "*****************************************\n\n\n"
            print msg
            LOG.write(msg)
            continue;
        
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
            msg = "\t best mutual hit:  %s  %s\n\n" % (gene_name, protein)
            print msg
            LOG.write(msg)
            found = 1;
            found_in_known = 1;
        else:
            msg = "\t back search returns  %s as the best hit\n\n" % back_gene_name
            print msg
            LOG.write(msg)
            msg = "\t no mutual best found in known -- try ab initio \n\n"
            print msg
            LOG.write(msg)
            
    if(found !=  1):
        
    # We didn't find the protein in the "known" set of genes, or it wasn the mutual best
    # Can we find it in the abinitio annotated genome?

        ab_init_forward_ids = find_by_blasting("%s/%s/%s" % (abinit_genome_path, genome, genome), orig_seq_file, blast, outfile)
        if (len(back_ids) == 0):
            found = 0
            not_found.append(genome)
            msg = "\t\t %s not found in %s, \"ab initio\" sequences\n\n" % (seqfile, genome)
            print msg
            LOG.write(msg)
        else:
            split_info = re.split(" ", ab_init_forward_ids[0])
            [protein, gene_name, transcript, gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
            output_protein = protein; output_gene_location = gene_location
            extract_first_seq(outfile, protein, "%s.fasta" % genome)
            
            # now back with this one
            
            ab_init_back_ids = find_by_blasting("%s/%s/%s" % (genome_path, orig_genome, orig_genome), "%s.fasta" % genome, blast, outfile)
            
            
            split_info = re.split(" ", ab_init_back_ids[0])
            [back_protein, back_gene_name, back_transcript, back_gene_location] = [split_info[0], split_info[1], split_info[2], split_info[4]]
        
            if (back_gene_name == orig_gene_name):
                found = 1
                found_in_abinitio = 1
                msg = "\t\t best mutual hit:  %s  %s\n\n" % (gene_name, protein)
                print msg
                LOG.write(msg)
            else:
                found = 0
                not_found.append(genome)
                
                forward_ids_dict[genome] = forward_ids
                ab_init_forward_ids_dict[genome] = ab_init_forward_ids

                msg = "\t\t no mutual best found in ab initio either\n\n\n"
                print msg
                LOG.write(msg)

    msg = "*****************************************\n\n"
    print msg
    LOG.write(msg)

    if (found == 0):
        system("rm %s.fasta" % genome)
        continue
    

    sequence = popen("cat %s.fasta" % genome).read()
    lines = re.split("\n", sequence)
    spec_name = genome2name(genome)
    
    if(seen.has_key(spec_name)):
        LOG.write("name %s appears twice in %s\n\n" % (spec_name, genome))
        system("rm %s.fasta" % genome)
        continue
    else:
        seen[spec_name] = 1
    
    # determine status of the protein
    if (found_in_known == 1):
        stats = 'KNOWN'

    if (found_in_abinitio == 1):
        stats = 'ABINITIO'

    descr_f.write("%s\n%s\n" % (genome, spec_name))
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
    
    system("rm %s.fasta" % genome)


fasta_f.close()
descr_f.close()
    
msg = "not found:\n"
print msg
LOG.write(msg)

for item in not_found:
    msg = "\t%s\n \tknown\n\n" % item
    print msg
    LOG.write(msg)
    for hit in forward_ids_dict[item]:
        hit = re.sub("\n", "", hit)
        msg = "\t\t %s\n" % hit
        print msg
        LOG.write(msg)
    
    msg = "\t\t %s\n" % hit
    print msg
    LOG.write(msg)
    for hit in ab_init_forward_ids_dict[item]:
        hit = re.sub("\n", "", hit)
        msg = "\t\t %s\n" % hit
        print msg
        LOG.write(msg)
    msg =  "\n\n"
    print msg
    LOG.write(msg)
    
###############################
###############################
#cleanup
remove(orig_seq_file)
#remove(outfile)
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
