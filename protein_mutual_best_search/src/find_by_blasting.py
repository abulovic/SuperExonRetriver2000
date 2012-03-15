'''
Created on Oct 6, 2011

@author: Mario
'''

from os import system, popen, remove;
import re;
import ConfigParser;

###############################
###############################
#Tools configuration
config_file = "../../config.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
blastp = config.get('Blast cfg', 'blastp')
fastacmd = config.get('Fastacmd cfg', 'fastacmd')
###############################
###############################
#Finds 10 best matches using BLAST
#in: proteome db path, query sequence path, results path
#out: "proteinID, geneID, transcriptID, location, novelty" list for the results
def find_by_blasting (proteome_db_f, query_sequence_f, result_seq_f):
    cmd = "%s -d %s -i %s -o tmp_blastout" % (blastp, proteome_db_f, query_sequence_f)
    system(cmd)
    tmp_blastout_line = open('tmp_blastout').read()
    if (len(tmp_blastout_line) == 0):
        return []
    
    #get protein IDs from 10 best hits
    system("head -n 10 tmp_blastout | awk \'{print $2}\' > tmp_prot_ids")
    #fastacmd to extract seqs
    print "%s -d %s -i tmp_prot_ids > %s" % (fastacmd, proteome_db_f, result_seq_f)
    system("%s -d %s -i tmp_prot_ids > %s" % (fastacmd, proteome_db_f, result_seq_f))
  
    headers = popen("grep \'>\' %s" % result_seq_f).read().split('\n')
    ids = [];
    for header in headers:
        header_str = parse_header(header) 
        if header_str == "Unrecognized header":
            print "Unrecognized header format:\n %s " % header
            return [];      
        if (header_str != "      "):
            ids.append(header_str); 
    remove('tmp_prot_ids')
    remove('tmp_blastout')
    return ids;
###############################
###############################
#Parses header
def parse_header(header):
    hdr = re.sub('>', '', header)
    hdr = re.sub('\|', ':', hdr)
    fields = hdr.split()
    protein = gene = transcript = location = novelty = ""
    for field in fields:
        aux = re.split(':', field)
        if aux[0] == "":
            continue
        elif aux[0] == "lcl":
            protein = aux[1]
        elif aux[0] == "gene":
            gene = aux[1]
        elif aux[0] == "transcript":
            transcript = aux[1]
        elif aux[0] == "pep":
            novelty = aux[1]
        elif location == "":
            location = field
        else:
            #print field;
            #return "Unrecognized header"
            #print "Unrecognized header elements: %s" % field
            continue
    return "%s %s %s  %s  %s" % (protein, gene, transcript, location, novelty )
###############################
###############################
#Returns the first sequence from the file containing multiple sequences
def extract_first_seq(infile, name, working_results_f):
    inf = open(infile, 'r');
    # would it be better to use the longest transcript here, instead
    # of the first hit?
    out = "";
    reading = False;
    for line in inf.readlines():
        if(re.match('\S', line)):
                if(re.match('^>', line)):
                    if(reading) :
                        break;
                    reading = True;
                else:
                    out += line;
    
    outf = open(working_results_f, 'w');
    outf.write(">%s\n" % (name));
    outf.write(re.sub("\n", "", out));    
    outf.close();