'''
Created on Oct 6, 2011

@author: Mario
'''

from os import system, popen;
import re;
import ConfigParser;

###############################
###############################
#Tools configuration
config_file = "marioot.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
genome_path = config.get('Database path', 'genome_path')
fastacmd = "fastacmd"
blast = "blastall"
blast += " -p blastp -e 1.e-2  -m 8 "; # 8 == output in a tabular format
###############################
###############################
#
def find_by_blasting (proteome, query_sequence, blast, temp_seq_file_name):
    cmd = "%s -d %s -i %s -o tmp_blastout" % (blast, proteome, query_sequence); #| head -n 10 | awk \'{print \$2}\' > tmp_ids ";
    system(cmd);
    
    system("head -n 10 tmp_blastout | awk \'{print $2}\' > tmp_ids");
    
    #system("%s -db %s -query tmp_ids > %s" % (fastacmd, proteome, temp_seq_file_name));
    system("%s -d %s -i tmp_ids > %s" % (fastacmd, proteome, temp_seq_file_name));
  
    headers = popen("grep \'>\' %s" % temp_seq_file_name).read().split('\n');
    
    ids = [];
    
    for hdr in headers:
        #parse header
        hdr1 = re.sub('>', '', hdr);
        hdr1 = re.sub('\|', ':', hdr1);
        
        fields = re.split(' ', hdr1);
        
        parsed = 1;
        
        protein = gene = transcript = location = novelty = "";
        
        for field in fields:
            aux = re.split(':', field);
            
            if aux[0] == "":
                continue;
            elif aux[0] == "lcl":
                protein = aux[1];
            elif aux[0] == "gene":
                gene = aux[1];
            elif aux[0] == "transcript":
                transcript = aux[1];
            elif aux[0] == "pep":
                novelty = aux[1];
            elif location == "":
                location = field;
            else:
                print field;
                parsed = 0
                break;
            
        if parsed == 0:
            print "unrecognized header format:\n %s " % hdr;
            print fields;
            print "%s | %s | %s | %s | %s\n" % (protein, gene, transcript, location, novelty);
            return [];
            
        #print "%s\n   %s %s %s %s %s\n" % (hdr, protein, gene, transcript, location, novelty );
        ids.append("%s %s %s  %s  %s" % (protein, gene, transcript, location, novelty ));
        
    return ids;


def extract_first_seq(infile, name, outfile):

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
    
    outf = open(outfile, 'w');
    outf.write(">%s\n" % (name));
    outf.write(re.sub("\n", "", out));    
    outf.close();
    
    
