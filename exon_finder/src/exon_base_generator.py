'''
Created on Mar 6, 2012

@author: marioot
'''
import sys
import re
import ConfigParser
from os import system
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
###############################
###############################
#  Tools configuration and setting paths
config_file = "../../marioot.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
# session resource
session_resource_path = config.get('Session resource', 'session_resource_path')

###############################
###############################
#  Use wise to find exons in <dna_f> for <protein_f>
def get_exon_list_wise(protein_f, dna_f):
    # wise config
    wise = config.get('Wise cfg', 'wise')
    wise_flags = config.get('Wise cfg', 'flags')
    wise_out = session_resource_path + config.get('Wise cfg', 'outfile')

    cmd = "%s %s %s %s > %s" % (wise, protein_f, dna_f, wise_flags, wise_out)
    system(cmd)
    return wise_out

###############################
###############################
#  With information gathered from wise
#  create the exon database.
def create_exon_db(wise_f, dna_f):
    # database config
    out_f = session_resource_path + config.get('Exon database path', 'exons') + "exons.fa"
    out_f_individual_exon = session_resource_path + config.get('Exon database path', 'exons') + "individual_exons/exon"
    formatdb = config.get('Formatdb cfg', 'formatdb')
    formatdb_flgs = config.get('Formatdb cfg', 'dna_flgs')
    # Read DNA file
    input_DNA_f = open(dna_f, 'r');
    input_DNA_f.readline();
    DNA = re.sub("\n", "", input_DNA_f.read());
    
    DNA = Seq( DNA, IUPAC.unambiguous_dna )
    input_DNA_f.close()
    # Read wise2_out file
    wise_txt = open(wise_f, 'r');
    # Create exon file
    out = open(out_f, 'w')
    # Fill the exon file
    exon_cnt = 0
    regex = re.compile(r'  Exon (\d+) (\d+) phase \d+')
    for line in wise_txt:
        m = re.match(regex, line)
        if m is not None:
            exon_cnt += 1
            lower = int(m.groups()[0]) - 1
            upper = int(m.groups()[1])
            record = SeqRecord(Seq(str(DNA[lower:upper]),
                       IUPAC.unambiguous_dna),
                       id=str(exon_cnt), description="exon length {0}".format(upper - lower))
            SeqIO.write(record, out, "fasta")
            
            individual_exon_f = open(out_f_individual_exon + "_{0}.fa".format(exon_cnt), 'w')
            SeqIO.write(record, individual_exon_f, "fasta")
            individual_exon_f.close()
    out.close()
    # Format database
    cmd = "%s -i %s %s" % (formatdb, out_f, formatdb_flgs)
    system(cmd)
    return exon_cnt
###############################
###############################
#  Call parameters are protein file and a base species DNA file (usually human)
def main(argv=None):
    if argv is None:
        argv = sys.argv
    protein_f = argv[1]
    dna_f = argv[2]
    
    wise_f = get_exon_list_wise(protein_f, dna_f)
    number_of_exons = create_exon_db(wise_f, dna_f)
    
    return number_of_exons
    
if __name__ == '__main__':
    sys.exit(main())