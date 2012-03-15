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
config_file = "../../config.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
wise_out = ""
exons_out = ""
###############################
###############################
#  Use wise to find exons in <dna_f> for <protein_f>
def get_exon_list_wise(protein_f, dna_f):
    # wise config
    wise = config.get('Wise cfg', 'wise')
    wise_flags = config.get('Wise cfg', 'flags')

    cmd = "{0} {1} {2} {3} > {4}".format(wise, protein_f, dna_f, wise_flags, wise_out)
    system(cmd)
    return wise_out

###############################
###############################
#  With information gathered from wise
#  create the exon database.
def create_exon_db(wise_f, dna_f):
    # database config
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
    out = open(exons_out, 'w')
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
    out.close()
    # Format database
    cmd = "%s -i %s %s" % (formatdb, exons_out, formatdb_flgs)
    system(cmd)
    return exon_cnt
###############################
###############################
#  Call parameters are protein file, a base species DNA file (usually human), session resource path
def main(argv=None):
    if argv is None:
        argv = sys.argv
    [protein_f, dna_f, session_resource_path] = [argv[1], argv[2], argv[3]]
    # session resource
    global wise_out
    global exons_out 
    wise_out = "{0}/{1}".format(session_resource_path, config.get('Wise cfg', 'outfile'))
    exons_out = "{0}/{1}/exons.fa".format(session_resource_path, config.get('Exon database path', 'exons_path'))
    ###############################
    wise_f = get_exon_list_wise(protein_f, dna_f)
    number_of_exons = create_exon_db(wise_f, dna_f)
    return number_of_exons
    
if __name__ == '__main__':
    sys.exit(main())