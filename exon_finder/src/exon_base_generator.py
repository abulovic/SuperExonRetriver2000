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
protein_exon_regions_out = ""
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
    exon_outfile = open(exons_out, 'w')
    protein_exon_mapping_file = open(protein_exon_regions_out, 'w')
    # Fill the exon file
    exon_cnt = 0
    regex_exon = re.compile(r'  Exon (\d+) (\d+) phase \d+')
    regex_protein = re.compile(r'     Supporting \d+ \d+ (\d+) (\d+)')
    for line in wise_txt:
        exon_match = re.match(regex_exon, line)
        if exon_match is not None:
            exon_cnt += 1
            lower = int(exon_match.groups()[0]) - 1
            upper = int(exon_match.groups()[1])
            record = SeqRecord(Seq(str(DNA[lower:upper]),
                       IUPAC.unambiguous_dna),
                       id=str(exon_cnt), description="exon length {0}".format(upper - lower))
            SeqIO.write(record, exon_outfile, "fasta")
        protein_match = re.match(regex_protein, line)
        if protein_match is not None:
            lower = int(protein_match.groups()[0])
            upper = int(protein_match.groups()[1])
            protein_exon_mapping_file.write("exon {0} {1} {2}\n".format(str(exon_cnt), lower, upper))
    exon_outfile.close()
    protein_exon_mapping_file.close()
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
    print argv
    [protein_f, dna_f, session_resource_path] = [argv[1], argv[2], argv[3]]
    # session resource
    global wise_out
    global exons_out
    global protein_exon_regions_out
    wise_out = "{0}/{1}".format(session_resource_path, config.get('Wise cfg', 'outfile'))
    exons_out = "{0}/{1}/exons.fa".format(session_resource_path, config.get('Exon database path', 'exons_path'))
    protein_exon_regions_out = "{0}/{1}/protein_regions".format(session_resource_path, config.get('Exon database path', 'exons_path'))
    ###############################
    get_exon_list_wise(protein_f, dna_f)
    number_of_exons = create_exon_db(wise_out, dna_f)
    return number_of_exons
    
if __name__ == '__main__':
    sys.exit(main())