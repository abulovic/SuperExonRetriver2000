'''
Created on Mar 6, 2012

@author: marioot
'''
import sys, re;
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

###############################
#  input negotiation
#
#descr = "test.descr";
if(len(sys.argv) < 2):
    print "Usage: %s <wise2_out> <dna-file>\n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
wise_out = sys.argv[1]
dna_f = sys.argv[2]

input_f = open(dna_f, 'r');
output = input_f.readline();
input = re.sub("\n", "", input_f.read());

my_seq = Seq( input, IUPAC.unambiguous_dna );

input_f.close();
output_f = open(sys.argv[1], 'w');
output_f.write(output);
output_f.write(str(my_seq.reverse_complement()));