'''
Created on Oct 18, 2011

@author: Mario
'''


import sys, re;
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


input_f = open(sys.argv[1], 'r');
output = input_f.readline();
input = re.sub("\n", "", input_f.read());

my_seq = Seq( input, IUPAC.unambiguous_dna );

input_f.close();
output_f = open(sys.argv[1], 'w');
output_f.write(output);
output_f.write(str(my_seq.reverse_complement()));
