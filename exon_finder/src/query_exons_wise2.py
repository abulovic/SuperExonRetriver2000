'''
Created on Mar 6, 2012

@author: marioot
'''

import re
import sys
import ConfigParser
from os import system

###############################
#  input negotiation
#
#descr = "test.descr";
if(len(sys.argv) < 2):
    print "Usage: %s <protein-file> <dna-file>\n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
protein_f = sys.argv[1]
dna_f = sys.argv[2]
###############################
###############################
#  Tools configuration and setting paths
config_file = "marioot.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)
# session resource
session_resource_path = config.get('Session resource', 'session_resource_path')
# wise config
wise = config.get('Wise cfg', 'wise')
wise_flags = config.get('Wise cfg', 'flags')
wise_out = session_resource_path + config.get('Wise cfg', 'outfile')
###############################
###############################
#  Get exon information in a session_resource/wise2_out
cmd = "%s %s %s %s > %s" % (wise, protein_f, dna_f, wise_flags, wise_out)
print cmd
system(cmd)

# Create Exon database

