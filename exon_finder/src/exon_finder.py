'''
Created on Mar 5, 2012

@author: marioot
'''

from os import system
from os import path
import re;
import sys;
import ConfigParser;

###############################
###############################
#  Tools configuration and setting paths
config_file = "marioot.cfg"
config = ConfigParser.RawConfigParser()
config.read(config_file)

#session resource
session_resource_path = config.get('Session resource', 'session_resource_path')
#resulting proteins path
proteins_path = session_resource_path + config.get('Found proteins path', 'proteins')
regions_path = session_resource_path + config.get('Gene regions path', 'regions')
log_path = session_resource_path + config.get('LOG', 'exon_finder')

blastn = config.get('Blast cfg', 'blastn')
tblastn = config.get('Blast cfg', 'tblastn')
database = "/Users/Mario/Desktop/project/ExonFinder/ExonFinder/src/human_parf_exons.fasta"

###############################
#  input negotiation
#
#  descr = "test.descr";
if(len(sys.argv) < 1):
    print "Usage: %s <descr file>\n" % re.split("/", sys.argv[0])[len(re.split("/", sys.argv[0])) - 1];
    exit;
descr = sys.argv[1] + ".descr"
descr = session_resource_path + descr

