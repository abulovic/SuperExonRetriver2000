'''
Created on May 2, 2012

@author: intern
'''
import os
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from utilities.Logger import Logger
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from data_analysis.base.GenewiseExon import GenewiseExon

class GenewiseExons(object):
    '''
    classdocs
    '''


    def __init__(self, data_map_key, ref_species):
        '''
        Constructor
        '''
        self.ref_protein_id = data_map_key[0]
        self.species        = data_map_key[1]
        self.ref_species    = ref_species
        self.exons          = {}
        
    def load_exons(self):
        
        dc = DirectoryCrawler()
        logger = Logger.Instance()
        container_logger = logger.get_logger('containters')
        
        exon_file_path = dc.get_exon_genewise_path(self.ref_protein_id)
        exon_file_path += "/%s.fa" % self.species
        
        if os.path.isfile(exon_file_path):
            container_logger.error ("{0},{1},genewise,no fasta file for genewise exons.".format(self.ref_protein_id, self.species))
            return False
        
        exon_file = open(exon_file_path, 'r')
        seq_records = SeqIO.parse(exon_file, "fasta", unambiguous_dna)
        
        for seq_record in seq_records:
            (num,ir1,ir2,data) = seq_record.split()
            num = int(num)
            (length, start, stop) = data.split('|')
            
            exon = GenewiseExon((self.ref_protein_id, self.species), num, start, stop, seq_record.seq)
            self.exons[num] = exon
        