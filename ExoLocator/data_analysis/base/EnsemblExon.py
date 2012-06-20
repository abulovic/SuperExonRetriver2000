'''
Created on Apr 30, 2012

@author: ana
'''

class EnsemblExon(object):
    '''
    Class representing one ensembl exon.
    '''
    
    def __init__(self, data_map_key, exon_id, start, stop, strand, sequence):
        '''
        @param data_map_key: (ref_protein_id, species)
        '''
        self.ref_protein_id = data_map_key[0]
        self.species        = data_map_key[1]
        self.exon_id        = exon_id
        self.start          = int(start)
        self.stop           = int(stop)
        self.length         = abs(self.stop - self.start)
        self.strand         = int(strand)
        self.sequence       = sequence
        
    def set_exon_ordinal (self, ordinal):
        self.ordinal = ordinal
        
    def set_frame (self, frame):
        self.frame = frame
        
    def set_relative_start(self, start):
        self.relative_start = start
        
    def set_relative_stop (self, stop):
        self.relative_stop = stop
        
    
        