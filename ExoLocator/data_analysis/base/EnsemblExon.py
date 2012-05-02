'''
Created on Apr 30, 2012

@author: intern
'''

class EnsemblExon(object):
    '''
    classdocs
    '''


    def __init__(self, data_map_key, exon_id, start, stop, strand, sequence):
        '''
        Constructor
        '''
        self.ref_protein_id = data_map_key[0]
        self.species = data_map_key[1]
        self.exon_id = exon_id
        self.start = int(start)
        self.stop = int(stop)
        self.length = abs(self.stop - self.start)
        self.strand = int(strand)
        self.sequence = sequence
        
    def set_exon_ordinal (self, ordinal):
        self.ordinal = ordinal
        
    
        