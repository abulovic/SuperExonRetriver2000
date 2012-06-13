'''
Created on May 2, 2012

@author: intern
'''

class GenewiseExon(object):
    '''
    Class representing one exon as predicted by GeneWise
    '''


    def __init__(self, data_map_key, exon_ordinal, start, stop, sequence):
        '''
        @param data_map_key: (ref_protein_id, species)
        '''
        self.ref_protein_id = data_map_key[0]
        self.species        = data_map_key[1]
        self.ordinal        = exon_ordinal
        self.start          = int(start)
        self.stop           = int(stop)
        self.sequence       = sequence