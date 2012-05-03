'''
Created on May 2, 2012

@author: intern
'''

class GenewiseExon(object):
    '''
    classdocs
    '''


    def __init__(self, data_map_key, exon_ordinal, start, stop, sequence):
        '''
        Constructor
        '''
        self.ref_protein_id = data_map_key[0]
        self.species        = data_map_key[1]
        self.ordinal        = exon_ordinal
        self.start          = int(start)
        self.stop           = int(stop)
        self.sequence       = sequence