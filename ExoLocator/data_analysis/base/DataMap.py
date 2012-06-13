'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import ConfigurationReader

class DataMap(object):
    '''
    Class that represents data for a particular (protein_id, species) value.
    '''

    def __init__(self, protein_id, transcript_id, gene_id, data_map_key, location_type, location_id, strand, start, stop, ab_initio = False):
        '''
        Essentially, contents of the description file 
        All id-s related to the species they belong to!
        '''
        self.protein_id         = protein_id
        self.transcript_id      = transcript_id
        self.gene_id            = gene_id
        self.ab_initio          = ab_initio
        
        self.ref_protein_id     = data_map_key[0]
        self.species            = data_map_key[1]
        
        self.location_type      = location_type
        self.location_id        = location_id
        self.strand             = int(strand)
        self.start              = int (start)
        self.stop               = int(stop)
        
        
    def get_expanded_start(self):
        conf_reader = ConfigurationReader.Instance()
        return max(1, int(self.start - int(conf_reader.get_value('gene_expansion', 'expand'))))
