'''
Created on Apr 26, 2012

@author: marioot
'''

class DataMap(object):
    '''
    Class that represents data for a particular (protein_id, species) value.
    '''

    def __init__(self, protein_id, transcript_id, gene_id, data_map_key, ab_initio = False):
        '''
        Essentially, contents of the description file 
        All id-s related to the species they belong to!
        '''
        self.protein_id = protein_id
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.ab_initio = ab_initio
        
        self.ref_protein = data_map_key[0]
        self.species     = data_map_key[1]
        
    def get_blastn(self):
        pass
    
    def get_tblastn(self):
        pass
    
    def get_SW_gene(self):
        pass
    
    def get_SW_exon(self):
        pass