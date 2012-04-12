'''
Created on Apr 12, 2012

@author: intern
'''

class DirectoryCrawler(object):
    '''
    classdocs
    '''


    def __init__(self, protein_id = None):
        '''
        Constructor
        '''
        
    def set_protein_id (self, protein_id):
        self.protein_id = protein_id
        
    def get_gene_path (self, protein_id = None):
        pass
        
    def get_expanded_gene_path(self, protein_id = None):
        pass
    
    def get_protein_path (self, protein_id = None):
        pass
    
    def get_exon_ensembl_path (self, protein_id = None):
        pass
    
    def get_exon_genewise_path (self, protein_id = None):
        pass
    
    def get_blastn_path (self, protein_id = None):
        pass
    
    def get_tblastn_path (self, protein_id = None):
        pass
    
    def get_SW_gene_path (self, protein_id = None):
        pass
    
    def get_SW_exon_path (self, protein_id = None):
        pass
    
    def get_genewise_path (self, protein_id = None):
        pass
    
    def get_misc_path (self, protein_id = None):
        pass
    
    def get_protein_description_file_path (self, protein_id = None):
        pass