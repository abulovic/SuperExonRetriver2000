'''
Created on Apr 13, 2012

@author: marioot
'''
from pipeline.DirectoryCrawler import DirectoryCrawler
import os

class AlignmentTargetGenerator(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        konstructor
        '''
        self.crawler = DirectoryCrawler()
    
    def get_blastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with blastn for that protein
        '''       
        path = self.crawler.get_blastn_path(protein_id)
        return self.get_species_list(path)
        
    def get_tblastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with tblastn for that protein
        '''       
        path = self.crawler.get_tblastn_path(protein_id)
        return self.get_species_list(path)
        
    def get_SW_gene_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_gene for that protein
        '''       
        path = self.crawler.get_SW_gene_path(protein_id)
        return self.get_species_list(path)
        
    def get_SW_exon_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_exon for that protein
        '''       
        path = self.crawler.get_SW_exon_path(protein_id)
        return self.get_species_list(path)
        
    def get_species_list(self, path):
        '''
        @param path: returns the list of species in the .status file, if the file doesn't exist, it returns the list contained in species.txt
        '''
        if (os.path.isfile("{0}/.status".format(path))):
            return open('{0}/.status'.format(path), 'r').readlines()
        else:
            return open('../species.txt', 'r').readlines()