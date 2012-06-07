'''
Created on Apr 13, 2012

@author: marioot
'''

import os
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from utilities.DescriptionParser import DescriptionParser

class AlignmentTargetGenerator(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        konstructor
        '''
        self.crawler = DirectoryCrawler()
        self.description_parser = DescriptionParser()
    
    def get_blastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species (RBS) not aligned with blastn for that protein
        '''       
        path = self.crawler.get_blastn_path(protein_id)
        return self._get_species_list(protein_id, path)
    
    def set_failed_blastn_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_blastn_path(protein_id)
        self._write_failed_species_to_status(failed_species_list, path)
        
    def get_tblastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with tblastn for that protein
        '''       
        path = self.crawler.get_tblastn_path(protein_id)
        return self._get_species_list(protein_id, path)
        
    def set_failed_tblastn_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_tblastn_path(protein_id)
        self._write_failed_species_to_status(failed_species_list, path)
        
    def get_SW_gene_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_gene for that protein
        '''       
        path = self.crawler.get_SW_gene_path(protein_id)
        return self._get_species_list(protein_id, path)
    
    def set_failed_SW_gene_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_SW_gene_path(protein_id)
        self._write_failed_species_to_status(failed_species_list, path)
    
    def get_SW_exon_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_exon for that protein
        '''       
        path = self.crawler.get_SW_exon_path(protein_id)
        return self._get_species_list(protein_id, path)
    
    def set_failed_SW_exon_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_SW_exon_path(protein_id)
        self._write_failed_species_to_status(failed_species_list, path)
        
    def _get_species_list(self, protein_id, path):
        '''
        @param path: returns the list of species in the .status file, if the file doesn't exist, it returns the list parsed from protein
                     description file.
        '''
        print path
        if (os.path.isfile("{0}/.status".format(path))):
            species_list = open('{0}/.status'.format(path), 'r').readlines()
            for species in species_list:
                species = species.strip()
        else:
            species_list = self.description_parser.get_species(protein_id)
    
    def _write_failed_species_to_status(self, failed_species_list, path):
        '''
        @param failed_species_list: writes the list of failed species to path/.status file
        @param path: path to the current protein/operation file 
        '''
        status = open('{0}/.status'.format(path), 'w')
        for species in failed_species_list:
            status.write("{0}\n".format(species))
        status.close()
    
def main ():
    atg = AlignmentTargetGenerator()
    print atg.get_blastn_targets("ENSP00000298743")
    print atg.get_tblastn_targets("ENSP00000298743")
    print atg.get_SW_gene_targets("ENSP00000298743")
    print atg.get_SW_exon_targets("ENSP00000298743")
    
if __name__ == '__main__':
    main()