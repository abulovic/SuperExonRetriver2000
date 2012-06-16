'''
Created on Apr 13, 2012

@author: marioot
'''

# utilities imports
from utilities.DescriptionParser    import DescriptionParser
from utilities.DirectoryCrawler     import DirectoryCrawler
from utilities.FileUtilities import write_failed_species_to_status,\
    get_species_list


class AlignmentTargetGenerator(object):
    '''
    Class used to retrieve the list of all possible targets from certain alignment
    '''
    def __init__(self):
        self.crawler = DirectoryCrawler()
        self.description_parser = DescriptionParser()
    
    def get_blastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species (RBS) not aligned with blastn for that protein
        '''       
        path = self.crawler.get_blastn_path(protein_id)
        return get_species_list(protein_id, path)
    
    
    def set_failed_blastn_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_blastn_path(protein_id)
        write_failed_species_to_status(failed_species_list, path)
        
    def get_tblastn_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with tblastn for that protein
        '''       
        path = self.crawler.get_tblastn_path(protein_id)
        return get_species_list(protein_id, path)
        
    def set_failed_tblastn_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_tblastn_path(protein_id)
        write_failed_species_to_status(failed_species_list, path)
        
    def get_SW_gene_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_gene for that protein
        '''       
        path = self.crawler.get_SW_gene_path(protein_id)
        return get_species_list(protein_id, path)
    
    def set_failed_SW_gene_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_SW_gene_path(protein_id)
        write_failed_species_to_status(failed_species_list, path)
    
    def get_SW_exon_targets(self, protein_id):
        '''
        @param protein_id: retrieves the list of species not aligned with SW_exon for that protein
        '''       
        path = self.crawler.get_SW_exon_path(protein_id)
        return get_species_list(protein_id, path)
    
    def set_failed_SW_exon_targets(self, protein_id, failed_species_list):
        path = self.crawler.get_SW_exon_path(protein_id)
        write_failed_species_to_status(failed_species_list, path)
        
        
    
    
def main ():
    atg = AlignmentTargetGenerator()
    print atg.get_blastn_targets("ENSP00000298743")
    print atg.get_tblastn_targets("ENSP00000298743")
    print atg.get_SW_gene_targets("ENSP00000298743")
    print atg.get_SW_exon_targets("ENSP00000298743")
    
if __name__ == '__main__':
    main()