'''
Created on Apr 25, 2012

@author: marioot
'''
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

class Gene(object):
    '''
    classdocs
    '''
    def __init__(self, gene_id, data_map_key, ref_species):
        '''
        Constructor
        '''
        self.gene_id     = gene_id
        self.ref_protein = data_map_key[0]
        self.species     = data_map_key[1]
        self.ref_species = ref_species
        
    def get_expanded_gene_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_expanded_gene_path(self.ref_protein), self.species)
    
    def get_gene_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_gene_path(self.ref_protein), self.species)
        
def main ():
    protein = Gene("ENSAMEG00000021110", ("ENSP00000372410", "Ailuropoda_melanoleuca"), "Homo_sapiens")
    print protein.get_expanded_gene_file_path()
    print protein.get_gene_file_path()
    
if __name__ == '__main__':
    main()