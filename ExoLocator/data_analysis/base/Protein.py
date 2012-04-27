'''
Created on Apr 25, 2012

@author: marioot
'''
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

class Protein(object):
    '''
    classdocs
    '''

    def __init__(self, protein_id, data_map_key, ref_species):
        '''
        Protein class constructor
        
        WARNING: possibility of multiple ref_proteins, this case is not handled!
        '''
        self.protein_id  = protein_id
        self.ref_protein = data_map_key[0]
        self.species     = data_map_key[1]
        self.ref_species = ref_species
        
    def get_protein_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_protein_path(self.ref_protein), self.species)
    
    
def main ():
    protein = Protein("ENSAMEP00000021110", ("ENSP00000372410", "Ailuropoda_melanoleuca"), "Homo_sapiens")
    print protein.get_protein_file_path()
    
if __name__ == '__main__':
    main()