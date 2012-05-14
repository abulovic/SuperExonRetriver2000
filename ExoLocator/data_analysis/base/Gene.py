'''
Created on Apr 25, 2012

@author: marioot
'''
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from utilities import FileUtilities

class Gene(object):
    '''
    classdocs
    '''
    def __init__(self, gene_id, data_map_key, ref_species):
        '''
        Constructor
        '''
        self.gene_id     = gene_id
        self.ref_protein_id = data_map_key[0]
        self.species     = data_map_key[1]
        self.ref_species = ref_species
        
    def get_expanded_gene_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_expanded_gene_path(self.ref_protein_id), self.species)
    
    def get_gene_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_gene_path(self.ref_protein_id), self.species)
    
    def get_sequence_record (self, gene_type):
        '''
        Retrieves the gene sequence record dependent on the gene type
        @param gene_type: gene type can be normal and extended
        @return sequence of SeqRecord type
        '''
        if gene_type == "normal":
            try:
                return self.sequence
            except AttributeError:
                self.sequence = FileUtilities.load_fasta_single_record(self.get_gene_file_path())
                    
        else:
            try:
                return self.extended_sequence
            except AttributeError:
                self.extended_sequence = FileUtilities.load_fasta_single_record(self.get_expanded_gene_file_path())
        
def main ():
    protein = Gene("ENSAMEG00000021110", ("ENSP00000372410", "Ailuropoda_melanoleuca"), "Homo_sapiens")
    print protein.get_expanded_gene_file_path()
    print protein.get_gene_file_path()
    
if __name__ == '__main__':
    main()