'''
Created on May 31, 2012

@author: intern
'''
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from utilities import FileUtilities

class TranslatedProtein (object):
    
    def __init__ (self, protein_id, species):
        self.protein_id = protein_id
        self.species = species
        
    def get_protein_file_path(self):
        return "{0}/{1}.fa".format(DirectoryCrawler().get_assembled_protein_path(self.protein_id), self.species)
        
    def get_sequence_record (self):
        '''
        Tries to get the 'sequence' attribute.
        If there is no such attribute, generate it by loading the appropriate fasta file.
        '''
        try:
            return self.sequence
        except AttributeError:
            self.sequence = FileUtilities.load_fasta_single_record(self.get_protein_file_path(), "protein")
            
        return self.sequence