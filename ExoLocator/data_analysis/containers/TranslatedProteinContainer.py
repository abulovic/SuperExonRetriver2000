'''
Created on May 31, 2012

@author: intern
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.TranslatedProtein import TranslatedProtein


@Singleton
class TranslatedProteinContainer(object):

    def __init__(self):
        '''
        Constructor
        '''
        self._translated_protein_container = {}
        
    def add (self, protein_id, species, translated_protein):
        prot_key = (protein_id, species)
        
        if type(translated_protein) is not TranslatedProtein: 
            raise TypeError('TranslatedProteinContainer, TypeError, {0}, {1}'.format(protein_id, species))
        
        if self._translated_protein_container.has_key(prot_key):
            raise KeyError('TranslatedProtein, KeyError, {0}, {1}'.format(protein_id, species))
        
        self._translated_protein_container[prot_key] = translated_protein
        
    def get (self, protein_id, species):
        return self._translated_protein_container[(protein_id, species)]
                           
                        

        