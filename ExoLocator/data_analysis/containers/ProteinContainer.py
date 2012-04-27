'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.Protein import Protein

@Singleton
class ProteinContainer(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        key: protein_id
        '''
        self._protein_container = {}
    
    def add(self, protein_id, protein):
        if type(protein) is not Protein:
            raise TypeError('Protein, TypeError, {0}'.format(protein_id))
        
        if self._protein_container.has_key(protein_id):
            raise KeyError('Protein, KeyError, {0}'.format(protein_id))
        
        self._protein_container[protein_id] = protein

        
    def get(self, protein_id):
        return self._protein_container[protein_id]