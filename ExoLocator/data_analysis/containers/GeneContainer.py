'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.Gene import Gene

@Singleton
class GeneContainer(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._gene_container = {}
    
    def add(self, gene_id, gene):
        if type(gene) is not Gene:
            raise TypeError('Gene, TypeError, {0}'.format(gene_id))
        
        if self._gene_container.has_key(gene_id):
            raise KeyError('Gene, KeyError, {0}'.format(gene_id))
        
        self._gene_container[gene_id] = gene

        
    def get(self, gene_id):
        return self._gene_container[gene_id]