'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.Exons import Exons
from data_analysis.base.EnsemblExons import EnsemblExons

@Singleton
class ExonContainer(object):
    '''
    exon key = (ref_protein_id, species, exon_type)
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._exon_container = {}
        
    def add(self, exon_type, data_map_key, exons):
        exons_key = (data_map_key[0], data_map_key[1], exon_type)
        if type(exons) is not Exons and type(exons) is not EnsemblExons: 
            raise TypeError('Exons, TypeError, {0}, {1}, {2}'.format(exon_type, data_map_key[0], data_map_key[1]))
        
        if self._exon_container.has_key(exons_key):
            raise KeyError('Exons, KeyError, {0}, {1}, {2}'.format(exon_type, data_map_key[0], data_map_key[1]))
        
        self._exon_container[exons_key] = exons

        
    def get(self, exon_key):
        return self._exon_container[exon_key]
    
    def update (self, exon_key, exons):
        self._exon_container[exon_key] = exons