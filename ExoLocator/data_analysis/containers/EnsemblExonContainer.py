'''
Created on Apr 30, 2012

@author: intern
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.EnsemblExon import EnsemblExon

@Singleton
class EnsemblExonContainer(object):
    '''
    Contains instances of EnsemblExon class
    '''


    def __init__(self):
        
        self._ensembl_exon_container = {}
        
    def add(self, exon_id, exon):
        if type(exon) is not EnsemblExon:
            raise TypeError('Gene, TypeError, {0}'.format(exon_id))
        
        if self._ensembl_exon_container.has_key(exon_id):
            raise KeyError('Gene, KeyError, {0}'.format(exon_id))
        
        self._ensembl_exon_container[exon_id] = exon
        
    def get(self, exon_id):
        try:
            return self._ensembl_exon_container[exon_id]
        except KeyError:
            raise KeyError("No ensembl exon id in container: %s" % exon_id)
        