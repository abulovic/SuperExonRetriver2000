'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton

@Singleton
class ExonContainer(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._exon_container = {}