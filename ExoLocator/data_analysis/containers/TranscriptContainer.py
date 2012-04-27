'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.Transcript import Transcript

@Singleton
class TranscriptContainer(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._transcript_container = {}
    
    def add(self, transcript_id, transcript):
        if type(transcript) is not Transcript:
            raise TypeError
        
        if self._transcript_container.has_key(transcript_id):
            raise KeyError
        
        self._transcript_container[transcript_id] = transcript

        
    def get(self, transcript_id):
        return self._transcript_container[transcript_id]