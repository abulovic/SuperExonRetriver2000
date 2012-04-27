'''
Created on Apr 25, 2012

@author: marioot
'''

class Transcript(object):
    '''
    classdocs
    '''
    def __init__(self, transcript_id, data_map_key, ref_species):
        '''
        Constructor
        '''
        self.transcript_id = transcript_id
        self.ref_protein   = data_map_key[0]
        self.species       = data_map_key[1]
        self.ref_species   = ref_species
        