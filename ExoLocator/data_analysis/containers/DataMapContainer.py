'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton
from data_analysis.base.DataMap import DataMap

@Singleton
class DataMapContainer(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        key: (protein_id, species)
        '''
        self._data_map_container = {}
        
    def add(self, key, data_map):
        '''
        @param key: (protein_id, species)
        '''
        if type(data_map) is not DataMap:
            raise TypeError
        
        if self._data_map_container.has_key(key):
            raise KeyError
        
        self._data_map_container[key] = data_map

        
    def get(self, key):
        '''
        @param key: (protein_id, species)
        '''
        return self._data_map_container[key]
    
def main():
    dmc = DataMapContainer.Instance()

    dmc.add("id1", "HS", DataMap("id1"))
    #dmc.add("id2", "HS", "DataMap string")  #type err
    #dmc.add("id1", "HS", DataMap("id1"))    # key err
    

if __name__ == '__main__':
    main()