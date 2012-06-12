'''
Created on Apr 15, 2012

@author: intern
'''
# Python imports
import os, sys, re

# utilities imports
import ConfigParser 
from utilities.Singleton import Singleton

       
@Singleton
class ConfigurationReader:
    '''
    Loads configuration files from the cfg directory in the Exolocator source directory
    '''

    def __init__ (self):
        
        self.cfg_dir = self.get_cfg_dir()
        
        self.cfg_files = []
        self.cfg_params = {}
        self.checkout_cfg_files()
        
    def get_cfg_dir(self):
        '''
        Find the absolute path to the configuration directory
        '''
        ex_path = sys.path[0]
        m = re.match("(.*ExoLocator).*", ex_path)
        proj_root_dir = m.groups()[0]
        return proj_root_dir + "/cfg/"
        

    def load_cfg_file(self, cfg_file_path):
        '''
        Load a configuration file and add the
        key, value pairs to the internal configuration dictionary
        '''
        cf = ConfigParser.ConfigParser()
        cf.read(cfg_file_path)
        sections = cf.sections()
        for section in sections:
            for cfg_item in cf.items(section):
                (key, value) = cfg_item
                self.cfg_params.setdefault(section, {}).update({key:value})
        
        
    
    def checkout_cfg_files (self):
        '''
        Check if all the configuration files have been loaded
        '''
        for cfg_file in os.listdir(self.cfg_dir):
            if cfg_file not in self.cfg_files and cfg_file != "logging.cfg" and cfg_file.endswith(".cfg"):
                self.cfg_files.append(cfg_file)
                self.load_cfg_file(self.cfg_dir + cfg_file)
                
    def get_value (self, cfg_section, cfg_item):
        '''
        Try to get the configuration value. If value is not mapped in the internal
        configuration dictionary, try to re-load the configuration directory
        '''
        if cfg_section in self.cfg_params:
            section = self.cfg_params[cfg_section]
            if cfg_item in section:
                return section[cfg_item]
            else:
                self.checkout_cfg_files()
        else:
            self.checkout_cfg_files()
            
        if cfg_section not in self.cfg_params or cfg_item not in self.cfg_params[cfg_section]:
            raise KeyError ('Configuration dictionary has no value mapped to the section %s, item %s' % (cfg_section, cfg_item))
        
        
        
        
if __name__ == '__main__':
    cr = ConfigurationReader.Instance()
    print cr.get_value('wise', 'flags')
    print cr.get_value('tralala', 'lala')