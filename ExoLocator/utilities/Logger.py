'''
Created on Apr 15, 2012

@author: intern
'''

import logging.config
import sys, re
from utilities.ConfigurationReader import Singleton

@Singleton
class Logger(object):
    '''
    classdocs
    '''


    def __init__(self):
        
        cfg_file = self.get_cfg_dir()+"/logging.cfg"
        logging.config.fileConfig(cfg_file)
        
    def get_logger (self, logger_name):
        return logging.getLogger(logger_name)
        
    def get_cfg_dir(self):
        '''
        Auxiliary function, gets the configuration directory
        '''
        ex_path = sys.path[0]
        m = re.match("(.*ExoLocator).*", ex_path)
        proj_root_dir = m.groups()[0]
        return proj_root_dir + "/cfg/"
        
    
        
if __name__ == '__main__':
    logger = Logger.Instance()
    l = logger.get_logger('mutual_best')
    l.debug('tralala')