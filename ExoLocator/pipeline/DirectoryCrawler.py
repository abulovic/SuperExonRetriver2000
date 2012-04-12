'''
Created on Apr 12, 2012

@author: intern
'''

import ConfigParser
import os

class DirectoryCrawler(object):
    '''
    classdocs
    '''


    def __init__(self, protein_id = None):
        '''
        Loads the configuration
        '''
        if (not os.path.isfile("../pipeline.cfg")):
            raise IOError("There is no pipeline.cfg file present in the project directory.")
        
        config = ConfigParser.RawConfigParser()
        config.read("../pipeline.cfg")
        
        # absolute project and session directories
        self.project_root_dir   = config.get('root', 'project_dir')
        self.sessions_dir       = config.get('root', 'session_dir')
        
        # sequence databases
        self.sequence_root  = config.get('sequence', 'root')
        self.gene           = "%s/%s" % (self.sequence_root, config.get('sequence', 'gene'))
        self.expanded_gene  = "%s/%s" % (self.sequence_root, config.get('sequence', 'exp_gene'))
        self.protein        = "%s/%s" % (self.sequence_root, config.get('sequence', 'protein'))
        self.exon_ensembl   = "%s/%s" % (self.sequence_root, config.get('sequence', 'exon_ens'))
        self.exon_genewise  = "%s/%s" % (self.sequence_root, config.get('sequence', 'exon_wise'))
        
        # alignment outputs
        self.alignment_root = config.get('alignment', 'root')
        self.blastn_output  = "%s/%s" % (self.alignment_root, config.get('alignment', 'blastn'))
        self.tblastn_output = "%s/%s" % (self.alignment_root, config.get('alignment', 'tblastn'))
        self.SW_gene        = "%s/%s" % (self.alignment_root, config.get('alignment', 'SW_gene'))
        self.SW_exon        = "%s/%s" % (self.alignment_root, config.get('alignment', 'SW_exon'))
        
        # annotation output
        self.annotation_root = config.get('annotation', 'root')
        self.genewise_output = "%s/%s" % (self.annotation_root, config.get('annotation', 'wise'))
        
        # log directory
        self.log_root = config.get('log', 'root')
        self.mutual_best_log = "%s/%s" % (self.log_root, config.get('log', 'mutual_best'))
        
        
    
    def set_protein_id (self, protein_id):
        self.protein_id = protein_id
        
        # sequence absolute paths
        self.sequence_root_abs  = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.sequence_root)
        self.gene_abs           = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.gene)
        self.expanded_gene_abs  = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.expanded_gene)
        self.protein_abs        = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.protein)
        self.exon_ensembl_abs   = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.exon_ensembl)
        self.exon_genewise_abs  = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.exon_genewise)
        
        # alignment outputs absolute paths
        self.alignment_root_abs = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.alignment_root)
        self.blastn_output_abs  = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.blastn_output)
        self.tblastn_output_abs = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.tblastn_output)
        self.SW_gene_abs        = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.SW_gene)
        self.SW_exon_abs        = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.SW_exon)
        
        # annotation output absolute path
        self.annotation_root_abs = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.annotation_root)
        self.genewise_output_abs = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.genewise_output)
        
        # log directory absolute path
        self.log_root_abs           = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.log_root)
        self.mutual_best_log_abs    = "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.mutual_best_log)
        
        
    def get_gene_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.gene)
        else:
            return self.gene_abs
        
    def get_expanded_gene_path(self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.expanded_gene)
        else:
            return self.expanded_gene_abs
    
    def get_protein_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.protein)
        else:
            return self.protein_abs
    
    def get_exon_ensembl_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.exon_ensembl)
        else:
            return self.exon_ensembl_abs
    
    def get_exon_genewise_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.exon_genewise)
        else:
            return self.exon_genewise_abs
    
    def get_blastn_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.blastn_output)
        else:
            return self.blastn_output_abs
    
    def get_tblastn_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.tblastn_output)
        else:
            return self.tblastn_output_abs
    
    def get_SW_gene_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.SW_gene)
        else:
            return self.SW_gene_abs
    
    def get_SW_exon_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.SW_exon)
        else:
            return self.SW_exon_abs
    
    def get_genewise_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}/{1}/{2}".format(self.sessions_dir, protein_id, self.genewise_output)
        else:
            return self.genewise_output_abs
    
    def get_misc_path (self, protein_id = None):
        pass
    
    def get_protein_description_file_path (self, protein_id = None):
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            return "{0}.descr".format(protein_id)
        else:
            return "{0}.descr".format(self.protein_id)
    
    def generate_directory_tree(self, protein_id = None):
        '''
        Generates directory tree for specified protein id
        @param protein_id: If specified, directory tree is generated for that protein.
                           If not, it is expected that the class attribute protein_id is set.
                           If not one of them is set, AttributeException is raised.
        '''
        if (self._generate_absolute_path(self.protein_id, protein_id)):
            working_id = protein_id
        else:
            working_id = self.protein_id
        
        tmp_dir = self.get_gene_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
        
        tmp_dir = self.get_expanded_gene_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_protein_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_exon_ensembl_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_exon_genewise_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_blastn_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_tblastn_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_SW_gene_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_SW_exon_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
            
        tmp_dir = self.get_genewise_path(working_id) 
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
        
    def _generate_absolute_path (self, p_self, p_extern):
        '''
        If p_extern is set, then new absolute path need be generated.
        @return: If p_extern is set, return True, otherwise False 
        '''
        if (not p_self and not p_extern):
            raise AttributeError('DirectoryCrawler instance has not attribute protein_id, cannot generate appropriate protein directory tree')
        if p_extern:
            return True
        else:
            return False
        
    
def main ():
    dc = DirectoryCrawler()
    dc.set_protein_id("ENSP00000311134") 
    print dc.protein_id
    dc.generate_directory_tree()
    
if __name__ == '__main__':
    main()