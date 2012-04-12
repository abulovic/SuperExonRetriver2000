'''
Created on Apr 12, 2012

@author: intern
'''

class AlignmentCommandGenerator(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
    def generate_fastacmd_command (self, sequence_id, species_name, sequence_type, location_type, output_file_path, strand, sequence_start, sequence_stop):
        pass
    
    def generate_blastn_command (self, database, input_file, output_file):
        pass
    
    def generate_tblastn_command (self, database, input_file, output_file):
        pass
    
    def generate_SW_command (self, query_sequence_file, target_fasta_db_file, output_file, supress_stdout = True):
        pass
    
    def generate_genewise_command (self, protein_file, dna_file, output_file, additional_flags = True):
        pass  
    
    def generate_formatdb_command (self, input_db_file, sequence_type, additional_flags = True):
        pass