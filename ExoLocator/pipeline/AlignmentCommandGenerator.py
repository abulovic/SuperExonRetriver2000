'''
Created on Apr 12, 2012

@author: intern
'''

import ConfigParser
import os, re

class AlignmentCommandGenerator(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Loads the utils configuration (utils.cfg)
        '''
        if (not os.path.isfile("../utils.cfg")):
            raise IOError("There is no utils.cfg file present in the project directory.")
        
        config = ConfigParser.RawConfigParser()
        config.read("../utils.cfg")
        
        # blast tools
        self.blast_e_value  = config.get('blast', 'expectation')
        
        self.blastn         = config.get('blast', 'blastn')
        self.blastn         = self.blastn % self.blast_e_value
        
        self.tblastn        = config.get('blast', 'tblastn')
        self.tblastn        = self.tblastn % self.blast_e_value
        
        self.blastp         = config.get('blast', 'blastp')
        self.blastp         = self.blastp % self.blast_e_value
        
        # ensembl database
        self.ensembldb      = config.get('local_ensembl', 'ensembldb')
        self.gene_expansion = int (config.get('local_ensembl', 'expansion'))
        self.dna_masked     = int (config.get('local_ensembl', 'masked'))
        
        # Smith-Waterman
        self.sw_sharp       = config.get('sw#', 'sw#')
        
        # genewise
        self.genewise       = config.get('wise', 'wise')
        self.genewise_flags = config.get('wise', 'flags')
        
        
    def generate_fastacmd_command (self, sequence_id, 
                                   species_name, 
                                   sequence_type, 
                                   location_type, 
                                   output_file_path, 
                                   masked,
                                   strand = None, 
                                   sequence_start = None, sequence_stop = None):
        
        database = "-d %s" % self._generate_genedb_file_name(species_name, sequence_type, sequence_id, masked)
        
        seq_id_cmd = "-s %s" % sequence_id
        
        if (sequence_type == "protein"):
            data_type_cmd = "-p T"
        else:
            data_type_cmd = "-p F"
            
        if (strand == None or int(strand) == 1):
            strand_cmd = "-S 1"
        else:
            strand_cmd = "-S 2"
            
        if (sequence_start and sequence_stop):
            location_cmd = "-L %d,%d" % (sequence_start, sequence_stop)
        else:
            location_cmd = ""
            
        output_cmd = "-o %s" % output_file_path
        
        return "fastacmd {0} {1} {2} {3} {4} {5}".format(database, seq_id_cmd, data_type_cmd, strand_cmd, location_cmd, output_cmd)
        
    
    
    
    def generate_blastn_command (self, database, input_file, output_file):
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.blastn, database, input_file, output_file)
        return cmd
    
    def generate_tblastn_command (self, database, input_file, output_file):
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.tblastn, database, input_file, output_file)
        return cmd
    
    def generate_SW_command (self, query_sequence_file, target_fasta_db_file, output_file, supress_stdout = True):
        cmd = "{0} -i {1} -j {2} --out {3}".format(self.sw_sharp, query_sequence_file, target_fasta_db_file, output_file)
        if supress_stdout:
            cmd += " > .sw_stdout_supressed"
        return cmd
    
    def generate_genewise_command (self, protein_file, dna_file, output_file, additional_flags = True):
        pass
    
    def generate_formatdb_command (self, input_db_file, sequence_type, additional_flags = True):
        
        if sequence_type == "protein" or sequence_type == "P":
            cmd = "formatdb -i {0} -p T".format(input_db_file)
        else:
            cmd = "formatdb -i {0} -p T".format(input_db_file)
        if additional_flags:
            cmd += " -a F -o F"
            
        return cmd
    
    
    def _generate_genedb_file_name (self, species, sequence_type, sequence_id, masked):
        '''
        @param species: species name
        @param sequence_type: scaffold / chromosome...
        @param sequence_id: ensembl sequence ID
        @param masked: 0 if dna should not be masked, 1 if it should
        '''
        file_name = "{0}/{1}/dna/".format(self.ensembldb, species.lower())
        # get the template name (dependent on the assembly)
        tmp_file=""
        for f in os.listdir(file_name):
            if (f != "README"):
                tmp_file = f
                break
        m = re.findall ('(.*).dna', tmp_file)   
        if (masked != 0):
            file_name = "%s/%s.dna_rm." % (file_name, m[0])
        else :
            file_name = "%s/%s.dna." % (file_name, m[0])
        if (sequence_type == 'chromosome'):
            file_name = "%schromosome.%s.fa" % (file_name, sequence_id)
        else :
            file_name = "%stoplevel.fa" % (file_name)
        return file_name
       
       
    def _generate_proteindb_file_name (self, species, protein_type):
        '''
        @param species: species name (ensembl)
        @param protein_type: all / abinitio
        @return: protein database name
        '''
        file_name = "%s/%s/pep" % (self.ensembldb, species.lower())
        tmp_file=""
        for f in os.listdir(file_name):
            if (f != "README"):
                tmp_file = f 
                break
        m = re.findall ('(.*).pep', tmp_file)
        if (protein_type == "all"):
            file_name = "%s/%s.pep.all.fa" % (file_name, m[0])
        else:
            file_name = "%s/%s.pep.abinitio.fa" % (file_name, m[0])
            
        return file_name
    
def main():
    acg = AlignmentCommandGenerator()
    cmd = acg.generate_SW_command("query.fa", "target.fa", "output", True)
    print cmd
    
    
if __name__ == '__main__':
    main()
