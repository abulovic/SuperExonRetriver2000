'''
Created on Apr 12, 2012

@author: intern
'''

import os, re
from utilities.ConfigurationReader import ConfigurationReader

class AlignmentCommandGenerator(object):
    '''
    Generates commands for utilities that are used (blast, sw, genewise, fastacmd, formatdb)
    '''

    def __init__(self):
        '''
        Loads the utils configuration (utils.cfg)
        
        if (not os.path.isfile("../utils.cfg")):
            raise IOError("There is no utils.cfg file present in the project directory.")
        
        config = ConfigParser.RawConfigParser()
        config.read("../utils.cfg")
        '''
        
        self.configReader = ConfigurationReader.Instance()
        
        # blast tools
        self.blast_e_value  = self.configReader.get_value('blast', 'expectation')
        
        self.blastn         = self.configReader.get_value('blast', 'blastn')
        self.blastn         = self.blastn % self.blast_e_value
        
        self.tblastn        = self.configReader.get_value('blast', 'tblastn')
        self.tblastn        = self.tblastn % self.blast_e_value
        
        self.blastp         = self.configReader.get_value('blast', 'blastp')
        self.blastp         = self.blastp % self.blast_e_value
        
        # ensembl database
        self.ensembldb      = self.configReader.get_value('local_ensembl', 'ensembldb')
        self.gene_expansion = int (self.configReader.get_value('local_ensembl', 'expansion'))
        self.dna_masked     = int (self.configReader.get_value('local_ensembl', 'masked'))
        
        # Smith-Waterman
        self.sw_sharp       = self.configReader.get_value('sw#', 'sw#')
        
        # genewise
        self.genewise       = self.configReader.get_value('wise', 'wise')
        self.genewise_flags = self.configReader.get_value('wise', 'flags')
        
        self.mafft          = self.configReader.get_value('mafft', 'mafft')
        
        
    def generate_fastacmd_gene_command (self, sequence_id, 
                                   species_name, 
                                   location_type,
                                   output_file_path, 
                                   masked = 0,
                                   strand = None, 
                                   sequence_start = None, sequence_stop = None):
        
        
        
        if (location_type != "chromosome"):
            seq_id_cmd = "-s %s" % sequence_id
        else:
            seq_id_cmd = "-s chrom%s" % sequence_id
        

        data_type_cmd = "-p F"
        database = "-d %s" % self._generate_genedb_file_name(species_name, location_type, sequence_id, masked)
            
        if (strand == None or int(strand) == 1):
            strand_cmd = "-S 1"
        else:
            strand_cmd = "-S 2"
            
        if (sequence_start and sequence_stop):
            location_cmd = "-L %s,%s" % (sequence_start, sequence_stop)
        else:
            location_cmd = ""
            
        output_cmd = "-o %s" % output_file_path
        
        return "fastacmd {0} {1} {2} {3} {4} {5}".format(database, seq_id_cmd, data_type_cmd, strand_cmd, location_cmd, output_cmd)
        
    
    def generate_fastacmd_protein_command (self, protein_id, species_name, protein_type, output_file_path):
        
        data_type_cmd = "-p T"
        prot_id_cmd = "-s %s" % protein_id
        database = "-d %s" % self._generate_proteindb_file_name(species_name, protein_type)
        output_cmd = "-o %s" % output_file_path
        return "fastacmd {0} {1} {2} {3}".format(prot_id_cmd, data_type_cmd, database, output_cmd)
    
    def generate_blastn_command (self, database, input_file, output_file):
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.blastn, database, input_file, output_file)
        return cmd
    
    def generate_tblastn_command (self, database, input_file, output_file):
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.tblastn, database, input_file, output_file)
        return cmd
    
    def generate_blastp_command (self, database, input_file, output_file):
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.blastp, database, input_file, output_file)
        return cmd
    def generate_blastp_command_for_species(self, species_name, input_file, output_file, protein_type):
        database = self._generate_proteindb_file_name(species_name, protein_type)
        return self.generate_blastp_command(database, input_file, output_file)
    
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
            cmd = "formatdb -i {0} -p F".format(input_db_file)
        if additional_flags:
            cmd += " -a F -o F"
            
        return cmd
    
    def generate_mafft_command (self, input_file, output_file):
        return "{0} {1} > {2}".format(self.mafft, input_file, output_file)
    
    
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
