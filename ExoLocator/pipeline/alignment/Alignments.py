'''
Created on Apr 13, 2012

@author: marioot
'''

from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from pipeline.alignment.AlignmentCommandGenerator import AlignmentCommandGenerator
from pipeline.alignment.AlignmentTargetGenerator  import AlignmentTargetGenerator
import subprocess
import os

def generate_blastn_alignments(protein_id, species_list = None):
    '''
        Runs the blastn program for a specified protein and list of species
        @param protein_id
        @param species_list: if provided, runs blastn for this list of species, \
                             otherwise runs for species that are missing the blastn output
    '''
    if (not species_list):
        species_list    = AlignmentTargetGenerator().get_blastn_targets(protein_id)
    
    crawler             = DirectoryCrawler()
    command_generator   = AlignmentCommandGenerator()
    for species in species_list:
        output_file     = "{0}/{1}.blastout".format(crawler.get_blastn_path(protein_id), species.strip())
        input_file      = "{0}/{1}.fa".format(crawler.get_expanded_gene_path(protein_id), species.strip())
        database        = "{0}/Homo_sapiens.fa".format(crawler.get_exon_ensembl_path(protein_id))
        
        #subprocess.call(command_generator.generate_blastn_command(database, input_file, output_file))
         
def generate_tblastn_alignments(protein_id, species_list = None):
    '''
        Runs the tblastn program for a specified protein and list of species
        @param protein_id
        @param species_list: if provided, runs tblastn for this list of species, \
                             otherwise runs for species that are missing the tblastn output
    '''
    if (not species_list):
        species_list    = AlignmentTargetGenerator().get_tblastn_targets(protein_id)
    
    crawler             = DirectoryCrawler()
    command_generator   = AlignmentCommandGenerator()
    for species in species_list:
        output_file     = "{0}/{1}.blastout".format(crawler.get_tblastn_path(protein_id), species.strip())
        input_file      = "{0}/{1}.fa".format(crawler.get_protein_path(protein_id), species.strip())
        database        = "{0}/Homo_sapiens.fa".format(crawler.get_exon_ensembl_path(protein_id))
        
        #subprocess.call(command_generator.generate_tblastn_command(database, input_file, output_file))
        
def generate_SW_gene_alignments(protein_id, species_list = None):
    '''
        Runs the SW program for a specified protein and list of species, using the expanded gene region.
        @param protein_id
        @param species_list: if provided, runs SW for this list of species, \
                             otherwise runs for species that are missing the SW output
    '''
    if (not species_list):
        species_list         = AlignmentTargetGenerator().get_SW_gene_targets(protein_id)
    crawler                  = DirectoryCrawler()
    command_generator        = AlignmentCommandGenerator()
    for species in species_list:
        output_file          = "{0}/{1}.swout".format(crawler.get_SW_gene_path(protein_id), species.strip())
        query_sequence_file  = "{0}/{1}.fa".format(crawler.get_expanded_gene_path(protein_id), species.strip())
        target_fasta_db_file = "{0}/Homo_sapiens.fa".format(crawler.get_exon_ensembl_path(protein_id))
        
        #subprocess.call(command_generator.generate_SW_command(query_sequence_file, target_fasta_db_file, output_file))
            
    os.remove(".sw_stdout_supressed")
        
def generate_SW_exon_alignments(protein_id, species_list = None):
    '''
        Runs the SW program for a specified protein, species and individual exons retrived from ensembl.
        @param protein_id
        @param species_list: if provided, runs SW for this list of species, \
                             otherwise runs for species that are missing the SW output
    '''
    if (not species_list):
        species_list    = AlignmentTargetGenerator().get_SW_exon_targets(protein_id)
    
    crawler             = DirectoryCrawler()
    command_generator   = AlignmentCommandGenerator()
    for species in species_list:
        ensembl_exons_path  = "{0}/{1}.fa".format(crawler.get_exon_ensembl_path(protein_id), species.strip())
        target_fasta_db_file = "{0}/Homo_sapiens.fa".format(crawler.get_exon_ensembl_path(protein_id))
        # Isolate each exon to a temporary file to be processed with SW
        tmp_file = "tmp_sw.fa"
        ensembl_exons_file = open(ensembl_exons_path, 'r')
        exon_seq = ""
        exon_header = ""
        exon_counter = 0
        for line in ensembl_exons_file.readlines():
            if line.startswith('>'):
                if exon_seq is not "":
                    exon_counter += 1
                    output_file = "{0}/{1}/exon{2}.swout".format(crawler.get_SW_exon_path(protein_id), species.strip(), exon_counter)
                    
                    tmp = open(tmp_file, 'w')
                    tmp.write("{0}\n{1}\n".format(exon_header, exon_seq))
                    tmp.close()
                    
                    #subprocess.call(command_generator.generate_SW_command(tmp_file, target_fasta_db_file, output_file))
                    print command_generator.generate_SW_command(tmp_file, target_fasta_db_file, output_file)
                exon_header = line.strip()
                exon_seq = ""
            else:
                exon_seq += line.strip()
        
        exon_counter += 1
        output_file = "{0}/{1}/exon{2}.swout".format(crawler.get_SW_exon_path(protein_id), species.strip(), exon_counter)
        
        tmp = open(tmp_file, 'w')
        tmp.write("{0}\n{1}\n".format(exon_header, exon_seq))
        tmp.close()
        
        #subprocess.call(command_generator.generate_SW_command(tmp_file, target_fasta_db_file, output_file))
        print command_generator.generate_SW_command(tmp_file, target_fasta_db_file, output_file)
        ensembl_exons_file.close()
        os.remove(tmp_file)
        os.remove(".sw_stdout_supressed")
  
def main ():
    generate_SW_exon_alignments("ENSP00000311134")
    
if __name__ == '__main__':
    main()