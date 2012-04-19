'''
Created on Apr 19, 2012

@author: marioot
'''
from subprocess                          import Popen, PIPE, STDOUT
from utilities.ConfigurationReader       import ConfigurationReader
from utilities.DescriptionParser         import DescriptionParser
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
import re, os
from utilities.Logger import Logger
from pipeline.utilities.CommandGenerator import CommandGenerator

def populate_sequence_exon_genewise(protein_id):
    '''
    Populates the "/PROTEIN_ID/sequence/exon/genewisel/<species>.fa" 
    folder with fasta files containing a list of all the exons for
    a particular transcript. The data is acquired using the genewise
    program.
    This is used for the proteins found with an ab_initio method, that
    dont have a list of exons on ensembl.
    '''
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('data_retrieval')
    directory_crawler   = DirectoryCrawler()
    command_generator   = CommandGenerator()
    
    try:
        (proteins_known, proteins_abinitio) = DescriptionParser().parse_descr_file(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, {1}, , {2}".format(protein_id, 'GENEWISE', e))
        return
    
    failed_species_list = []
    for (species, data) in proteins_abinitio.items():
        protein_file    = "{0}/{1}.fa".format(directory_crawler.get_protein_path(protein_id), species)
        gene_file       = "{0}/{1}.fa".format(directory_crawler.get_gene_path(protein_id), species)
        wise_file       = "{0}/{1}.genewise".format(directory_crawler.get_exon_genewise_path(protein_id), species)
        exon_file       = "{0}/{1}.fa".format(directory_crawler.get_exon_genewise_path(protein_id), species)
        
        wise_command    = command_generator.generate_genewise_command(protein_file, gene_file, wise_file)
        
        Popen(wise_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        #LOGGING
        wise_file_line  = open(wise_file, 'r').readline()
        wise_file.close()
        
        error_pattern = re.compile("Warning Error")
        if re.search(error_pattern, wise_file_line) is None:
            continue
        alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, 'ENSEMBL', species.strip(), wise_file_line.strip()))
        os.remove(wise_file)
        failed_species_list.append(species.strip())
        
        
        
def main ():
    populate_sequence_exon_genewise("ENSP00000311134")

if __name__ == '__main__':
    main()