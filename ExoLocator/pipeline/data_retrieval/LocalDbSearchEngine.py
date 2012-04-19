'''
Created on Apr 15, 2012

@author: intern
'''
from subprocess                                     import Popen, PIPE, STDOUT
from pipeline.utilities.AlignmentCommandGenerator   import AlignmentCommandGenerator
from pipeline.utilities.DirectoryCrawler            import DirectoryCrawler
from utilities.Logger                               import Logger
from utilities.DescriptionParser import DescriptionParser
from utilities.ConfigurationReader import ConfigurationReader
import re, os

def _populate_gene (protein_id, gene_path, region_expansion):
    '''
    Function that populates the respective path with genomic DNA.
    Serves as both regular and expanded search.
    '''
    logger                      = Logger.Instance()
    alignment_logger            = logger.get_logger('data_retrieval')

    alignment_command_generator = AlignmentCommandGenerator()

    logger_name                 = 'GENE'
    if region_expansion > 0:
        logger_name             = 'EXPANDED GENE'

    try:
        (genes_known, genes_abinitio) = DescriptionParser().get_gene_regions(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, {1}, , {2}".format(protein_id, logger_name, e))
        return
    
    failed_species_list = []
    for species, gene_information in genes_known.items():
        (location_type, assembly, location_id, seq_begin, seq_end, strand) = gene_information
        (seq_begin, seq_end) = _expand_region(seq_begin, seq_end, region_expansion)
        gene_fasta      = "{0}/{1}.fa".format(gene_path, species)
        fastacmd        = alignment_command_generator.generate_fastacmd_gene_command(location_id, species, location_type, gene_fasta, 0, strand, seq_begin, seq_end)
        command_return  = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            #Ignoring sequence location.
            invalid_extension_pattern = re.compile("Ignoring sequence location.")
            if re.search(invalid_extension_pattern, output) is not None:
                continue
            alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, logger_name, species.strip(), output.strip()))
            os.remove(gene_fasta)
            failed_species_list.append(species.strip())
             
    for species, gene_information in genes_abinitio.items():
        (location_type, assembly, location_id, seq_begin, seq_end, strand) = gene_information
        (seq_begin, seq_end) = _expand_region(seq_begin, seq_end, region_expansion)
        gene_fasta      = "{0}/{1}.fa".format(gene_path, species)
        fastacmd        = alignment_command_generator.generate_fastacmd_gene_command(location_id, species, location_type, gene_fasta, 0, strand, seq_begin, seq_end)
        command_return  = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, logger_name, species.strip(), output.strip()))
            failed_species_list.append(species)
            
    if failed_species_list:
        _write_failed_species(gene_path, failed_species_list)

def populate_sequence_expanded_gene (protein_id):
    '''
    Populates the "/PROTEIN_ID/sequence/expanded_gene/<species>.fa" 
    folder with fasta files of expanded DNA regions of all the species 
    registered by the Reciprocal Best Search.
    '''
    expanded_gene_path  = DirectoryCrawler().get_expanded_gene_path(protein_id)
    region_expansion    = int(ConfigurationReader.Instance().get_value('local_ensembl', 'expansion'))
    _populate_gene(protein_id, expanded_gene_path, region_expansion)

def populate_sequence_gene(protein_id):
    '''
    Populates the "/PROTEIN_ID/sequence/gene/<species>.fa" folder
    with fasta files containing DNA sequences of all the species 
    registered by the Reciprocal Best Search.
    '''
    gene_path  = DirectoryCrawler().get_gene_path(protein_id)
    region_expansion    = 0
    _populate_gene(protein_id, gene_path, region_expansion)

def populate_sequence_protein (protein_id):
    '''
    Populates the "/PROTEIN_ID/sequence/protein/<species>.fa" 
    folder with fasta files containing protein sequence for
    all the species registered by the Reciprocal Best Search.
    '''
    logger                      = Logger.Instance()
    alignment_logger            = logger.get_logger('data_retrieval')
    
    alignment_command_generator = AlignmentCommandGenerator()
    directory_crawler           = DirectoryCrawler()
    
    try:
        (proteins_known, proteins_abinitio) = DescriptionParser().get_protein_ids(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, PROTEIN, , {2}".format(protein_id, e))
        return
    
    failed_species_list = []
    for species, orthologous_protein_id in proteins_known.items():
        protein_fasta   = "{0}/{1}.fa".format(directory_crawler.get_protein_path(protein_id), species)
        fastacmd        = alignment_command_generator.generate_fastacmd_protein_command(orthologous_protein_id, species, 'all', protein_fasta)
        command_return  = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            alignment_logger.warning("{0}, PROTEIN, {1}, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species)
          
    for species, orthologous_protein_id in proteins_abinitio.items():
        protein_fasta   = "{0}/{1}.fa".format(directory_crawler.get_protein_path(protein_id), species)
        fastacmd        = alignment_command_generator.generate_fastacmd_protein_command(orthologous_protein_id, species, 'abinitio', protein_fasta)
        command_return  = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            alignment_logger.warning("{0}, PROTEIN, {1}, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species)
    if failed_species_list:
        _write_failed_species(directory_crawler.get_protein_path(protein_id), failed_species_list)

def _write_failed_species(path, failed_species_list):
    '''
    Creates or updates the .status file for a respective sequence
    type with the names of species for which sequences haven't been found.
    '''
    status_file = open("{0}/.status".format(path), 'w')
    for species in failed_species_list:
        status_file.write("{0}\n".format(species))
    status_file.close()
   
def _expand_region(seq_begin, seq_end, region_expansion):
    '''
    Helper function that expands the region of the DNA.
    '''
    if region_expansion == 0: 
        return (seq_begin, seq_end)
    seq_begin_exp   = int(seq_begin) - region_expansion
    seq_end_exp     = int(seq_end) + region_expansion
    if seq_begin_exp < 0:
        return (0, seq_end_exp)
    else:
        return (seq_begin_exp, seq_end_exp)
    
def main ():
    populate_sequence_gene("ENSP00000311134")
    populate_sequence_expanded_gene("ENSP00000311134")
    populate_sequence_protein("ENSP00000311134")
    
if __name__ == '__main__':
    main()