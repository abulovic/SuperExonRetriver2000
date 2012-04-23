'''
Created on Apr 19, 2012

@author: marioot
'''
from subprocess                          import Popen, PIPE, STDOUT
from utilities.DescriptionParser         import DescriptionParser
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
import re, os, time
from utilities.Logger import Logger
from pipeline.utilities.CommandGenerator import CommandGenerator

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
    exon_genewise_path  = directory_crawler.get_exon_genewise_path(protein_id)
    try:
        (proteins_known, proteins_abinitio) = DescriptionParser().parse_descr_file(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, {1}, , {2}".format(protein_id, 'GENEWISE', e))
        return
    
    status_species_list = _read_failed_species(exon_genewise_path)
    
    failed_species_list = []
    for (species, data) in proteins_abinitio.items():
        if species.strip() not in status_species_list and status_species_list != []:
            continue
        protein_file    = "{0}/{1}.fa".format(directory_crawler.get_protein_path(protein_id), species)
        gene_file       = "{0}/{1}.fa".format(directory_crawler.get_gene_path(protein_id), species)
        wise_file       = "{0}/{1}.genewise".format(directory_crawler.get_genewise_path(protein_id), species)
        exon_file       = "{0}/{1}.fa".format(directory_crawler.get_exon_genewise_path(protein_id), species)
        
        wise_command    = command_generator.generate_genewise_command(protein_file, gene_file, wise_file)
        print wise_command
        wise_command_output = Popen(wise_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        time.sleep(1)
        if wise_command_output != "":
            print wise_command_output
            alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, 'GENEWISE', species.strip(), re.sub("\n", " ", wise_command_output)))
            os.remove(wise_file)
            failed_species_list.append(species.strip())
            continue
        #LOGGING
        while (True):
            if (os.path.getsize(wise_file)):
                break
        wise_file_line  = open(wise_file, 'r').readline()
        
        error_pattern = re.compile("Warning Error")
        if re.search(error_pattern, wise_file_line) is not None:
            print wise_file_line
            alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, 'GENEWISE', species.strip(), wise_file_line.strip()))
            os.remove(wise_file)
            failed_species_list.append(species.strip())
            continue
        try:
            _analyse_wise_file(protein_id, wise_file, protein_file, gene_file, exon_file)
        except IOError, e:
            alignment_logger.error("{0}, {1}, {2}, {3}".format(protein_id, 'GENEWISE', species.strip(), e))    
            failed_species_list.append(species.strip())
    if failed_species_list:
        _write_failed_species(exon_genewise_path, failed_species_list)
        return False
    return True
        
def _analyse_wise_file(protein_id, wise_file, protein_file, gene_file, exon_file):    
    exon_pattern = re.compile(r'  Exon (\d+) (\d+) phase \d+')
    gene_location_pattern = re.compile(r'>.*:(\d+)\-(\d+).*')
    
    # Reading the gene sequence
    input_DNA_f = open(gene_file, 'r');
    gene_header = input_DNA_f.readline();
    DNA = re.sub("\n", "", input_DNA_f.read());
    DNA = Seq( DNA, IUPAC.unambiguous_dna )
    input_DNA_f.close()
    #
    exon_out = open(exon_file, 'w')
    wise_out    = open(wise_file, 'r')
    exon_cnt = 0
    for line in wise_out.readlines():
        exon_pattern_match = re.match(exon_pattern, line)
        if exon_pattern_match is not None:
            exon_cnt += 1
            relative_lower_coord = int(exon_pattern_match.groups()[0]) - 1
            relative_upper_coord = int(exon_pattern_match.groups()[1])
            
            location_match = re.match(gene_location_pattern, gene_header)
            absolute_lower_coord = int(location_match.groups()[0]) + relative_lower_coord
            absolute_upper_coord = int(location_match.groups()[0]) + relative_upper_coord - 1
            
            record = SeqRecord(Seq(str(DNA[relative_lower_coord:relative_upper_coord]),
                                   IUPAC.unambiguous_dna),
                                   id=str(exon_cnt), 
                                   description="exon length {0}|{1}|{2}".format(relative_upper_coord - relative_lower_coord, absolute_lower_coord, absolute_upper_coord))
            SeqIO.write([record], exon_out, "fasta")
    exon_out.close()
    return True

def _write_failed_species(path, failed_species_list):
    '''
    Creates or updates the .status file for a respective sequence
    type with the names of species for which sequences haven't been found.
    '''
    status_file = open("{0}/.status".format(path), 'w')
    for species in failed_species_list:
        status_file.write("{0}\n".format(species))
    status_file.close()
    
    if failed_species_list == []:
        os.remove("{0}/.status".format(path))
                  
def _read_failed_species(path):
    '''
    Reads the contents of the .status file and returns the list of species
    that have failed before.
    '''
    failed_species = []
    if os.path.exists("{0}/.status".format(path)):
        status_file = open("{0}/.status".format(path), 'r')
        for line in status_file.readlines():
            failed_species.append(line.strip())
        status_file.close()
    return failed_species
  
def main ():
    populate_sequence_exon_genewise("ENSP00000311134")

if __name__ == '__main__':
    main()