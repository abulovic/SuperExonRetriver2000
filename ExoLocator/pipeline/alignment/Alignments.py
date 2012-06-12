'''
Created on Apr 13, 2012

@author: marioot
'''
# Python imports
import os
from subprocess                                   import Popen, PIPE, STDOUT

# BioPython imports
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# utilities imports
from utilities.DirectoryCrawler                   import DirectoryCrawler
from utilities.Logger                             import Logger
from utilities.DescriptionParser                  import DescriptionParser

# pipeline imports
from pipeline.utilities.CommandGenerator          import CommandGenerator
from pipeline.alignment.AlignmentTargetGenerator  import AlignmentTargetGenerator

# data analysis imports
from data_analysis.containers.ExonContainer       import ExonContainer


'''
Hepler functions for generating alignments (blastn, tblastn, SW on gene, SW on cDNA)
'''

def generate_blastn_alignments(protein_id, species_list = None, referenced_species = "Homo_sapiens"):
    '''
        Runs the blastn program for a specified protein and list of species
        @param protein_id
        @param species_list: if provided, runs blastn for this list of species, \
                             otherwise runs for species that are missing the blastn output \
                             who are determined by .status file in the blastn folder.
    '''
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('alignment')
    
    crawler             = DirectoryCrawler()
    
    command_generator   = CommandGenerator()
    alignment_generator = AlignmentTargetGenerator()
    
    failed_species_list = []
    
    # retrieve the blastn targets
    if (not species_list):
        species_list    = alignment_generator.get_blastn_targets(protein_id)

    for species in species_list:
        
        ############# MOVE TO ANOTHER FNC
        output_file     = "{0}/{1}.blastout".format(crawler.get_blastn_path(protein_id), species.strip())
        input_file      = "{0}/{1}.fa".format(crawler.get_expanded_gene_path(protein_id), species.strip())
        database        = "{0}/{1}.fa".format(crawler.get_database_path(protein_id), referenced_species)
        
        command         = command_generator.generate_blastn_command(database, input_file, output_file)
        command_return  = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            os.remove(output_file)
            alignment_logger.warning("{0}, {1}, BLASTN, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species.strip())
            
    if failed_species_list:
        alignment_generator.set_failed_blastn_targets(protein_id, failed_species_list)
        return False
    return True

def generate_tblastn_alignments(protein_id, species_list = None, referenced_species = "Homo_sapiens"):
    '''
        Runs the tblastn program for a specified protein and list of species
        @param protein_id
        @param species_list: if provided, runs tblastn for this list of species, \
                             otherwise runs for species that are missing the tblastn output \
                             who are determined by .status file in the tblastn folder.
    '''
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('alignment')
    
    alignment_generator = AlignmentTargetGenerator()
    crawler             = DirectoryCrawler()
    command_generator   = CommandGenerator()
    
    if (not species_list):
        species_list    = alignment_generator.get_tblastn_targets(protein_id)
    
    failed_species_list = []
    for species in species_list:
        
        ############## MOVE
        output_file     = "{0}/{1}.blastout".format(crawler.get_tblastn_path(protein_id), species.strip())
        input_file      = "{0}/{1}.fa".format(crawler.get_protein_path(protein_id), species.strip())
        database        = "{0}/{1}.fa".format(crawler.get_database_path(protein_id), referenced_species)
        
        command         = command_generator.generate_tblastn_command(database, input_file, output_file)
        command_return  = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            os.remove(output_file)
            alignment_logger.warning("{0}, {1}, TBLASTN, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species.strip())
      
    if failed_species_list:        
        alignment_generator.set_failed_tblastn_targets(protein_id, failed_species_list)
        return False
    return True


def generate_SW_gene_alignments(protein_id, species_list = None, referenced_species = "Homo_sapiens"):
    '''
        Runs the SW program for a specified protein and list of species, using the expanded gene region.
        @param protein_id
        @param species_list: if provided, runs SW for this list of species, \
                             otherwise runs for species that are missing the SW output \
                             who are determined by .status file in the /SW/gene folder.
    '''       
    logger                   = Logger.Instance()
    alignment_logger         = logger.get_logger('alignment')
     
    alignment_generator      = AlignmentTargetGenerator()
    crawler                  = DirectoryCrawler()
    command_generator        = CommandGenerator()
    
    if (not species_list):
        species_list         = alignment_generator.get_SW_gene_targets(protein_id)

    failed_species_list = []
    for species in species_list:
        
        ########### MOVE
        output_file          = "{0}/{1}.swout".format(crawler.get_SW_gene_path(protein_id), species.strip())
        query_sequence_file  = "{0}/{1}.fa".format(crawler.get_expanded_gene_path(protein_id), species.strip())
        target_fasta_db_file = "{0}/{1}.fa".format(crawler.get_database_path(protein_id), referenced_species)
        
        command              = command_generator.generate_SW_command(query_sequence_file, target_fasta_db_file, output_file)
        command_return       = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output               = command_return.stdout.read()
        if output != "":
            #LOGGING
            alignment_logger.warning("{0}, {1}, SW GENE, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species.strip())
    os.remove(".sw_stdout_supressed")
    
    if failed_species_list: 
        alignment_generator.set_failed_SW_gene_targets(protein_id, failed_species_list)
        return False
    return True
  

def _merge_exons_from_fasta(exon_fasta_file, species):
    '''
    Auxiliary function. Merges all the exons to a single "cDNA" sequence
    and writes down the original exon locations
    '''
    
    exon_fasta = open(exon_fasta_file, 'r')
    merged_exons_seq=""
    exon_locations = {}
    start = 1
    end = 1
    exon_id = 1
    for seq_record in SeqIO.parse(exon_fasta, "fasta"):
        end += len(seq_record)-1
        exon_locations[exon_id] = (start, end)
        exon_id += 1
        start = end+1
        end = start
        merged_exons_seq += seq_record.seq

    cdna_seq = SeqRecord(seq=merged_exons_seq, id=species, description="")
    exon_fasta.close()
    return (cdna_seq, exon_locations)
    


def _write_locations_to_file(locations_file_path, exon_locations_original, exon_locations_target):
    '''
    Auxiliary function. Writes the exon locations to a file.
    '''
    
    locations_file = open(locations_file_path, 'w')
    
    locations_file.write ('# format: R/T (reference / target speacies) exon_id start end\n')
    for exon_num in range(1, len(exon_locations_original)+1):
        (start,end) = exon_locations_original[exon_num]
        locations_file.write ("R\t%d\t%d\t%d\n" % (exon_num, int(start), int(end)))
    for exon_num in range(1, len(exon_locations_target)+1):
        (start,end) = exon_locations_target[exon_num]
        locations_file.write ("T\t%d\t%d\t%d\n" % (exon_num, int(start), int(end)))    
        
    locations_file.close()


######## THIS IS THE ONE
def generate_SW_exon_alignments2 (protein_id, species_list = None, referenced_species = "Homo_sapiens"):
    
    # utilities
    alignment_generator = AlignmentTargetGenerator()
    crawler             = DirectoryCrawler()
    command_generator   = CommandGenerator()
    
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('alignment')
    
    exon_container      = ExonContainer.Instance()
    
    tmp_fasta_target_path   = "tmp_target.fa"
    ref_exons_fasta = "%s/%s.fa" % (crawler.get_database_path(protein_id), referenced_species)
    
    #ref_exons = exon_container.get((protein_id, referenced_species, "ensembl"))
    #ref_exons.export_coding_exons_to_fasta(tmp_ref_exons_fasta_path)
    
    failed_species_list = []
    
    if (not species_list):
        species_list    = alignment_generator.get_SW_exon_targets(protein_id)
        
    try:
        (proteins_known, proteins_abinitio) = DescriptionParser().parse_descr_file(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, , SW cDNA_EXONS, {2}".format(protein_id, e))
        return False
             
    if not species_list:
        return False  
    for species in species_list:
        
        if species.strip() in proteins_known:
            exon_type = "ensembl"  
        else:
            exon_type = "genewise"
        
        try:
            exons = exon_container.get((protein_id, species, exon_type))
            
        except KeyError:
            print "{0}, {1}, SW cDNA_EXONS, {2}".format(protein_id, species.strip(), "Target species exon file missing")
            alignment_logger.warning("{0}, {1}, SW cDNA_EXONS, {2}".format(protein_id, species.strip(), "Target species exon file missing"))
            failed_species_list.append(species.strip())
            continue
            
        
        cDNA = exons.get_coding_cDNA()
        cDNA_record = SeqRecord(cDNA, id = species, description="coding cDNA, no UTR")
        
        
        # write merged sequence to file  
        try:
            tmp_fasta_target = open(tmp_fasta_target_path, 'w')
            SeqIO.write([cDNA_record], tmp_fasta_target, "fasta")
            tmp_fasta_target.close()
        except TypeError, e:
            alignment_logger.error("{0}, {1}, SW cDNA_EXONS, {2}".format(protein_id, species.strip(), e))
            failed_species_list.append(species.strip())    
            continue
        
        swout_file_path = "{0}/{1}.swout".format(crawler.get_SW_exon_path(protein_id), species.strip())
        #command         = command_generator.generate_SW_command(tmp_fasta_target_path, original_fasta_db_file, swout_file_path)
        
        command         = command_generator.generate_SW_command(tmp_fasta_target_path, ref_exons_fasta, swout_file_path)
        print command
        
        command_return  = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output          = command_return.stdout.read()
        if output != "":
            #LOGGING
            alignment_logger.warning("{0}, {1}, SW cDNA_EXONS, {2}".format(protein_id, species.strip(), output.strip()))
            failed_species_list.append(species.strip())    
            continue

    try:    
        os.remove(tmp_fasta_target_path)
        os.remove(".sw_stdout_supressed")
    except:
        pass
    
    if failed_species_list: 
        alignment_generator.set_failed_SW_exon_targets(protein_id, failed_species_list)
        return False
    return True



def generate_referenced_species_database(protein_id, referenced_species):
    '''
        Creates a database for a referenced species and protein, using formatdb
        @param protein_id
        @param referenced_species
    ''' 
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('alignment')
    
    command_generator   = CommandGenerator()
    crawler             = DirectoryCrawler()
    
    exon_container      = ExonContainer.Instance()
    
    input_exons = exon_container.get((protein_id, referenced_species, "ensembl"))
    
    #source_exon_file    = "{0}/{1}.fa".format(crawler.get_exon_ensembl_path(protein_id), referenced_species)
    input_db_file       = "{0}/{1}.fa".format(crawler.get_database_path(protein_id), referenced_species)
    sequence_type       = "Nucleotide"
    
    input_exons.export_coding_exons_to_fasta(input_db_file)
    
    command             = command_generator.generate_formatdb_command(input_db_file, sequence_type)
    command_return      = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output              = command_return.stdout.read()
    if output != "":
        #LOGGING
        alignment_logger.warning("{0}, {1}, REF SPECIES DB, {2}".format(protein_id, referenced_species.strip(), output.strip()))      
        return False
    return True
    
def main ():
    generate_referenced_species_database("ENSP00000311134", "Homo_sapiens")
    
if __name__ == '__main__':
    main()