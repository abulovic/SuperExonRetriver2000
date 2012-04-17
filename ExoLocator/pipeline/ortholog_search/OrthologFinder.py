'''
Created on Apr 16, 2012

@author: intern
'''

import os, re

from Bio.Blast import NCBIXML

from utilities.ConfigurationReader import ConfigurationReader
from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from utilities.FileUtilities import get_species_list



def find_ortholog_by_RBH (original_species, target_species, original_protein_fasta, original_protein_id, descr_file, mutual_best_logger = None):

    #(species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    #(species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    
    protein_info_known = _search_for_ortholog_in_database(original_species, target_species, original_protein_fasta, original_protein_id, "all")
    if protein_info_known:
        (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, gene_id, transcript_id) = protein_info_known
        mutual_best_logger.info("%s,%s,match found in all: location %s %s,%s" % 
                                (target_species.upper(), original_protein_id, location_type, location_id, protein_id))
        print "%s,%s,match found in all: location %s %s,%s" % (target_species.upper(), protein_id, location_type, location_id, protein_id)
                                
        
        descr_file.write("{0}\t{1}\t{2}\t{3}\t{4}:{5}:{6}:{7}:{8}:{9}\n".format(target_species, 
                                                                              protein_id, 
                                                                              gene_id, 
                                                                              transcript_id, 
                                                                              location_type, 
                                                                              assembly, 
                                                                              location_id, 
                                                                              seq_start, 
                                                                              seq_end, 
                                                                              strand))
    else:
        mutual_best_logger.info("%s,%s, match not found in all - trying abinitio..." % (target_species.upper(), original_protein_id))
        print "%s, match not found in all - trying abinitio" % target_species.upper()
        protein_info_abinitio = _search_for_ortholog_in_database(original_species, target_species, original_protein_fasta, original_protein_id, "abinitio")
        if protein_info_abinitio:
            (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, transcript_id) = protein_info_abinitio
            print "%s, match found in abinitio: %s type %s, location type %s" % (target_species.upper(), protein_id, protein_type, location_type)
            mutual_best_logger.info("%s,%s,match found in abinitio: location %s %s,%s" % 
                                    (target_species.upper(), original_protein_id, location_type, location_id, protein_id))
            descr_file.write("{0}\t{1}\t{2}:{3}:{4}:{5}:{6}:{7}\n".format(target_species, 
                                                                              protein_id, 
                                                                              location_type, 
                                                                              assembly, 
                                                                              location_id, 
                                                                              seq_start, 
                                                                              seq_end, 
                                                                              strand))
        else:
            print "%s, protein not found in abinitio either." % target_species.upper()
            mutual_best_logger.info("%s,%s,protein not found in abinitio either." % (target_species.upper(), original_protein_id))
    
 
def _search_for_ortholog_in_database(original_species, target_species, original_protein_fasta, original_protein_id, db_type):
    '''
    @param db_type: all / abinitio
    '''
    acg = AlignmentCommandGenerator()
    
    if (db_type == "all"):
        protein_id_pattern = re.compile("lcl\|(.*)\spep:(.*)\s(.*):(.*):(.*):(.*):(.*):(.*)\sgene:(.*)\stranscript:(.*)\s.*\s.*")
    else:
        protein_id_pattern = re.compile("lcl\|(.*)\spep:(.*)\s(.*):(.*):(.*):(.*):(.*):(.*)\stranscript:(.*)\s.*")
    
    output_file = "tmp.xml"
    forward_blastp_cmd = acg.generate_blastp_command_for_species(target_species, original_protein_fasta, output_file, db_type)
    os.system(forward_blastp_cmd)
    
    result_handle = open(output_file)  
    blast_records = NCBIXML.parse(result_handle)
    
    best_forward_hit = blast_records.next()
    if (best_forward_hit.alignments):
        bfh_alignment = best_forward_hit.alignments[0]
        protein_match = re.match(protein_id_pattern, bfh_alignment.title)
    else:
        return None
    
    if (db_type == "all"):
        (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, gene_id, transcript_id) = protein_match.groups()
    else:
        (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, transcript_id) = protein_match.groups()
    
    fasta_input_file = "species.fasta"
    fastacmd = acg.generate_fastacmd_protein_command(protein_id, target_species, db_type, fasta_input_file)
    os.system(fastacmd)
    
    backward_blastp_cmd = acg.generate_blastp_command_for_species(original_species, fasta_input_file, output_file, "all")
    os.system(backward_blastp_cmd)
    
    result_handle = open(output_file)  
    blast_records = NCBIXML.parse(result_handle)
    
    best_backward_hit = blast_records.next()
    if (best_backward_hit.alignments):
        bbh_alignment = best_backward_hit.alignments[0]
        protein_match_b = re.match(protein_id_pattern, bbh_alignment.title)
        protein_id_b = protein_match_b.groups()[0]
    else:
        return None
    
    os.remove(output_file)
    os.remove(fasta_input_file)
    
    if (original_protein_id == protein_id_b):
        if (db_type == "all"):
            return (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, gene_id, transcript_id)
        else:
            return (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, transcript_id)
    else:
        return None
    
    
if __name__ == '__main__':
    protein_id = "ENSP00000311134"
    acg = AlignmentCommandGenerator()
    dc = DirectoryCrawler()
    
    dc.generate_directory_tree(protein_id)
    descr_file_path = dc.get_protein_description_file_path(protein_id)
    descr_file = open(descr_file_path, 'w')
    
    output_file_path = dc.get_protein_path(protein_id) + "/" + "Homo_sapiens.fasta"
    
    fastacmd = acg.generate_fastacmd_protein_command(protein_id, "Homo_sapiens", "all", output_file_path)
    os.system(fastacmd)
    
    for species in get_species_list():
        find_ortholog_by_RBH("Homo_sapiens", species, output_file_path, protein_id)
        
    descr_file.close()