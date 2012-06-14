'''
Created on Apr 16, 2012

@author: ana, mario
'''

# Python imports
import os, re
# BioPython imports
from Bio.Blast                                      import NCBIXML

# utilities imports
from utilities.FileUtilities                        import get_default_species_list, execute_command_and_log
from utilities.DirectoryCrawler                     import DirectoryCrawler

# pipeline utilities imports
from pipeline.utilities.AlignmentCommandGenerator   import AlignmentCommandGenerator





def find_ortholog_by_RBH (original_species, target_species, original_protein_fasta, original_protein_id, descr_file, mutual_best_logger = None):
    
    '''
    Searches for the protein orthologue of a reference species protein.
    Writes the found orthologue data to the description file.
    @param original_species:           latin name of the reference species
    @param target_species:             latin name of  the target species
    @param original_protein_fasta:     fasta file containing the reference species protein
    @param original_protein_id:        reference species protein ensembl ID
    @param descr_file:                 description file handle
    @param mutual_best_logger:         mutual best logger (writes in a file)
    
    Available data for each protein:
    - known / novel proteins:
        (species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    - ab initio (genscan)
        (species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
    '''
    
    # first try to find the orthologue among the known/novel proteins:
    protein_info_known = _search_for_ortholog_in_database(original_species, target_species, 
                                                          original_protein_fasta, original_protein_id, 
                                                          "all", mutual_best_logger)
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
        # now try to find the orthologue among the ab initio proteins
        mutual_best_logger.info("%s,%s, match not found in all - trying abinitio..." % (target_species.upper(), original_protein_id))
        print "%s, match not found in all - trying abinitio" % target_species.upper()
        
        protein_info_abinitio = _search_for_ortholog_in_database(original_species, target_species, 
                                                                 original_protein_fasta, original_protein_id, 
                                                                 "abinitio", mutual_best_logger)
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
            
            
    
 
def _search_for_ortholog_in_database(original_species, target_species, original_protein_fasta, original_protein_id, db_type, logger):
    '''
    Takes the reference species protein and makes a BLASTp query on the target_species database.
    Provided that there has been at least one hit, the best BLAST hit is queried against the
    reference species protein database. If the best hit from the second query is the original protein,
    then the protein data is returned. The queries as referred to as the forward and the backward
    hit respectively.  
    @param db_type: all / abinitio
    '''
    acg = AlignmentCommandGenerator()
    
    if (db_type == "all"):
        protein_id_pattern = re.compile("lcl\|(.*)\spep:(.*)\s(.*):(.*):(.*):(.*):(.*):(.*)\sgene:(.*)\stranscript:(.*)\s.*\s.*")
    else:
        protein_id_pattern = re.compile("lcl\|(.*)\spep:(.*)\s(.*):(.*):(.*):(.*):(.*):(.*)\stranscript:(.*)\s.*")
    
    output_file = "tmp.xml"
    forward_blastp_cmd = acg.generate_blastp_command_for_species(target_species, original_protein_fasta, output_file, db_type)
    print forward_blastp_cmd
    
    execute_command_and_log(logger, forward_blastp_cmd, (target_species, original_protein_id))
    
    result_handle = open(output_file)  
    blast_records = NCBIXML.parse(result_handle)
    
    try:
        best_forward_hit = blast_records.next()
    except (ValueError):
        logger.error("%s,%s,XML file empty - no forward blast results" % (target_species, original_protein_id))
        return None
        
    if (best_forward_hit.alignments):
        bfh_title = _get_best_alignment(original_protein_id, best_forward_hit)
        protein_match = re.match(protein_id_pattern, bfh_title)
    else:
        return None
    
    if (db_type == "all"):
        (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, gene_id, transcript_id) = protein_match.groups()
    else:
        (protein_id, protein_type, location_type, assembly, location_id, seq_start, seq_end, strand, transcript_id) = protein_match.groups()
    
    fasta_input_file = "species.fasta"
    fastacmd = acg.generate_fastacmd_protein_command(protein_id, target_species, db_type, fasta_input_file)
    
    execute_command_and_log(logger, fastacmd, (target_species, original_protein_id))
    
    backward_blastp_cmd = acg.generate_blastp_command_for_species(original_species, fasta_input_file, output_file, "all")
    print backward_blastp_cmd
    
    execute_command_and_log(logger, backward_blastp_cmd, (target_species, original_protein_id))
    
    result_handle = open(output_file)  
    blast_records = NCBIXML.parse(result_handle)
    
    try:
        best_backward_hit = blast_records.next()
    except (ValueError):
        logger.error("%s,%s,XML file empty - no backward blast results" % (target_species, original_protein_id))
        return None
    
    if (best_backward_hit.alignments):
        bbh_title = _get_best_alignment(original_protein_id, best_backward_hit)
        protein_match_b = re.match(protein_id_pattern, bbh_title)
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
  
  
    
def _get_best_alignment (protein_id, blast_record):
    '''
    Auxiliary function - because of the presence of different haplotype 
    proteins which can have the same BLAST score (or very similar), not
    only the first BLAST hit is taken into account, but first few.
    '''
    
    best_score = blast_record.descriptions[0].score
    best_alignments = [blast_record.descriptions[0].title]
    i = 1
    '''
    Get the best alignments
    '''
    while i < len(blast_record.descriptions):
        
        if (blast_record.descriptions[i].score >= best_score):
            best_alignments.append(blast_record.descriptions[i].title)
        else:
            break
        i += 1
    
    pattern = re.compile("lcl\|(.*)\spep::*")
    for title in best_alignments:
        prot_match = re.match(pattern, title)
        if prot_match.groups()[0] == protein_id:
            return title
    return best_alignments[0]
    
    
    
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
    
    for species in get_default_species_list():
        find_ortholog_by_RBH("Homo_sapiens", species, output_file_path, protein_id)
        
    descr_file.close()