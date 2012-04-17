'''
Created on Apr 16, 2012

@author: intern
'''

import os
from subprocess import Popen, PIPE, STDOUT

from utilities.FileUtilities import get_protein_list, get_species_list,\
    get_protein_ids, read_status_file, update_entry_in_status_file
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator
from pipeline.ortholog_search.OrthologFinder import find_ortholog_by_RBH
from utilities.Logger import Logger

def main():
    
    reference_species = "Homo_sapiens"
    
    dc = DirectoryCrawler()
    acg = AlignmentCommandGenerator()
    
    logger = Logger.Instance()
    mutual_best_logger = logger.get_logger('mutual_best')
    
    protein_list = get_protein_list()
    species_list = get_species_list()
    failed_proteins = []
    
    for (protein_id, num_of_exons) in protein_list:
        
        known_dict = {}
        abinitio_dict = {}
        print protein_id
        dc.generate_directory_tree(protein_id)
        
        descr_file_path = dc.get_protein_description_file_path(protein_id)
        print descr_file_path
        status_file_path = dc.get_mutual_best_status_file_path(protein_id)
        
        if (os.path.isfile(status_file_path) and os.path.getsize(status_file_path)):
            print get_protein_ids(protein_id)
            
            status_dict = read_status_file(protein_id)
            if (status_dict.has_key('MUTUAL_BEST')):
                if status_dict['MUTUAL_BEST'] == 'OK':
                    mutual_best_logger.info('-,%s,mutual_best already exists for this protein - moving to the next one' % protein_id)
                else :
                    mutual_best_logger.error('-,%s,mutual_best has failed for this protein (no orthologs found) - moving on the next one' % protein_id)
                    failed_proteins.append(protein_id)
            continue
        
        
        descr_file = open(descr_file_path, 'w')
        
        ref_species_pep =  dc.get_protein_path(protein_id) + "/" + reference_species + ".fasta"
        fastacmd = acg.generate_fastacmd_protein_command(protein_id, "Homo_sapiens", "all", ref_species_pep)
        
        p = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if output:
            mutual_best_logger.error("%s,fastacmd error" % protein_id)
             
        
        for species in species_list:
            find_ortholog_by_RBH("Homo_sapiens", species, ref_species_pep, protein_id, descr_file, mutual_best_logger)
            
        descr_file.close()
        
        mutual_best_logger.info("\n\n")
        
        # check what we've found out, whether this protein has any orthologs
        (known_dict, abinitio_dict) = get_protein_ids(protein_id)
        if (not abinitio_dict and (not known_dict or (len(known_dict.keys()) == 1 and known_dict.keys()[0] == reference_species))):
            mutual_best_logger.info ("-,%s, mutual best failed for this protein." % protein_id)
            update_entry_in_status_file(protein_id, "MUTUAL_BEST", "FAILED")
            failed_proteins.append(protein_id)
            
        else:
            update_entry_in_status_file(protein_id, "MUTUAL_BEST", "OK")
            
            
    for failed_protein_id in failed_proteins:
        print failed_protein_id

if __name__ == '__main__':
    main()