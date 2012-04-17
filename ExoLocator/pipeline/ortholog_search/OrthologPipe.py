'''
Created on Apr 16, 2012

@author: intern
'''

import os
from subprocess import Popen, PIPE, STDOUT

from utilities.FileUtilities import get_protein_list, get_species_list,\
    get_protein_ids
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
        
        dc.generate_directory_tree(protein_id)
        
        descr_file_path = dc.get_protein_description_file_path(protein_id)
        if (os.path.isfile(descr_file_path) and os.path.getsize(descr_file_path)):
            (known_dict, abinitio_dict) = get_protein_ids(protein_id)
            # mutual best contains only the reference species protein - means there are no orthologs according
            # to the rbh criterion 
            if (not abinitio_dict and (len(known_dict.keys()) == 1 and known_dict.keys()[0] == reference_species)):
                mutual_best_logger.info ("%s, mutual best failed for this protein." % protein_id)
                failed_proteins.append(protein_id)
            else:
                mutual_best_logger.info("%s,Mutual best already exists, moving on to the next protein." % protein_id)
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
        

if __name__ == '__main__':
    main()