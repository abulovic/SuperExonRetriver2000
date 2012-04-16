'''
Created on Apr 16, 2012

@author: intern
'''
from utilities.FileUtilities import get_protein_list, get_species_list
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
import os
from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator
from pipeline.ortholog_search.OrthologFinder import find_ortholog_by_RBH
from utilities.Logger import Logger

def main():
    
    dc = DirectoryCrawler()
    acg = AlignmentCommandGenerator()
    logger = Logger.Instance()
    mutual_best_logger = logger.get_logger('mutual_best')
    
    protein_list = get_protein_list()
    species_list = get_species_list()
    
    for (protein_id, num_of_exons) in protein_list:
        
        dc.generate_directory_tree(protein_id)
        mutual_best_logger.info("\tPROTEIN\t\t%s" % protein_id)
        
        descr_file_path = dc.get_protein_description_file_path(protein_id)
        if (os.path.isfile(descr_file_path) and os.path.getsize(descr_file_path)):
            mutual_best_logger.info("Mutual best already exists, moving on to the next protein.")
            continue
        
        descr_file = open(descr_file_path, 'w')
        
        ref_species_pep =  dc.get_protein_path(protein_id) + "/" + "Homo_sapiens.fasta"
        fastacmd = acg.generate_fastacmd_protein_command(protein_id, "Homo_sapiens", "all", ref_species_pep)
        os.system(fastacmd)
             
        
        for species in species_list:
            find_ortholog_by_RBH("Homo_sapiens", species, ref_species_pep, protein_id, descr_file, mutual_best_logger)
            
        descr_file.close()
        mutual_best_logger.info("\n\n")
        

if __name__ == '__main__':
    main()