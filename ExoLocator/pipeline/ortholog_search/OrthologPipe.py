'''
Created on Apr 16, 2012

@author: ana, mario
'''

# Python import
import os
from subprocess                                     import Popen, PIPE, STDOUT

# script specific import
from pipeline.ortholog_search.OrthologFinder        import find_ortholog_by_RBH

# utilities import
from utilities.FileUtilities                        import get_protein_list, get_default_species_list,\
                                                        read_status_file, update_entry_in_status_file
from utilities.DescriptionParser                    import DescriptionParser
from utilities.Logger                               import Logger
from utilities.DirectoryCrawler                     import DirectoryCrawler

# pipeline utilites import
from pipeline.utilities.AlignmentCommandGenerator   import AlignmentCommandGenerator




def main():
    
    '''
    Retrieves the list of all the proteins from reference species.
    For each ref species protein, it tries to find orthologues for all the species (from the species list)
    and generates the description file accordingly. If the description file already exists, it checks
    the status (OK/PARTIAL/FAILED).
    '''
    
    reference_species = "Homo_sapiens"
    
    dc = DirectoryCrawler()
    acg = AlignmentCommandGenerator()
    
    logger = Logger.Instance()
    mutual_best_logger = logger.get_logger('mutual_best')
    
    protein_list = get_protein_list()
    species_list = get_default_species_list()
    failed_proteins = []
    
    for (protein_id, num_of_exons) in protein_list:
        
        known_dict = {}
        abinitio_dict = {}
        print protein_id
        
        # generate all the directories for the protein
        dc.generate_directory_tree(protein_id)
        
        descr_file_path = dc.get_protein_description_file_path(protein_id)
        status_file_path = dc.get_mutual_best_status_file_path(protein_id)
        
        if (os.path.isfile(status_file_path) and os.path.getsize(status_file_path)):
            print DescriptionParser().get_protein_ids(protein_id)
            
            status_dict = read_status_file(protein_id)
            if (status_dict.has_key('MUTUAL_BEST')):
                if status_dict['MUTUAL_BEST'] == 'OK':
                    mutual_best_logger.info('-,%s,mutual_best already exists for this protein - moving to the next one' % protein_id)
                else :
                    mutual_best_logger.error('-,%s,mutual_best has failed for this protein (no orthologs found) - moving on the next one' % protein_id)
                    failed_proteins.append(protein_id)
            continue
        
        
        # create the description file
        descr_file = open(descr_file_path, 'w')
        # reference protein file
        ref_species_pep =  dc.get_protein_path(protein_id) + "/" + reference_species + ".fasta"
        fastacmd = acg.generate_fastacmd_protein_command(protein_id, reference_species, "all", ref_species_pep)
        
        p = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if output:
            mutual_best_logger.error("%s,fastacmd error" % protein_id)
             
        # find orthologues for all species
        for species in species_list:
            find_ortholog_by_RBH (reference_species, species, ref_species_pep, protein_id, descr_file, mutual_best_logger)
            
        descr_file.close()
        
        mutual_best_logger.info("\n\n")
        
        # check what we've found out, whether this protein has any orthologs
        (known_dict, abinitio_dict) = DescriptionParser().get_protein_ids(protein_id)
        if (not abinitio_dict and (not known_dict or (len(known_dict.keys()) == 1 and known_dict.keys()[0] == reference_species))):
            mutual_best_logger.info ("-,%s, mutual best failed for this protein." % protein_id)
            update_entry_in_status_file(protein_id, "MUTUAL_BEST", "FAILED")
            failed_proteins.append(protein_id)
            
        else:
            update_entry_in_status_file(protein_id, "MUTUAL_BEST", "OK")
            
    print "Failed proteins: "        
    for failed_protein_id in failed_proteins:
        print failed_protein_id

if __name__ == '__main__':
    main()