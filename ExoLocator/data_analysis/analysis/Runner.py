'''
Created on Jun 14, 2012

@author: anana
'''
from utilities.FileUtilities                        import get_protein_list, check_status_file, update_entry_in_status_file, read_status_file
    
from data_analysis.utilities.generate_structure     import fill_all_containers
from data_analysis.containers.ExonContainer         import ExonContainer
from data_analysis.analysis.TranscriptionMachinery  import translate_alignment_exons_for_protein
from utilities.DirectoryCrawler import DirectoryCrawler
from data_analysis.analysis.AlignmentStatistics import create_protein_statistics

def translate_alignment_exons ():
    '''
    For the protein in the protein list file does the following:
        - check if the status file is ok. If not, it writes the failed status of translation
        - if the status is ok, it checks if the translation status is already OK
        - if the translation status is OK, then it just continues to the next protein
        - if the status is FAILED or PARTIAL, it tries to translate exons to proteins 
          for all the species for which it is necessary (meaning the translated
          protein hasn't already been generated). 
    '''

    protein_list = get_protein_list()
    
    for (protein_id, exon_num) in protein_list:
        
        if not check_status_file(protein_id):
            print "ABORTING {0} TRANSLATION: some resources have FAILED stats!".format(protein_id)
            update_entry_in_status_file(protein_id, 'EXON_TRANSLATION', 'FAILED')
            continue
        try:
            if read_status_file(protein_id)['EXON_TRANSLATION'] == 'OK':
                print "SKIPPING {0} TRANSLATION: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "TRANSLATING EXONS: {0}".format(protein_id)
        failed_species = translate_alignment_exons_for_protein(protein_id, exon_num)
        if failed_species:
            update_entry_in_status_file(protein_id, 'EXON_TRANSLATION', 'PARTIAL')
        else:
            update_entry_in_status_file(protein_id, 'EXON_TRANSLATION', 'OK')
            
def create_statistics(protein_list):
    dc = DirectoryCrawler()

    for (protein_id, exon_num) in protein_list:
    
        stat_file = "%s/stats.csv" % dc.get_root_path(protein_id)
        if not check_status_file(protein_id):
            continue
        create_protein_statistics(protein_id, stat_file)
        
def main ():
    fill_all_containers(True)
    ec = ExonContainer.Instance()
    #translate_alignment_exons()
    
    create_statistics(get_protein_list())
    

if __name__ == '__main__':
    main()
