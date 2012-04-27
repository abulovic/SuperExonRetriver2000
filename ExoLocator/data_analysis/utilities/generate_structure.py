'''
Created on Apr 26, 2012

@author: marioot
'''
from utilities                                      import FileUtilities
from utilities.Logger                               import Logger
from utilities.DescriptionParser                    import DescriptionParser
from data_analysis.containers.DataMapContainer      import DataMapContainer
from data_analysis.containers.ProteinContainer      import ProteinContainer
from data_analysis.containers.GeneContainer         import GeneContainer
from data_analysis.containers.TranscriptContainer   import TranscriptContainer
from data_analysis.base.DataMap                     import DataMap
from data_analysis.base.Protein                     import Protein
from data_analysis.base.Gene                        import Gene
from data_analysis.base.Transcript                  import Transcript


def load_protein_configuration(protein_id, ref_species_dict = None):
    '''
    Loads the data from a description file, and calls the containers generating functions to create basic objects.
    @param protein_id: id of a single protein
    '''
    if not check_status_file(protein_id):
        return False
    
    if ref_species_dict is None:
        ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    data_map_container  = DataMapContainer.Instance()
    protein_container   = ProteinContainer.Instance()
    gene_container      = GeneContainer.Instance()
    transcript_container= TranscriptContainer.Instance()
    
    (known_proteins, abinitio_proteins) = DescriptionParser().parse_description_file_general_info(protein_id)
    
    for species_data in known_proteins:
        (species_name, 
         spec_protein_id, 
         spec_gene_id, 
         spec_transcript_id, 
         location_type, 
         assembly, 
         location_id, 
         seq_begin, 
         seq_end, 
         strand) = species_data
        ab_initio = False
        
        data_map_key = (protein_id, species_name)
        data_map     = DataMap(spec_protein_id, spec_transcript_id, spec_gene_id, data_map_key, ab_initio)
        #TODO: ERR handling!
        data_map_container.add(data_map_key, data_map)

        protein     = Protein(spec_protein_id, data_map_key, ref_species_dict[species_name])
        gene        = Gene(spec_gene_id, data_map_key, ref_species_dict[species_name])
        transcript  = Transcript(spec_transcript_id, data_map_key, ref_species_dict[species_name])
        
        protein_container.add(protein.protein_id, protein)
        gene_container.add(gene.gene_id, gene)
        transcript_container.add(transcript.transcript_id, transcript)
    
    for species_data in abinitio_proteins:
        #create objects for abinitio proteins
        (species_name, 
         spec_protein_id, 
         location_type, 
         assembly, 
         location_id, 
         seq_begin, 
         seq_end, 
         strand) = species_data
        ab_initio = True
        
        data_map_key = (protein_id, species_name)
        data_map = DataMap(spec_protein_id, "", "", data_map_key, ab_initio)
        #TODO: ERR handling!
        data_map_container.add(data_map_key, data_map)
        
        protein     = Protein(spec_protein_id, data_map_key, ref_species_dict[species_name])
        
        protein_container.add(protein.protein_id, protein)
    return True
    
    
def load_protein_configuration_batch(protein_id_list):
    '''
    @param protein_id_list: list of protein id's
    '''
    ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('containers')
    
    folders_loaded_cnt  = 0
    try:
        for protein_id in protein_id_list:
            if load_protein_configuration(protein_id, ref_species_dict) == True:
                folders_loaded_cnt += 1
    except (KeyError, TypeError), e:
        alignment_logger.warning("{0}, {1}".format(protein_id, e))
            
    return folders_loaded_cnt
    
def check_status_file(protein_id):
    try:
        status_file_keys = FileUtilities.get_status_file_keys()
        protein_status  = FileUtilities.read_status_file(protein_id)
        for key in status_file_keys:
            if protein_status[key] == 'FAILED':
                print "{0}: Loading protein data failed due to .status file!".format(protein_id)
                return False
    except KeyError:
        print "{0}: Loading protein data failed due missing key in .status file!".format(protein_id)
        return False
    return True

def main ():
    protein_list_raw = FileUtilities.get_protein_list()
    protein_list = []
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
        
    print load_protein_configuration_batch(protein_list)
    
    dmc = DataMapContainer.Instance()
    pc  = ProteinContainer.Instance()
    gc  = GeneContainer.Instance()
    tc  = TranscriptContainer.Instance()
    
    pass
    
if __name__ == '__main__':
    main()