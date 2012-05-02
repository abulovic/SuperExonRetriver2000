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
from data_analysis.base.EnsemblExons                import EnsemblExons
from data_analysis.base.Exons						import Exons
from data_analysis.containers.EnsemblExonContainer  import EnsemblExonContainer
from data_analysis.containers.ExonContainer import ExonContainer



def load_protein_configuration(protein_id, ref_species_dict = None):
    '''
    Loads the data from a description file, and calls the containers generating functions to create basic objects.
    @param protein_id: id of a single protein
    '''
    if not check_status_file(protein_id):
        return False
    
    if ref_species_dict is None:
        ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('containers')
    
    data_map_container  = DataMapContainer.Instance()
    protein_container   = ProteinContainer.Instance()
    gene_container      = GeneContainer.Instance()
    transcript_container= TranscriptContainer.Instance()
    ens_exon_container  = EnsemblExonContainer.Instance()
    
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
        try:
            data_map_key = (protein_id, species_name)
            data_map     = DataMap(spec_protein_id, spec_transcript_id, spec_gene_id, data_map_key, ab_initio)
            data_map_container.add(data_map_key, data_map)
            
            protein     = Protein(spec_protein_id, data_map_key, ref_species_dict[species_name])
            gene        = Gene(spec_gene_id, data_map_key, ref_species_dict[species_name])
            transcript  = Transcript(spec_transcript_id, data_map_key, ref_species_dict[species_name])
            ens_exons   = EnsemblExons(data_map_key, ref_species_dict[species_name])
            ens_exons.load_exons()
            
            for exon in ens_exons.exons:
                ens_exon_container.add (exon.exon_id, exon)
            
            protein_container.add(protein.protein_id, protein)
            gene_container.add(gene.gene_id, gene)
            transcript_container.add(transcript.transcript_id, transcript)
            
        except (KeyError, TypeError), e:
            alignment_logger.warning("{0}, {1}, {2}".format(protein_id, species_name, e.args[0]))
    
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
        
        try:
            data_map_key = (protein_id, species_name)
            data_map = DataMap(spec_protein_id, "", "", data_map_key, ab_initio)
            #TODO: ERR handling!
            data_map_container.add(data_map_key, data_map)
            
            protein     = Protein(spec_protein_id, data_map_key, ref_species_dict[species_name])
            
            protein_container.add(protein.protein_id, protein)
        except (KeyError, TypeError), e:
            alignment_logger.warning("{0}, {1}, {2}".format(protein_id, species_name, e.args[0]))
    return True
  

def load_exon_configuration2 (ref_protein_id, ref_species_dict=None):
    
    if not check_status_file(ref_protein_id):
        return False
    
    if not ref_species_dict:
        ref_species_dict = FileUtilities.get_reference_species_dictionary()
        
    logger = Logger.Instance()
    alignment_logger = logger.get_logger('containers')
    
    exon_container = ExonContainer.Instance()
    
    (known_proteins, abinitio_proteins) = DescriptionParser().parse_description_file_general_info(ref_protein_id)
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
         
        ref_species = ref_species_dict[species_name]
        exons = Exons((ref_protein_id, species_name), ref_species, "blastn")
        exon_dict = exons.load_exons()
        if not exon_dict:
            continue
        exons.set_reference_exons()
        
        data_map_key = [ref_protein_id, species_name]
        exon_container.add("blastn", data_map_key, exons)



def load_protein_configuration_batch(protein_id_list):
    '''
    @param protein_id_list: list of protein id's
    '''
    ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    folders_loaded_cnt  = 0
    for protein_id in protein_id_list:
        if load_protein_configuration(protein_id, ref_species_dict) == True:
            folders_loaded_cnt += 1
    return folders_loaded_cnt


def load_exon_configuration_batch(protein_id_list, ens_exon_container):
    
    ref_species_dict = FileUtilities.get_reference_species_dictionary()
   
    folders_loaded_cnt  = 0
    for protein_id in protein_id_list:
        if load_exon_configuration2(protein_id, ref_species_dict) == True:
            folders_loaded_cnt += 1
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
        
    ens_exon_container = load_protein_configuration_batch(protein_list)
    if ens_exon_container:
        load_exon_configuration_batch (protein_list, ens_exon_container)
    
    dmc = DataMapContainer.Instance()
    pc  = ProteinContainer.Instance()
    gc  = GeneContainer.Instance()
    tc  = TranscriptContainer.Instance()
    eec = EnsemblExonContainer.Instance()
    ec  = ExonContainer.Instance()
    
    for exons in ec._exon_container.values():
        print "EXONS...", exons.exon_type, exons.species
        for exon in exons.alignment_exons.values():
            for ref_exon_id in exon.reference_exons:
                al_values = exon.reference_exons[ref_exon_id]
                print "  EXON ID: ", ref_exon_id
                for al_value in al_values:
                    print "\tALIGNMENT: ", al_value["start"], al_value["stop"], "ordinal:",
                    ref_exon = al_value["ref_exon"]
                    print ref_exon.ordinal
    
    pass
    
if __name__ == '__main__':
    main()