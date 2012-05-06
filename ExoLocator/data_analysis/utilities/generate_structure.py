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
from data_analysis.containers.ExonContainer         import ExonContainer
from data_analysis.base.GenewiseExons               import GenewiseExons
from data_analysis.analysis.AlignmentStatistics     import create_protein_statistics,\
    produce_statistics_for_alignment
from data_analysis.analysis.TranscriptionMachinery import transcribe_aligned_exons
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from data_analysis.analysis.AlignmentPostprocessing import annotate_spurious_alignments_batch
from utilities.FileUtilities import check_status_file
from utilities.ConfigurationReader import ConfigurationReader



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
            
            for exon in ens_exons.exons.values():
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
  

def load_exon_configuration (ref_protein_id, ref_species_dict, alignment_type):
    
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
        if alignment_type != "genewise":
            if alignment_type == "ensembl":
                exons = EnsemblExons ((ref_protein_id, species_name), ref_species)
                exon_dict = exons.load_exons()
            else:
                exons = Exons((ref_protein_id, species_name), ref_species, alignment_type)
            exon_dict = exons.load_exons()
            if not exon_dict:
                continue
        
            if (alignment_type != "ensembl"):
                exons.set_exon_ordinals()
            data_map_key = [ref_protein_id, species_name]
            exon_container.add(alignment_type, data_map_key, exons)
        
    for species_data in abinitio_proteins:
        (species_name,
         spec_protein_id,
         location_type,
         assembly,
         location_id,
         seq_begin,
         seq_end,
         strand) = species_data
         
        ref_species = ref_species_dict [species_name]
        if alignment_type != "ensembl":
            if alignment_type == "genewise":
                exons = GenewiseExons ((ref_protein_id, species_name), ref_species)
            else:
                exons = Exons((ref_protein_id, species_name), ref_species, alignment_type)
            exon_dict = exons.load_exons()
            if not exon_dict:
                continue
            
            if alignment_type != "genewise":
                exons.set_exon_ordinals()
            data_map_key = [ref_protein_id, species_name]
            exon_container.add(alignment_type, data_map_key, exons)



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


def load_exon_configuration_batch(protein_id_list, alignment_type):
    
    ref_species_dict = FileUtilities.get_reference_species_dictionary()
    ens_exon_container = EnsemblExonContainer.Instance()
   
    folders_loaded_cnt  = 0
    for protein_id in protein_id_list:
        if load_exon_configuration(protein_id, ref_species_dict, alignment_type) == True:
            folders_loaded_cnt += 1
    return folders_loaded_cnt
    





def main ():
    
    cr = ConfigurationReader.Instance()
    
    protein_list_raw = FileUtilities.get_protein_list()
    protein_list = []
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
        
    ens_exon_container = load_protein_configuration_batch(protein_list)
    if ens_exon_container:
        
        load_exon_configuration_batch(protein_list, "ensembl")
        load_exon_configuration_batch(protein_list, "genewise")
        load_exon_configuration_batch (protein_list, "blastn")
        load_exon_configuration_batch(protein_list, "tblastn")
        load_exon_configuration_batch(protein_list, "sw_gene")
        load_exon_configuration_batch(protein_list, "sw_exon")
        
    algorithms = ["blastn", "tblastn", "sw_gene", "sw_exon"]
    annotate_spurious_alignments_batch(protein_list, algorithms)
    #transcribe_aligned_exons (protein_list, algorithms)
    
        
    
    dmc = DataMapContainer.Instance()
    pc  = ProteinContainer.Instance()
    gc  = GeneContainer.Instance()
    tc  = TranscriptContainer.Instance()
    eec = EnsemblExonContainer.Instance()
    ec  = ExonContainer.Instance()
    
    produce_statistics_for_alignment( ("ENSP00000365108", "Homo_sapiens"), "sw_exon")    
    '''
    for protein_id in protein_list:
        if check_status_file(protein_id):
            statistics_file_path = "%s/%s/stats.csv" % ( cr.get_value("root", "session_dir"), protein_id)
            create_protein_statistics(protein_id, statistics_file_path)
            '''
    
if __name__ == '__main__':
    main()