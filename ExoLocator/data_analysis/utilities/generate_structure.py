'''
Created on Apr 26, 2012

@author: marioot
'''
# Python imports
from itertools import chain

# utilities imports
from utilities                                      import FileUtilities
from utilities.Logger                               import Logger
from utilities.DescriptionParser                    import DescriptionParser
from utilities.DirectoryCrawler                     import DirectoryCrawler
from utilities.FileUtilities                        import check_status_file, check_status_file_no_alignment

# data analysis imports
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
from data_analysis.analysis.AlignmentPostprocessing import annotate_spurious_alignments_batch,\
     remove_overlapping_alignments_batch




def load_protein_configuration(protein_id, ref_species_dict = None):
    '''
    Loads the data from a description file, and calls the containers generating functions to create basic objects.
    @param protein_id: id of a single protein
    '''
    if not check_status_file_no_alignment(protein_id):
        return False
    
    if ref_species_dict is None:
        ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    logger                  = Logger.Instance()
    containers_logger       = logger.get_logger('containers')
    
    data_map_container      = DataMapContainer.Instance()
    protein_container       = ProteinContainer.Instance()
    gene_container          = GeneContainer.Instance()
    transcript_container    = TranscriptContainer.Instance()
    ens_exon_container      = EnsemblExonContainer.Instance()
    
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

        # data map
        data_map_key    = (protein_id, species_name)
        data_map        = DataMap(spec_protein_id, spec_transcript_id, 
                                  spec_gene_id, data_map_key, location_type, 
                                  location_id, strand, seq_begin, seq_end, ab_initio)
        try:
            data_map_container.add(data_map_key, data_map)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding to datamap".format(protein_id, species_name, e.args[0]))
        
        # everything else - protein, transcript, gene, ensembl exons
        protein         = Protein(spec_protein_id, data_map_key, ref_species_dict[species_name])
        gene            = Gene(spec_gene_id, data_map_key, ref_species_dict[species_name])
        transcript      = Transcript(spec_transcript_id, data_map_key, ref_species_dict[species_name])
        ens_exons       = EnsemblExons(data_map_key, ref_species_dict[species_name])
        try:
            ens_exons.load_exons()
        except (Exception), e:
            containers_logger.error("{0}, {1}, {2}, error loading exons".format(protein_id, species_name, e.args[0]))
        
        # add ensembl exons separately to the ensembl exon container
        for exon in ens_exons.exons.values():
            try:
                ens_exon_container.add (exon.exon_id, exon)
            except (KeyError, TypeError), e:
                containers_logger.error("{0}, {1}, {2}, error adding ensembl exon {3}".format(protein_id, species_name, e.args[0], exon.exon_id))
            
        # add everything to its container
        try:
            protein_container.   add(protein.protein_id, protein)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding protein".format(protein_id, species_name, e.args[0]))
        try:
            gene_container.      add(gene.gene_id, gene)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding gene".format(protein_id, species_name, e.args[0]))
        try:
            transcript_container.add(transcript.transcript_id, transcript)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding transcript".format(protein_id, species_name, e.args[0]))
            
        
    
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
        

        data_map_key    = (protein_id, species_name)
        data_map        = DataMap(spec_protein_id, "", "", data_map_key, location_type, location_id, strand, seq_begin, seq_end, ab_initio)
        try:
            data_map_container.add(data_map_key, data_map)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding to datamap".format(protein_id, species_name, e.args[0]))
        
        
        try:
            protein_container.add(protein.protein_id, protein)
        except (KeyError, TypeError), e:
            containers_logger.error("{0}, {1}, {2}, error adding protein".format(protein_id, species_name, e.args[0]))
        

    return True
  

def load_exon_configuration (ref_protein_id, ref_species_dict, exon_type):
    '''
    Load exons of a particular type for all available species
    @param ref_protein_id: referent protein id
    @param exon_type: exon_type: ensembl, genewise, blatn, tblastn, sw_gene, sw_exon
    '''
    
    dc                  = DescriptionParser()
    exon_container      = ExonContainer.Instance()
    
    logger              = Logger.Instance()
    containers_logger   = logger.get_logger('containers')
    
    if exon_type == "ensembl" or exon_type == "genewise":
        if not check_status_file_no_alignment(ref_protein_id):
            containers_logger.info ("{0},exon_type:{1},check status file -> failed".format(ref_protein_id, exon_type))
            return False
    else:
        if not check_status_file(ref_protein_id):
            containers_logger.info ("{0},exon_type:{1},check status file -> failed".format(ref_protein_id, exon_type))
            return False
    
    if not ref_species_dict:
        ref_species_dict = FileUtilities.get_reference_species_dictionary()

    (known_species, abinitio_species) = dc.get_separated_species(ref_protein_id)
    
    for species in known_species:
         
        ref_species = ref_species_dict[species]
        if exon_type != "genewise":
            if exon_type == "ensembl":
                exons = EnsemblExons ((ref_protein_id, species), ref_species)
                try:
                    exon_dict = exons.load_exons()
                except Exception, e:
                    containers_logger.error("{0},{1},{2},error loading exons".format(ref_protein_id, species, exon_type))
                    continue
            else:
                exons = Exons((ref_protein_id, species), ref_species, exon_type)
            try:
                exon_dict = exons.load_exons()
            except Exception, e:
                    containers_logger.error("{0},{1},{2},error loading exons".format(ref_protein_id, species, exon_type))
                    continue
            if not exon_dict:
                continue
        
            if (exon_type != "ensembl"):
                exons.set_exon_ordinals()
            data_map_key = [ref_protein_id, species]
            exon_container.add(exon_type, data_map_key, exons)
        
    for species in abinitio_species:
         
        ref_species = ref_species_dict [species]
        if exon_type != "ensembl":
            if exon_type == "genewise":
                exons = GenewiseExons ((ref_protein_id, species), ref_species)
            else:
                exons = Exons((ref_protein_id, species), ref_species, exon_type)
            try:
                exon_dict = exons.load_exons()
            except Exception, e:
                    containers_logger.error("{0},{1},{2},error loading exons".format(ref_protein_id, species, exon_type))
                    continue
            if not exon_dict:
                continue
            
            if exon_type != "genewise":
                exons.set_exon_ordinals()
            data_map_key = [ref_protein_id, species]
            try:
                exon_container.add(exon_type, data_map_key, exons)
            except Exception, e:
                    containers_logger.error("{0},{1},{2},error adding exons".format(ref_protein_id, species, exon_type))



def load_protein_configuration_batch(protein_id_list):
    '''
    Loads data from .descr files from all the proteins in the protein list
    @param protein_id_list: list of protein id's
    '''
    ref_species_dict    = FileUtilities.get_reference_species_dictionary()
    
    folders_loaded_cnt  = 0
    for protein_id in protein_id_list:
        if load_protein_configuration(protein_id, ref_species_dict) == True:
            folders_loaded_cnt += 1
    return folders_loaded_cnt


def load_exon_configuration_batch(protein_id_list, alignment_type):
    '''
    Loads exons for all the proteins in the protein list 
    for a particular alignment type
    '''
    ref_species_dict = FileUtilities.get_reference_species_dictionary()
   
    folders_loaded_cnt  = 0
    for protein_id in protein_id_list:
        if load_exon_configuration(protein_id, ref_species_dict, alignment_type) == True:
            folders_loaded_cnt += 1
    return folders_loaded_cnt
    

def fill_all_containers (load_alignments):
    '''
    Fills all the containers with correspondent data.
    The containers are: data maps, proteins, genes, transcripts, ensembl exons, and all the alignment exons
    '''
    dc = DirectoryCrawler()
    
    protein_list_raw = FileUtilities.get_protein_list()
    # flatten the raw protein list and take every second element, which is a protein id
    protein_list = list(chain.from_iterable(protein_list_raw))[0::2]
    algorithms = ["blastn", "tblastn", "sw_gene", "sw_exon"]
    for protein_id in protein_list:
        dc.generate_directory_tree(protein_id)
        
        
    ens_exon_container = load_protein_configuration_batch(protein_list)
    if ens_exon_container:
        
        load_exon_configuration_batch(protein_list, "ensembl")
        load_exon_configuration_batch(protein_list, "genewise")
        if load_alignments:
            load_exon_configuration_batch (protein_list, "blastn")
            load_exon_configuration_batch(protein_list, "tblastn")
            load_exon_configuration_batch(protein_list, "sw_gene")
            load_exon_configuration_batch(protein_list, "sw_exon") 
            remove_overlapping_alignments_batch(protein_list, ["blastn", "tblastn"])
            annotate_spurious_alignments_batch(protein_list, algorithms)
            

def main ():
    
    # fill all the data containers
    fill_all_containers(True)
    
    ec = ExonContainer.Instance()
    
    
    protein_list_raw = FileUtilities.get_protein_list()
    exon_number = {}
    for prot_id, exon_num in protein_list_raw:
        exon_number[prot_id] = int(exon_num)
    protein_list = list(chain.from_iterable(protein_list_raw))[0::2]
    '''
    # list of alignment algorithms    
    algorithms = ["blastn", "tblastn", "sw_gene", "sw_exon"]
    # POSTPROCESSING
    
    
    algorithms.append("genewise")
    
    # translate the alignment produced exons and populate the TranslatedProteinContainer
    #translate_alignment_exons_for_protein(protein_list, exon_number)
    
    ec = ExonContainer.Instance()
    for prot_id in protein_list:
        for spec in get_default_species_list():
            for alg in algorithms:
                try:
                    print alg, spec
                    exons = ec.get((prot_id, spec, alg))
                    ordered = exons.get_ordered_exons()
                    print alg
                    for exon in ordered:
                        print (exon.ordinal, exon.alignment_ordinal), exon.viability
                    
                except Exception:
                    pass
        
    exons = ec.get(("ENSP00000304822", "Equus_caballus", "sw_gene"))

    #new_exons = remove_UTR_ensembl_exons("ENSP00000253108", "Homo_sapiens", exons)
    #rint len(new_exons)
    #print exons.get_coding_cDNA()
    
    #create_statistics(protein_list)
   
    tp = TranslatedProteinContainer.Instance()
    ec = ExonContainer.Instance()
    #print tp.get("ENSP00000366061", "Ailuropoda_melanoleuca").get_sequence_record()

    #translate_ensembl_exons(protein_list)
    '''
    
    #translate_alignment_exons_for_protein(protein_list, exon_number)
    
    exons = ec.get(("ENSP00000293201", "Homo_sapiens", "ensembl"))
    coding_exons = exons.get_coding_exons()
    
    #for ce in coding_exons:
    #    print ce.start, ce.stop
    #    print ce.sequence
    
if __name__ == '__main__':
    main()