'''
Created on Apr 26, 2012

@author: marioot
'''
import os
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
from data_analysis.analysis.AlignmentPostprocessing import annotate_spurious_alignments_batch,\
    annotate_spurious_alignments, remove_overlapping_alignments_batch
from utilities.FileUtilities                        import check_status_file
from utilities.ConfigurationReader                  import ConfigurationReader
from data_analysis.analysis.AlignmentStatistics     import create_protein_statistics
from data_analysis.analysis.Exon_translation import Exon_translation
from data_analysis.analysis.TranscriptionMachinery import transcribe_exons,\
   chop_off_start_utr, chop_off_end_utr
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator



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
            data_map     = DataMap(spec_protein_id, spec_transcript_id, spec_gene_id, data_map_key, strand, seq_begin, seq_end, ab_initio)
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
            data_map = DataMap(spec_protein_id, "", "", data_map_key, strand, seq_begin, seq_end, ab_initio)
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
    dc = DirectoryCrawler()
    acg = AlignmentCommandGenerator()
    logger = Logger.Instance()
    data_logger = logger.get_logger("data_loading")
    
    protein_list_raw = FileUtilities.get_protein_list()
    protein_list = []
    exon_number = {}
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
        exon_number[protein_tuple[0]] = int(protein_tuple[1])
        dc.generate_directory_tree(protein_tuple[0])
        
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
    remove_overlapping_alignments_batch(protein_list, algorithms)
    #annotate_spurious_alignments(("ENSP00000311134", "Sus_scrofa", "sw_gene"))
    #transcribe_aligned_exons (protein_list, algorithms)
    
    
        
    
    dmc = DataMapContainer.Instance()
    pc  = ProteinContainer.Instance()
    gc  = GeneContainer.Instance()
    tc  = TranscriptContainer.Instance()
    eec = EnsemblExonContainer.Instance()
    ec  = ExonContainer.Instance()
    

    algorithms = ["sw_gene"]        
    for protein_id in protein_list :
        
        if not check_status_file(protein_id):
            continue
        species_list = DescriptionParser().get_species(protein_id)
        print species_list
        #species_list = ["Cavia_porcellus"]
        #species_list = ["Homo_sapiens"]
        for species in species_list:
            
            for alg in algorithms:
                
                try:
                    exon_key = (protein_id, species, alg)
                    target_prot = pc.get(protein_id)
                    target_prot_seq = target_prot.get_sequence_record().seq
                    
                    print protein_id, species
                    
                    try:
                        exons = ec.get(exon_key)
                    except KeyError:
                        print "NONONO"
                        continue
                    exons_for_transcription = []
                    
                    test_next_exon = False
                    coding_dna_started = False
                    last_exon = False
                    for al_exon in exons.get_ordered_exons():
                        al_exon = al_exon[0]
                        ref_exon = eec.get(al_exon.ref_exon_id)
                        trans_exon = Exon_translation(al_exon.ordinal, 
                                                      ref_exon.length, 
                                                      al_exon.alignment_info["query_seq"], 
                                                      al_exon.alignment_info["sbjct_seq"])
                        trans_exon.set_intervals(al_exon.alignment_info["query_start"], 
                                                 al_exon.alignment_info["query_end"], 
                                                 al_exon.alignment_info["sbjct_start"], 
                                                 al_exon.alignment_info["sbjct_end"])
                        trans_exon.set_identity(al_exon.alignment_info["identities"], al_exon.alignment_info["length"])
                        trans_exon.set_viablity(al_exon.viability)
    
                        if last_exon:
                            trans_exon.viability = False
                            
                        if trans_exon.viability:
                            (trans_exon, last_exon) = chop_off_start_utr(al_exon.ref_exon_id, trans_exon, target_prot_seq, exon_number[protein_id])
                            trans_exon = chop_off_end_utr (al_exon.ref_exon_id, trans_exon, target_prot_seq)
                        
                        exons_for_transcription.append(trans_exon)
                        
                        
                    sequences_for_fasta = [target_prot.get_sequence_record()]    
                        
                    
                    print "Target prot:", target_prot_seq
                    resulting_protein_seq = transcribe_exons(exons_for_transcription, target_prot_seq)
                    print "Query prot: ", resulting_protein_seq
                    seq_to_write = SeqRecord(Seq(resulting_protein_seq), species, description = "")
                    fasta = "%s/%s.fa" % (dc.get_assembled_protein_path(protein_id), species)
                    fasta_file = open(fasta, "w")
                    SeqIO.write(seq_to_write, fasta_file, "fasta")
                    fasta_file.close()
                    data_map = dmc.get((protein_id, species))
                    species_protein = pc.get(data_map.protein_id)
                    print "Ensem prot: ", species_protein.get_sequence_record().seq, "\n"
                    
                    sequences_for_fasta.append(seq_to_write)
                    ensembl_protein = species_protein.get_sequence_record()
                    sequences_for_fasta.append(ensembl_protein)
                    
                    msa_fasta = "%s/%s.fa" % (dc.get_mafft_path(protein_id), species)
                    msa_afa = "%s/%s.afa" % (dc.get_mafft_path(protein_id), species)
                    msa_fasta_file = open(msa_fasta, "w")
                    SeqIO.write(sequences_for_fasta, msa_fasta_file, "fasta")
                    msa_fasta_file.close()
                    
                    mafft_cmd = acg.generate_mafft_command(msa_fasta, msa_afa)
                    os.system(mafft_cmd)
                    
                except Exception, e:
                    data_logger.error("%s,%s,%s" % (protein_id, species, e))
    
    
if __name__ == '__main__':
    main()