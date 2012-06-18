'''
Created on May 2, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.containers.EnsemblExonContainer import EnsemblExonContainer
from utilities import FileUtilities
from utilities.DescriptionParser import DescriptionParser
import csv
from utilities.Logger import Logger
from data_analysis.analysis.AlignmentPostprocessing import remove_overlapping_alignments
from data_analysis.utilities.ExonUtils import remove_UTR_ensembl_exons





def produce_statistics_for_alignment (exons_key, alignment_type, reference_exons):
    '''
    Produces the straight-forward statistics for the alignment.
    For each alignment, it only calculates the coverage percentage. 
    Where there are multiple alignments for a particular exon, percentages are summed. 
    @param exons_key: (reference protein id, species)
    @param alignment_type: blastn, tblastn, sw_gene, sw_exon
    @return: list of similarity percentages in correct order starting from first reference exon to thel last
    '''
    print alignment_type, 
    (ref_protein_id, species) = exons_key
    print species
    
    # if the alignment type is tblastn, we have to multiply
    # the coverage by 3 because the length is expressed in AAs, not in NBs
    if alignment_type == "tblastn":
        coverage_constant = 3.
    else:
        coverage_constant = 1.
    
    exon_container          = ExonContainer.Instance()
    reference_species_dict  = FileUtilities.get_reference_species_dictionary()

    logger = Logger.Instance()
    containers_logger = logger.get_logger("containers")
    
    #reference_exons = exon_container.get((ref_protein_id, reference_species_dict[species], "ensembl"))
    try:
        alignment_exons = exon_container.get((ref_protein_id, species, alignment_type))
    except KeyError:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    perc_list = []
    
    # TODO
    
    remove_overlapping_alignments((ref_protein_id, species, alignment_type))
    for ref_exon in reference_exons:
        ref_exon_id = ref_exon.exon_id
        if ref_exon_id not in alignment_exons.alignment_exons:
            perc_list.append(0)
            #print "%s length: %d\n\tnot present in alignment" % (ref_exon_id, len(ref_exon.sequence))
        else:
            al_exons = alignment_exons.alignment_exons[ref_exon_id]
            internal_stat = 0.
            print ref_exon_id, "length: %d" % len(ref_exon.sequence)
            for al_exon in al_exons:
                if al_exon.viability:
                    #print al_exon.alignment_info["sbjct_start"], al_exon.alignment_info["sbjct_end"]
                    
                    internal_stat += coverage_constant * float(al_exon.alignment_info["identities"]) / len(ref_exon.sequence)
                    if internal_stat > 1:
                        print "Coverage cannot be larger than 1 (%s,%s,%s)" % (ref_protein_id, species, alignment_type)
                        print al_exon.alignment_info["sbjct_start"], len(ref_exon.sequence)
                        raise ValueError ("Coverage cannot be larger than 1 (%s,%s,%s,%s)" % (ref_protein_id, species, alignment_type,ref_exon_id))
                #print "\t%1.2f" % ( float(al_exon.alignment_info["length"] - al_exon.alignment_info["gaps"]) / len(ref_exon.sequence))
            perc_list.append(internal_stat)
    return perc_list
                

def create_protein_statistics (protein_id, statistics_file):
    '''
    Creates protein statistics for a protein
    It works only on exons which get translated to protein. The UTR regions are cropped.
    '''
    
    print "Creating protein statistics for %s, writing to %s" % (protein_id, statistics_file)
    species_list = DescriptionParser().get_species(protein_id)

    exon_container          = ExonContainer.Instance()
    ensembl_exon_containter = EnsemblExonContainer.Instance()
    
    header_row = ["Species name", "Alignment"]
    
    reference_exons = exon_container.get((protein_id, "Homo_sapiens", "ensembl"))
    reference_coding_exons = reference_exons.get_coding_exons()
    for ref_exon in reference_coding_exons:
        header_row.append(ref_exon.exon_id)

    csv_writer = csv.writer(open(statistics_file, 'wb'), delimiter=",")
    csv_writer.writerow(header_row)
    alignments = ["blastn", "tblastn", "sw_gene", "sw_exon"]
    
    for species in species_list :
        for al in alignments:
            new_row = [species, al]
            stats = produce_statistics_for_alignment((protein_id, species), al, reference_coding_exons)
            if stats:
                new_row.extend (stats)
            else:
                continue
            csv_writer.writerow(new_row)


        
