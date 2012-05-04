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
from timeit import itertools


def _calculate_total_score(valid_id_list, alignment_exons):
    score = 0.
    for exon in alignment_exons:
        if exon.ordinal in valid_id_list:
            score += exon.get_fitness()
    return score

def _is_sorted(in_list):
    for i in range (len(in_list)-1):
        if in_list[i] > in_list[i+1]:
            return False
    return True


def _find_best_orderred_subset(alignment_exons, reference_exons, strand):
    
    full_array          = []
    highest_score       = 0.
    best_combination    = []
    
    for exon in alignment_exons:
        #ref_exon = reference_exons.exons[exon.ref_exon_id]
        full_array.append (exon.ordinal)
    print full_array
    for i in range (len(full_array)):
        all_combinations = itertools.combinations(full_array, len(full_array)-i)
        for comb in all_combinations:
            if _is_sorted (comb):
                score = _calculate_total_score(comb, alignment_exons)
                if score > highest_score:
                    highest_score = score
                    best_combination = comb
        if highest_score > 0.:
            print best_combination
            return best_combination
    
    
            
def remove_spurious_alignments(exons_key, alignment_type):
    '''
    Removes all the alignments which are not in the correct order.
    (Supporting the assumption that all exons are in the correct, sequential order)
    @param exons_key: (reference protein id, species)
    @param alignment_type: blastn, tblastn, sw_gene, sw_exon
    '''
    (ref_protein_id, species)   = exons_key
    exon_container              = ExonContainer.Instance()
    reference_species_dict      = FileUtilities.get_reference_species_dictionary()
    
    logger              = Logger.Instance()
    containers_logger   = logger.get_logger("containers")
    
    # get the reference exons: (ref_prot_id, ref_species, ensembl)
    reference_exons     = exon_container.get((ref_protein_id, reference_species_dict[species], "ensembl"))
    # try to get the exons which are the product of specified alignment
    try:
        alignment_exons = exon_container.get((ref_protein_id, species, alignment_type))
    except KeyError:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    
    exon_list = []
    # flatten the exons (there may be multiple exons for one reference exon)
    for al_exons in alignment_exons.alignment_exons.values():
        for al_exon in al_exons:
            exon_list.append(al_exon)
            
    strands = DescriptionParser().get_strand_information(ref_protein_id)
    strand = strands[species]
    
    if strand == 1:
        correct_order_exons = _find_best_orderred_subset (sorted(exon_list, 
                                                          key = lambda al_exon: (al_exon.alignment_info["query_start"])), 
                                                          reference_exons,
                                                          strand)
    else:
        correct_order_exons = _find_best_orderred_subset (sorted(exon_list, 
                                                                 key = lambda al_exon: (al_exon.alignment_info["query_start"]),
                                                                 reverse = True), 
                                                          reference_exons,
                                                          strand)
    # set viability for all the exons
    # if exon is not in the correct order exons, set viability to False, True otherwise
    for exon in exon_list:
        if exon.ordinal in correct_order_exons:
            for al_exon in alignment_exons.alignment_exons[exon.ref_exon_id]:
                if al_exon.alignment_info["query_start"] == exon.alignment_info["query_start"]:
                    al_exon.set_viability(True)
        else:
            for al_exon in alignment_exons.alignment_exons[exon.ref_exon_id]:
                if al_exon.alignment_info["query_start"] == exon.alignment_info["query_start"]:
                    al_exon.set_viability(False)
    

def produce_statistics_for_alignment (exons_key, alignment_type):
    '''
    Produces the straight-forward statistics for the alignment.
    For each alignment, it only calculates the coverage percentage. 
    Where there are multiple alignments for a particular exon, percentages are summed. 
    @param exons_key: (reference protein id, species)
    @param alignment_type: blastn, tblastn, sw_gene, sw_exon
    @return: list of similarity percentages in correct order starting from first reference exon to thel last
    '''
    (ref_protein_id, species) = exons_key
    
    exon_container          = ExonContainer.Instance()
    ensembl_exon_containter = EnsemblExonContainer.Instance()
    reference_species_dict  = FileUtilities.get_reference_species_dictionary()

    logger = Logger.Instance()
    containers_logger = logger.get_logger("containers")
    
    reference_exons = exon_container.get((ref_protein_id, reference_species_dict[species], "ensembl"))
    try:
        alignment_exons = exon_container.get((ref_protein_id, species, alignment_type))
    except KeyError:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    perc_list = []
    for ref_exon in reference_exons.get_ordered_exons():
        ref_exon_id = ref_exon.exon_id
        if ref_exon_id not in alignment_exons.alignment_exons:
            perc_list.append(0)
            print "%s length: %d\n\tnot present in alignment" % (ref_exon_id, len(ref_exon.sequence))
        else:
            al_exons = alignment_exons.alignment_exons[ref_exon_id]
            internal_stat = 0.
            print ref_exon_id, "length: %d" % len(ref_exon.sequence)
            for al_exon in al_exons:
                internal_stat += float(al_exon.alignment_info["length"] - al_exon.alignment_info["gaps"]) / len(ref_exon.sequence)
                print "\t%1.2f" % ( float(al_exon.alignment_info["length"] - al_exon.alignment_info["gaps"]) / len(ref_exon.sequence))
            perc_list.append(internal_stat)
    return perc_list
                

def create_protein_statistics (protein_id, statistics_file):
    
    species_list = DescriptionParser().get_species(protein_id)

    exon_container          = ExonContainer.Instance()
    ensembl_exon_containter = EnsemblExonContainer.Instance()
    
    header_row = ["Species name", "Alignment"]
    
    reference_exons = exon_container.get((protein_id, "Homo_sapiens", "ensembl"))
    for ref_exon in reference_exons.get_ordered_exons():
        header_row.append(ref_exon.exon_id)

    csv_writer = csv.writer(open(statistics_file, 'wb'), delimiter=",")
    csv_writer.writerow(header_row)
    alignments = ["blastn", "tblastn", "sw_gene", "sw_exon"]
    
    for species in species_list :
        for al in alignments:
            new_row = [species, al]
            stats = produce_statistics_for_alignment((protein_id, species), al)
            if stats:
                new_row.extend (stats)
            else:
                continue
            csv_writer.writerow(new_row)
            
    
if __name__ == '__main__':
    #create_protein_statistics("ENSP00000341765", "/home/intern/test_stats.csv")
    remove_spurious_alignments(("ENSP00000311134", "Ailuropoda_melanoleuca"), "sw_gene")

        
