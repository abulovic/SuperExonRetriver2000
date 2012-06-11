'''
Created on May 6, 2012

@author: intern
'''

from timeit import itertools

from utilities                                      import FileUtilities
from utilities.DescriptionParser                    import DescriptionParser

from data_analysis.containers.ExonContainer         import ExonContainer
from utilities.Logger import Logger
from utilities.FileUtilities import check_status_file


def _calculate_total_score(valid_id_list, alignment_exons):
    '''
    Calculates the total score of ordered alignments.
    There may be multiple ordered subsets and this funcion serves to find the best subset.
    Score is determined by adding individual scores of exons (@Exon), which will
    favor longer alignments, and subsets with more exons
    '''
    score = 0.
    for exon in alignment_exons:
        if (exon.ordinal, exon.alignment_ordinal) in valid_id_list:
            score += exon.get_fitness()
    return score

def _is_sorted(exon_ordinal_list, strand):
    '''
    Checks whether list of exon ordinals is sorted
    '''
    sorted_list = sorted(exon_ordinal_list)
    if list(exon_ordinal_list) == sorted_list:
        return True
    return False


def _find_best_orderred_subset(alignment_exons, reference_exons, strand):
    
    full_array          = []
    highest_score       = 0.
    best_combination    = []
    
    # add all the exon ordinals to a list of all ordinals
    for exon in alignment_exons:
        full_array.append ((exon.ordinal, exon.alignment_ordinal))
        
    
    for i in range (len(full_array)):
        
        # find all combinations of ordinals of all possible lengths
        # ranging from 1 to number of exons
        all_combinations = itertools.combinations(full_array, len(full_array)-i)
        
        for comb in all_combinations:
            if _is_sorted (comb, strand):
                score = _calculate_total_score(comb, alignment_exons)
                if score > highest_score:
                    highest_score = score
                    best_combination = comb
        if highest_score > 0.:
            return best_combination
    
    
            
def annotate_spurious_alignments(exons_key):
    '''
    Annotates all the alignments which are not in the correct order.
    Annotation means their viability variable will be set to False.
    (Supporting the assumption that all exons are in the correct, sequential order)
    @param exons_key: (reference protein id, species)
    @param alignment_type: blastn, tblastn, sw_gene, sw_exon
    '''
    
    (ref_protein_id, 
     species, 
     alignment_type)            = exons_key
     
    print "Annotating spurious alignments %s,%s,%s" % (ref_protein_id, species, alignment_type)
     

    if not check_status_file(ref_protein_id):
        return None
     
    exon_container              = ExonContainer.Instance()
    reference_species_dict      = FileUtilities.get_reference_species_dictionary()
    
    # load logging utilities
    logger                      = Logger.Instance()
    containers_logger           = logger.get_logger("containers")
    
    # get the reference exons: (ref_prot_id, ref_species, ensembl)
    reference_exons     = exon_container.get((ref_protein_id, 
                                              reference_species_dict[species], 
                                              "ensembl"))
    # try to get the exons which are the product of specified alignment
    try:
        alignment_exons = exon_container.get((ref_protein_id, species, alignment_type))
    except KeyError:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    
    exon_list = []
    # flatten the exon list (there may be multiple exons for one reference exon)
    for al_exons in alignment_exons.alignment_exons.values():
        for al_exon in al_exons:
            exon_list.append(al_exon)
    
    # dependant on the strand, calculate the ordering of exons
    strands = DescriptionParser().get_strand_information(ref_protein_id)
    strand  = strands[species]
    
    if strand == 1:
        reverse_flag = False
    else:
        reverse_flag = True

    correct_order_exons = _find_best_orderred_subset (alignment_exons.get_ordered_exons(),
                                                      reference_exons,
                                                      strand)
    # set viability for all the exons
    # if exon is not in the correct order, set viability to False, True otherwise
    for exon in exon_list:
        ordinal = (exon.ordinal, exon.alignment_ordinal)
        if ordinal in correct_order_exons:
            for al_exon in alignment_exons.alignment_exons[exon.ref_exon_id]:
                if al_exon.alignment_info["query_start"] == exon.alignment_info["query_start"]:
                    al_exon.set_viability(True)
        else:
            for al_exon in alignment_exons.alignment_exons[exon.ref_exon_id]:
                if al_exon.alignment_info["query_start"] == exon.alignment_info["query_start"]:
                    if alignment_type == "sw_exon":
                        al_exon.set_viability(True)
                    else:
                        #print "Exon %f not viable: (%s,%s,%s)" % (al_exon.ordinal, ref_protein_id, species, alignment_type)
                        al_exon.set_viability(False)
                        
    exon_container.update(exons_key, alignment_exons)
                    
                    
def annotate_spurious_alignments_batch (protein_list, algorithms):
    
    
    for protein_id in protein_list :
        species_list = DescriptionParser().get_species(protein_id)
        
        for species in species_list:
        
            for alg in algorithms:
                
                exon_key = (protein_id, species, alg)
                #remove_overlapping_alignments(exon_key)
                annotate_spurious_alignments(exon_key)
                

def remove_overlapping_alignments_batch (protein_list, algorithms):
    for protein_id in protein_list :
        species_list = DescriptionParser().get_species(protein_id)
        
        for species in species_list:
        
            for alg in algorithms:
                
                exon_key = (protein_id, species, alg)
                #remove_overlapping_alignments(exon_key)
                remove_overlapping_alignments(exon_key)
    

def remove_overlapping_alignments (exons_key):
    (ref_protein_id, 
     species, 
     alignment_type)            = exons_key
    printin = False
    if printin: 
        print "Removing blastn overlaps (%s,%s,%s)..." % (ref_protein_id, species, alignment_type)

    if not check_status_file(ref_protein_id):
        return None
    
    exon_container              = ExonContainer.Instance()
    reference_species_dict      = FileUtilities.get_reference_species_dictionary()
    
    # load logging utilities
    logger                      = Logger.Instance()
    containers_logger           = logger.get_logger("containers")
    
    # get the reference exons: (ref_prot_id, ref_species, ensembl)
    reference_exons     = exon_container.get((ref_protein_id, 
                                              reference_species_dict[species], 
                                              "ensembl"))
    # try to get the exons which are the product of specified alignment
    try:
        alignment_exons = exon_container.get(exons_key)
    except KeyError:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    
    for ref_exon_id in alignment_exons.alignment_exons:
        al_exons = alignment_exons.alignment_exons[ref_exon_id]
        if printin:
            print ref_exon_id
        toplevel_start = 0
        toplevel_stop = 0
        #for al_exon in sorted(al_exons, key = lambda al_exon: al_exon.get_fitness(), reverse = True):
        for al_exon in al_exons:
            
            exon_start = al_exon.alignment_info["sbjct_start"]
            exon_stop = exon_start + al_exon.alignment_info["length"]
            
            if not al_exon.viability:
                if printin:
                    print "Exon %d not viable" % al_exon.ordinal
                continue
            
            if not toplevel_start:
                # if toplevel locations haven't been set, set them
                toplevel_start = exon_start
                toplevel_stop = exon_stop
                toplevel_exon = al_exon
                al_exon.set_viability(True)
                if printin:
                    print "First exon: %d - %d" % (exon_start, exon_stop)
                
            elif exon_start < toplevel_start and exon_stop > toplevel_stop:
                toplevel_exon.set_viability(False)
                toplevel_exon = al_exon
                toplevel_start = exon_start
                toplevel_stop = exon_stop
                al_exon.set_viability(True)
                if printin:
                    print "  New toplevel: %d - %d" % (exon_start, exon_stop)
                
            else:
                # what this wonderful if checks if one of the following cases:
                # if the exon is contained within the toplevel exon
                #          |----------------------|
                #               |------|
                # or the start is to the left of the toplevel, but they are still overlapping
                #      |----------|
                # or the end is to the right of the toplevel, but they are still overlapping
                #                        |--------------|
                if (exon_start >=toplevel_start and exon_stop <= toplevel_stop) or \
                (exon_start <= toplevel_start and (exon_stop >= toplevel_start and exon_stop <= toplevel_stop)) or \
                ((exon_start >= toplevel_start and exon_start <= toplevel_stop) and exon_stop >= toplevel_stop):
                    if printin:
                        print "   Bad exon: %d - %d" % (exon_start, exon_stop)
                    al_exon.set_viability(False)
                else:
                    if exon_start < toplevel_start:
                        toplevel_start = exon_start
                    if exon_stop > toplevel_stop:
                        toplevel_stop = exon_stop
                    if printin:
                        print "  Good exon: %d - %d" % (exon_start, exon_stop)
                        
    exon_container.update(exons_key, alignment_exons)
    
    
if __name__ == '__main__':
    ec = ExonContainer.Instance()
    

            
                    
                    
