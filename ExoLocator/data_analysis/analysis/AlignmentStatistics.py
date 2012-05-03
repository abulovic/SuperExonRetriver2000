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

def produce_statistics_for_alignment (exons_key, alignment_type):
    (ref_protein_id, species) = exons_key
    
    exon_container          = ExonContainer.Instance()
    ensembl_exon_containter = EnsemblExonContainer.Instance()
    reference_species_dict  = FileUtilities.get_reference_species_dictionary()
    
    logger = Logger.Instance()
    containers_logger = logger.get_logger("containers")
    
    reference_exons = exon_container.get((ref_protein_id, reference_species_dict[species], "ensembl"))
    try:
        alignment_exons = exon_container.get((ref_protein_id, species, alignment_type))
    except KeyError, e:
        containers_logger.error ("{0},{1},{2}".format(ref_protein_id, species, alignment_type))
        return None
    perc_list = []
    for ref_exon in reference_exons.get_ordered_exons():
        ref_exon_id = ref_exon.exon_id
        if ref_exon_id not in alignment_exons.alignment_exons:
            perc_list.append("")
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
    create_protein_statistics("ENSP00000341765", "/home/intern/test_stats.csv")

        
