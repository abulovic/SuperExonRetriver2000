'''
Created on Jun 24, 2012

@author: anana
'''
from utilities.FileUtilities import get_protein_list, get_species_list
from data_analysis.translation.BestProteinProduct import BestProteinProduct
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.translation.BestExonAlignmentContainer import BestExonAlignmentContainer
from data_analysis.utilities.generate_structure import fill_all_containers
from Bio.SeqRecord import SeqRecord
from utilities.DirectoryCrawler import DirectoryCrawler
from Bio import SeqIO

def main ():
    #ERROR FILE:::
    err_f = open('/home/marioot/Eclipse/err_status_monday.txt', 'w')

    fill_all_containers(True)
    
    protein_tuples = get_protein_list()
    ec = ExonContainer.Instance()
    beac = BestExonAlignmentContainer.Instance()
    dc = DirectoryCrawler()
    
    for (protein_id, exon_num) in protein_tuples:
        
        if int(exon_num) > 15:
            print "too big"
            continue
        
        species_list = get_species_list(protein_id, None)
        try:
            ref_exons = ec.get((protein_id, "Homo_sapiens", "ensembl"))
        except KeyError:
            print "ERROR: No protein %s" % protein_id
            continue
        
        for species in species_list:
            #try:
            print "\nBest_exon_al: %s, %s" % (protein_id, species)
            
            bpp = BestProteinProduct (protein_id, species, "Homo_sapiens")
            bpp.load_alignments()
            bpp.decide_on_best_exons()
            #bpp.patch_interexon_AAS()
            
            for ref_exon in ref_exons.get_coding_exons():
                
                best_exon_alignment = bpp.best_exons[ref_exon.exon_id]
                if best_exon_alignment:
                    beac.add(ref_exon.exon_id, species, best_exon_alignment)
                    print "%d. Exon status: %s (%s)" % (ref_exon.ordinal, best_exon_alignment.status, ref_exon.exon_id)
                    if best_exon_alignment.sw_gene_alignment:
                        print ref_exon.sequence[ref_exon.frame:].translate()
                        best_exon_alignment.sw_gene_alignment.create_cDNA()
                        print "\tAdded  %2d alignment pieces" % (len(best_exon_alignment.sw_gene_alignment.alignment_pieces))
                        for al_piece in best_exon_alignment.sw_gene_alignment.alignment_pieces:
                            print "\t\t%s:" % (al_piece.type),
                            if al_piece.type in ["coding", "insertion"]:
                                print "PROT: %d-%d, GENOME: %d-%d, %s" % (al_piece.ref_protein_start,
                                                                        al_piece.ref_protein_stop,
                                                                        al_piece.genomic_start, 
                                                                        al_piece.genomic_stop, 
                                                                        al_piece.sequence_id)
                                print "\t\t\tHUMAN:", al_piece.ref_protein_seq
                                print "\t\t\tSPEC :", al_piece.spec_protein_seq
                            else:
                                print
                                
            whole_prot =  bpp.get_spec_protein_translation()
            whole_prot_rec = SeqRecord(whole_prot, id = species, description = "assembled_protein")
            file_name = "%s/%s.fa" % (dc.get_assembled_protein_path(protein_id), species)
            SeqIO.write(whole_prot_rec, file_name, "fasta")
            #except Exception, e:
            #    print '{0} {1} \n'.format(protein_id, species)
            #    err_f.write('{0} {1} \n'.format(protein_id, species))
            


if __name__ == '__main__':
    main()
