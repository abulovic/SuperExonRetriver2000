'''
Created on Jun 25, 2012

@author: anana
'''
from utilities.DirectoryCrawler import DirectoryCrawler
from utilities.FileUtilities import get_protein_list, get_species_list,\
    check_status_file, write_seq_records_to_file, load_fasta_single_record
import os
from Bio import SeqIO
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer
from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator
from data_analysis.utilities.generate_structure import fill_all_containers
from Bio.Alphabet import IUPAC

def create_msa_alignments ():
    
    dc = DirectoryCrawler()
    pc = ProteinContainer.Instance()
    dmc = DataMapContainer.Instance()
    acg = AlignmentCommandGenerator()
    
    fill_all_containers(False)
    
    
    
    for (prot_id, exon_num) in get_protein_list():
        
        if not check_status_file(prot_id):
            continue
        
        ref_prot_rec = pc.get(prot_id).get_sequence_record()
        
        exoloc_proteins = []
        ensembl_proteins = []
        
        exoloc_proteins.append(ref_prot_rec)
        ensembl_proteins.append(ref_prot_rec)
        
        assembled_dir = dc.get_assembled_protein_path(prot_id)
        for fasta in sorted(os.listdir(assembled_dir)):
            if fasta == "Homo_sapiens.fa":
                continue
            abs_fasta = "%s/%s" % (assembled_dir, fasta)
            prot_rec = load_fasta_single_record(abs_fasta, IUPAC.protein)
            exoloc_proteins.append(prot_rec)
            
        species_list = get_species_list(prot_id, None)
        for species in species_list:
            if species == "Homo_sapiens":
                continue
            data_map = dmc.get((prot_id, species))
            prot_rec = pc.get(data_map.protein_id).get_sequence_record()
            prot_rec.id = species
            ensembl_proteins.append(prot_rec) 
            
        msa_exoloc_path = "%s/msa_exoloc.fa" % dc.get_mafft_path(prot_id)
        msa_ensembl_path = "%s/msa_ensembl.fa" % dc.get_mafft_path(prot_id)
        
        write_seq_records_to_file(msa_exoloc_path, exoloc_proteins)
        write_seq_records_to_file(msa_ensembl_path, ensembl_proteins)
        
        cmd = acg.generate_mafft_command(msa_exoloc_path, "%s/msa_exoloc.afa" % dc.get_mafft_path(prot_id))
        print cmd
        os.system(cmd)
        cmd = acg.generate_mafft_command(msa_ensembl_path, "%s/msa_ensembl.afa" % dc.get_mafft_path(prot_id))
        print  cmd
        os.system(cmd)
        
def create_species_msa_alignments ():
    
    dc = DirectoryCrawler()
    pc = ProteinContainer.Instance()
    dmc = DataMapContainer.Instance()
    acg = AlignmentCommandGenerator()
    
    fill_all_containers(False)
    

    for (prot_id, exon_num) in get_protein_list():
        
        if not check_status_file(prot_id):
            continue
        
        ref_prot_rec = pc.get(prot_id).get_sequence_record()
        ref_prot_rec.id = "Homo_sapiens"

        assembled_dir = dc.get_assembled_protein_path(prot_id)    
        species_list = get_species_list(prot_id, None)
        
        for species in species_list:
            
            protein_recs = []
            protein_recs.append(ref_prot_rec)
            
            if species == "Homo_sapiens":
                continue
            data_map = dmc.get((prot_id, species))
            prot_rec = pc.get(data_map.protein_id).get_sequence_record()
            prot_rec.id = species
            protein_recs.append(prot_rec)
            
            if "%s.fa" % species in os.listdir(assembled_dir):
                exoloc_protein_rec = load_fasta_single_record("%s/%s.fa" % (assembled_dir, species), IUPAC.protein)
                protein_recs.append(exoloc_protein_rec)

            
            msa_species_path = "%s/%s.fa" % (dc.get_mafft_path(prot_id), species)
            
            if len(protein_recs) == 1:
                continue
            write_seq_records_to_file(msa_species_path, protein_recs)
            
            cmd = acg.generate_mafft_command(msa_species_path, "%s/%s.afa" % (dc.get_mafft_path(prot_id), species))
            print cmd
            os.system(cmd)
            
            os.remove(msa_species_path)


if __name__ == '__main__':
    #create_msa_alignments()
    create_species_msa_alignments()
