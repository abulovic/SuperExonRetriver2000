'''
Created on Jun 25, 2012

@author: anana
'''
from utilities.DirectoryCrawler import DirectoryCrawler
from utilities.FileUtilities import get_protein_list, get_species_list,\
    check_status_file
import os
from Bio import SeqIO
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer
from pipeline.utilities.AlignmentCommandGenerator import AlignmentCommandGenerator
from data_analysis.utilities.generate_structure import fill_all_containers

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
            prot_rec = SeqIO.read(abs_fasta, "fasta")
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
        
        SeqIO.write(exoloc_proteins, msa_exoloc_path, "fasta")
        SeqIO.write(ensembl_proteins, msa_ensembl_path, "fasta")
        
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
                exoloc_protein_rec = SeqIO.read("%s/%s.fa" % (assembled_dir, species), "fasta")
                protein_recs.append(exoloc_protein_rec)

            
            msa_species_path = "%s/%s.fa" % (dc.get_mafft_path(prot_id), species)
            
            if len(protein_recs) == 1:
                continue
            SeqIO.write(protein_recs, msa_species_path, "fasta")
            
            cmd = acg.generate_mafft_command(msa_species_path, "%s/%s.afa" % (dc.get_mafft_path(prot_id), species))
            print cmd
            os.system(cmd)
            
            os.remove(msa_species_path)


if __name__ == '__main__':
    #create_msa_alignments()
    create_species_msa_alignments()
