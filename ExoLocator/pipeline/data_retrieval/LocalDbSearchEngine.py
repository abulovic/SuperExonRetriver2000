'''
Created on Apr 15, 2012

@author: intern
'''

from utils.ConfigurationReader import * 
from utils.FileUtils import *
from pipeline.utilities.AlignmentCommandGenerator import *
import os


acg = AlignmentCommandGenerator()
cr = ConfigurationReader.Instance()

def populate_proteins_for_prot_id (protein_id):
    protein_dir = cr.get_value('root', 'session_dir') + "/" + protein_id + "/" + cr.get_value('sequence', 'protein') + "/"
    print protein_dir
    
    (genes_known, genes_abinitio) = get_gene_regions(protein_id)
    
    for key, value in genes_known.items():
        (location_type, assembly, location_id, seq_begin, seq_end, strand) = value
        output_file = cr.get_value('root', 'session_dir') + "/" + protein_id + "/" + cr.get_value('sequence', 'root') + "/" + cr.get_value('sequence', 'gene') + "/" + key + ".fasta"
        fastacmd = acg.generate_fastacmd_command(location_id, key, "dna", location_type, output_file, 0, strand, int(seq_begin), int(seq_end))
        print fastacmd
        os.system(fastacmd)
        
    for key, value in genes_abinitio.items():
        (location_type, assembly, location_id, seq_begin, seq_end, strand) = value
        output_file = cr.get_value('root', 'session_dir') + "/" + protein_id + "/" + cr.get_value('sequence', 'root') + "/" + cr.get_value('sequence', 'gene') + "/" + key + ".fasta"
        fastacmd = acg.generate_fastacmd_command(location_id, key, "dna", location_type, output_file, 0, strand, int(seq_begin), int(seq_end))
        print fastacmd
        os.system(fastacmd)
        
    
    
    
if __name__ == '__main__':
    populate_proteins_for_prot_id("ENSP00000311134")
        
        