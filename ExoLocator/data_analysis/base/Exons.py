'''
Created on Apr 30, 2012

@author: intern
'''
from data_analysis.alignment.parsing.BlastParser import parse_blast_output
from data_analysis.containers.EnsemblExonContainer import EnsemblExonContainer

class Exons(object):
    '''
    Information about particular exons.
    Exons can be of 6 types: ensembl, genewise, blastn, tblastn, SW_gene and SW_exon.
    Instance of Exons class knows where its data file is located and how to read it.
    Contains exons in form of a dictionary. Key is reference exon id, and values are
    Exon instances which contain information on available alignments to that reference exon.
    '''

    def __init__(self, data_map_key, ref_species, exon_type):
        '''
        Constructor
        '''
        self.ref_protein    = data_map_key[0]
        self.species        = data_map_key[1]
        self.ref_species    = ref_species
        self.exon_type      = exon_type
        self.id             = self._generate_id ()
        
    def load_exons (self):
        exon_dict = parse_blast_output(self.ref_protein, self.species, self.exon_type)
        self.alignment_exons = exon_dict
        return exon_dict
        
    def set_reference_exons (self):
        ens_exon_container = EnsemblExonContainer.Instance()
        for ref_exon_id in self.alignment_exons:
            if ens_exon_container.get(ref_exon_id):
                self.alignment_exons[ref_exon_id].set_reference_exon(ens_exon_container.get(ref_exon_id))
            
        
    def _generate_id(self):
        (spec1, spec2) = self.species.split("_")
        spec_name = "%s%s" % (spec1[0:2], spec2[0:2])
        return "%s_%s_%s" % (self.species, self.ref_protein, self.exon_type)
     
    
        
        