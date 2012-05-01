'''
Created on Apr 30, 2012

@author: intern
'''

class Exons(object):
    '''
    Information about particular exons.
    Exons can be of 6 types: ensembl, genewise, blastn, tblastn, SW_gene and SW_exon.
    Instance of Exons class knows where its data file is located and how to read it.
    '''


    def _generate_id(self):
        (ref_prot_id, species) = self.data_map_key
        (spec1, spec2) = species.split()
        spec_name = "%s%s" % (spec1[0:2], spec2[0:2])
        return "%s_%s_%s" % (spec_name, ref_prot_id, self.exon_type)
        
    
    
    def __init__(self, data_map_key, ref_species, exon_type):
        '''
        Constructor
        '''
        self.data_map_key   = data_map_key
        self.ref_species    = ref_species
        self.exon_type      = exon_type
        self.id             = self._generate_id ()
        
    
        
        