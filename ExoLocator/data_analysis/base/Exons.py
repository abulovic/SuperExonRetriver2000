'''
Created on Apr 30, 2012

@author: intern
'''
from data_analysis.alignment.parsing.BlastParser    import parse_blast_output
from data_analysis.containers.EnsemblExonContainer  import EnsemblExonContainer
from data_analysis.alignment.parsing.SWParser       import parse_SW_output
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio import SeqIO
from data_analysis.containers.DataMapContainer import DataMapContainer
from Bio.Alphabet import IUPAC

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
        self.ref_protein_id = data_map_key[0]
        self.species        = data_map_key[1]
        self.ref_species    = ref_species
        self.exon_type      = exon_type
        self.id             = self._generate_id ()
        
    def load_exons (self):
        if self.exon_type == "blastn" or self.exon_type == "tblastn":
            exon_dict = parse_blast_output(self.ref_protein_id, self.species, self.exon_type)
        elif self.exon_type.lower() == "sw_gene" or self.exon_type.lower() == "sw_exon":
            exon_dict = parse_SW_output (self.ref_protein_id, self.species, self.exon_type)
        else:
            raise KeyError("No parsing option for exon_type %s" % self.exon_type)
        if not exon_dict:
            self.alignment_exons = None
        else:
            self.alignment_exons = exon_dict
        return exon_dict
    
    def set_exon_ordinals (self):
        '''
        Assign ordinals to alignment exons.
        This serves as an aid to correct ordering of exons and the 
        discarding of false positives later employed by the
        AlignmentStatistics
        '''
        ensembl_container = EnsemblExonContainer.Instance()
        for ref_exon_id, al_exons in self.alignment_exons.items():
            ref_exon = ensembl_container.get(ref_exon_id)
            # see how many alignment exons are there, and then give them foating ordinals
            num_of_exons = len(al_exons)
            if num_of_exons == 1:
                al_exons[0].set_ordinal(ref_exon.ordinal)
            else:
                tmp_ordinal = ref_exon.ordinal
                for al_exon in sorted(al_exons, key = lambda al_exon : al_exon.alignment_info["query_start"]):
                    al_exon.set_ordinal(tmp_ordinal)
                    tmp_ordinal = float(tmp_ordinal) + 1./num_of_exons
            
        pass
    
    def get_ordered_exons (self):
        
        self.set_exon_ordinals()            
        ordered_exons = sorted(self.alignment_exons.values(), 
                              key = lambda al_exons: al_exons[0].ordinal)
        
        return ordered_exons
            
        
        
        
    
    def export_to_fasta (self, file_name):
        
        seq_records = []
        fasta_file = open(file_name, 'w')
        
        for ref_exon_id, al_exons in self.alignment_exons.items():
            i = 1
            for al_exon in al_exons:
                record = SeqRecord(Seq(al_exon.alignment_info["sequence"], IUPAC.ambiguous_dna), 
                                   id = "%s_%d"%(ref_exon_id, i), 
                                   description="%d|%d"%(al_exon.alignment_info["query_start"], al_exon.alignment_info["query_end"]))
                i += 1
                seq_records.append(record)
                
        SeqIO.write(seq_records, fasta_file, "fasta")
            
    
    def get_flattened_exons (self):
        '''
        Creates a list with elements of type Exon, not of list containing possibly
        multiple Exons as is the case with alignment_exons attribute
        '''      
        exon_list = []
        for al_exons in self.alignment_exons:
            exon_list.extend(al_exons)
            
        return exon_list
        
    def _generate_id(self):
        (spec1, spec2) = self.species.split("_")
        spec_name = "%s%s" % (spec1[0:2], spec2[0:2])
        return "%s_%s_%s" % (self.species, self.ref_protein_id, self.exon_type)
     
    
        
        