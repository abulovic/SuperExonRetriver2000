'''
Created on Apr 30, 2012

@author: intern
'''

from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from data_analysis.base.EnsemblExon import EnsemblExon
from Bio.SeqRecord import SeqRecord

class EnsemblExons(object):
    '''
    classdocs
    '''


    def __init__(self, data_map_key, ref_species):
        '''
        Constructor
        '''
        self.ref_protein    = data_map_key[0]
        self.species        = data_map_key[1]
        self.ref_species    = ref_species
        self.exons          = []
        
    def get_exon_file_path (self):
        
        dc = DirectoryCrawler()
        return "{0}/{1}.fa".format(dc.get_exon_ensembl_path(self.ref_protein), self.species)
    
    def load_exons (self):
        fasta_path = self.get_exon_file_path()
        fasta = open(fasta_path, 'r')
        exon_list = []
        for seq_record in SeqIO.parse(fasta, "fasta", unambiguous_dna):
            (start, stop, transcript_id, exon_id, strand) = seq_record.id.split('|')
            if (int(strand) == 1):
                exon = EnsemblExon((self.ref_protein, self.species), exon_id, start, stop, strand, seq_record.seq)
            else:
                exon = EnsemblExon((self.ref_protein, self.species), exon_id, stop, start, strand, seq_record.seq)
            exon_list.append(exon)
        fasta.close()
        self.exons = exon_list
        return exon_list
    
    def get_ordered_exons (self):
        if self.exons[0].strand == 1:
            return sorted (self.exons, key = lambda exon: exon.start )
        else:
            return sorted (self.exons, key = lambda exon: exon.start, reverse = True)
        
    def get_cDNA (self):
        exons = self.get_ordered_exons()
        merged_exons_seq = ""
        exon_locations = {}
        start = 1
        end = 1
        exon_id = 1
        for exon in exons:
            end += len(exon.sequence)-1
            exon_locations[exon_id] = (start, end)
            exon_id += 1
            start = end+1
            end = start
            merged_exons_seq += exon.sequence
        cdna_seq = SeqRecord(seq=merged_exons_seq, id=self.species, description="")
        return (cdna_seq, exon_locations)
        
        
        
    
if __name__ == '__main__':
    ee = EnsemblExons(("ENSP00000341765", "Ailuropoda_melanoleuca"), "Homo_sapiens")
    print ee.get_exon_file_path()
    ee.load_exons()
    for exon in ee.get_ordered_exons():
        print exon.start, exon.stop
    (cdna_seq, exon_locations) = ee.get_cDNA()
    print exon_locations
    print cdna_seq.seq
    
        
        