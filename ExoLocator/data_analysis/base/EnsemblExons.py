'''
Created on Apr 30, 2012

@author: intern
'''

from pipeline.utilities.DirectoryCrawler import DirectoryCrawler

from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from data_analysis.base.EnsemblExon import EnsemblExon
from Bio.SeqRecord import SeqRecord
from utilities.Logger import Logger
from data_analysis.containers.DataMapContainer import DataMapContainer
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from utilities.FileUtilities import get_reference_species_dictionary
from data_analysis.utilities.ExonUtils import remove_UTR_ensembl_exons

class EnsemblExons(object):
    '''
    classdocs
    '''


    def __init__(self, data_map_key, ref_species=None):
        '''
        Constructor
        '''
        self.ref_protein_id    = data_map_key[0]
        self.species        = data_map_key[1]
        if not ref_species:
            spec_dict = get_reference_species_dictionary()
            ref_species = spec_dict[data_map_key[1]]
        self.ref_species    = ref_species
        self.exons          = {}
        
    def get_exon_file_path (self):
        
        dc = DirectoryCrawler()
        return "{0}/{1}.fa".format(dc.get_exon_ensembl_path(self.ref_protein_id), self.species)
    
    def load_exons (self):
        
        data_map_container = DataMapContainer.Instance()
        logger = Logger.Instance()
        containers_logger = logger.get_logger('containers')
        
        data_map = data_map_container.get((self.ref_protein_id, self.species))
        self.strand = data_map.strand
        
        fasta_path = self.get_exon_file_path()
        try:
            fasta = open(fasta_path, 'r')
        except IOError:
            containers_logger.error("%s,%s,%s" % (self.ref_protein_id, self.species, "Loading ensembl exons failed."))
            return None
         
        exon_list = []
        for seq_record in SeqIO.parse(fasta, "fasta", unambiguous_dna):
            (start, stop, transcript_id, exon_id, strand) = seq_record.id.split('|')
            if (int(strand) == 1):
                self.strand = 1
                exon = EnsemblExon((self.ref_protein_id, self.species), exon_id, start, stop, strand, seq_record.seq)
            else:
                self.strand = -1
                exon = EnsemblExon((self.ref_protein_id, self.species), exon_id, stop, start, strand, seq_record.seq)
            exon_list.append(exon)
        fasta.close()
        self.exons = dict([(exon.exon_id, exon) for exon in exon_list])
        
        # assign orinals to exons
        ordinal = 1
        if self.strand == 1:
            for exon in sorted (self.exons.values(), key = lambda exon: exon.start ):
                exon.set_exon_ordinal(ordinal)
                ordinal += 1
        else:
            for exon in sorted (self.exons.values(), key = lambda exon: exon.start, reverse = True):
                exon.set_exon_ordinal(ordinal)
                ordinal += 1
        
        return exon_list
        
    def get_ordered_exons (self):
        if self.strand == 1:
            return sorted (self.exons.values(), key = lambda exon: exon.start )
        else:
            return sorted (self.exons.values(), key = lambda exon: exon.start, reverse = True)
        
    def get_cDNA (self):
        exons = self.get_ordered_exons()
        merged_exons_seq = Seq("", IUPAC.ambiguous_dna)
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
    
    def get_coding_cDNA(self):
        '''
        Gets the cDNA without the UTR regions
        '''
        dmc = DataMapContainer.Instance()
        dm_key = dmc.get((self.ref_protein_id, self.species))
        new_exons = remove_UTR_ensembl_exons(self.ref_protein_id, self.species, self.get_ordered_exons())
        
        coding_cDNA = Seq("", IUPAC.ambiguous_dna)
        for exon in new_exons:
            coding_cDNA += exon.sequence
            
        return coding_cDNA
    
    def export_coding_exons_to_fasta (self, fasta_file):
        
        dmc = DataMapContainer.Instance()
        data_map = dmc.get((self.ref_protein_id, self.species))
        
        new_exons = remove_UTR_ensembl_exons(self.ref_protein_id, self.species, self.get_ordered_exons())
        exon_records = []
        " >969067|969174|ENSAMET00000013141|ENSAMEE00000125733|1"
        for exon in new_exons:
            exon_id = "%d|%d|%s|%s|%d" % (exon.start, exon.stop, data_map.transcript_id, exon.exon_id, exon.strand)
            record = SeqRecord(seq = exon.sequence, id = exon_id, description = "")
            exon_records.append(record)
            
        SeqIO.write(exon_records, fasta_file, "fasta")
        
        
            
        
        
        
    
if __name__ == '__main__':
    ee = EnsemblExons(("ENSP00000341765", "Ailuropoda_melanoleuca"), "Homo_sapiens")
    print ee.get_exon_file_path()
    ee.load_exons()
    for exon in ee.get_ordered_exons():
        print exon.start, exon.stop, exon.ordinal
    (cdna_seq, exon_locations) = ee.get_cDNA()
    print exon_locations
    print cdna_seq.seq
    
        
        