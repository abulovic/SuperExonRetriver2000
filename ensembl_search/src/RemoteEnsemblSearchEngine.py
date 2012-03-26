'''
Created on Mar 26, 2012

@author: Ana
'''

from cogent.db.ensembl import Species, Genome
from cogent.db.ensembl.region import *
from cogent.db.ensembl.region import Transcript

class RemoteEnsemblSearchEngine(object):
    '''
    Class provides set of methods to retrieve core Ensembl data for particular species. 
    The access is remote, not local. It is to be expected that this
    search will be much slower than the local search.
    '''

    def __init__(self, species, release=66):
        '''
        Creates a genome for the species
        @param release: release version (66 if not specified)
        @param species: Ensembl species name
        '''
        self.release = release
        self.species = species
        
        self.genome = Genome(Species = self.species, Release = self.relase, account = None)
        
    def getGene (self, geneStableID, fileName):
        '''
        Retrieves a gene for the species
        @param geneStableID: Ensembl gene ID (ENSG000001234)
        @param fileName: if specified, a name of the file to write the sequence to (fasta format)
        @return: A Gene object
        '''
        gene = self.species.getGeneByStableId (StableId = geneStableID)
        
        if (fileName != None):
            fasta_file = open(fileName, 'w');
            fasta_file.write(">%s %s\n%s", gene.StableId, gene.Status, gene.Seq);
            
        return gene
    
    def getGeneSequence (self, geneStableId, fileName):
        '''
        Retrieves the gene sequence for ensembl stable id
        @param geneStableId: Ensembl gene ID (ENSG000001234)
        @param fileName: if specified, a name of the file to write the sequence to (fasta format)
        @return: string containing gene sequence
        '''
        gene = self.getGene (geneStableId, None)
        return gene.Seq
    
    def getTranscript (self, geneStableId, gene, transcriptStableId, fileName):
        '''
        Retrieves the transcript for the gene
        @param geneStableId: if specified, first the gene is retrieved, then transcript
        @param gene: A Gene object, if specified, no need to specify geneStableId
        @param transcriptStableId: Ensembl transcript stable id
        @return: A Transcript object 
        '''
        
        if (gene == None and geneStableId == None):
            print ("ERROR: Neither gene nor geneStableID specified.\n")
            return None
        if (gene == None):
            gene = self.genome.getGeneByStableId(StableId=geneStableId)
        
        transcript = gene.getMember (StableId = transcriptStableId)
        
        if (fileName != None):
            fastaFile = open(fileName, 'w');
            fastaFile.write(">%s dna:%s %s\n", transcript.StableId, gene.StableId, transcript.Status, transcript.Seq)
            
        return transcript
    
    def getTranscriptSequence (self, geneStableId, gene, transcriptStableId, fileName):
        '''
        Retrieves the transcript sequence for ensembl stable id
        @param geneStableId: if specified, first the gene is retrieved, then transcript
        @param gene: A Gene object, if specified, no need to specify geneStableId
        @param transcriptStableId: Ensembl transcript stable id
        @param fileName: if specified, a name of the file to write the sequence to (fasta format)
        @return: string containing gene sequence
        '''
        transcript = self.getTranscript(geneStableId, gene, transcriptStableId, fileName)
        return transcript.Seq
    
    def getExonsForTranscript (self, transcript, gene, fileName, exon_type = "Untranslated"):
        '''
        Retrieves exon for a particular transcript 
        @param gene: Gene / String object 
        @param transcript: Transcript / String object
        @param exon_type: Translated / Untranslated
        @param fileName: if specified, name of fasta file to which the exons will be written
        @return: list of exons
        '''
        exons = []
        if (not isinstance(transcript, Transcript)):
            if (not isinstance (gene, Gene)):
                gene = self.genome.getGeneByStableId(StableId=gene)
            transcript = gene.getMember(StableId=transcript)
        
        # get exons dependent on their type
        if (exon_type == "Untranslated"):
            exons = list(transcript.Exons)
        else:
            exons = list(transcript.TranslatedExons)
            
        if (fileName != None):
            fastaFile = open(fileName, 'w')
            for exon in exons:
                fastaFile.write(">%s %s dna:%s %s\n" % (transcript.StableId, exon.StableId, gene.StableId, transcript.Status))
                fastaFile.write("%s" % exon.Seq)
                
        return exons
        
        
    
        