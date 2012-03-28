'''
Created on Mar 26, 2012

@author: Ana
'''

from cogent.db.ensembl import  Genome
from cogent.db.ensembl.region import Gene, Transcript

class RemoteEnsemblSearchEngine(object):
    '''
    Class provides set of methods to retrieve core Ensembl data for particular species. 
    The access is remote, not local. It is to be expected that this
    search will be much slower than the local search.
    NOTICE: If you get the gene before you want to get the transcript (
    or the transcript before you need the exons), be sure to pass
    the object instance (Gene / Transcript) instead of the gene ID, because it is much 
    faster (one inquiry less on the remote database).
    '''

    def __init__(self, species, release=66):
        '''
        Creates a genome for the species
        @param release: release version (66 if not specified)
        @param species: Ensembl species name
        '''
        self.release = release
        self.species = species
        
        self.genome = Genome(Species = self.species, Release = self.release, account = None)
        
    def setSessionFolder (self, sessionFolder):
        '''
        Sets the session folder
        '''
        
        self.sessionFolder = sessionFolder
        
    def getGene (self, geneStableID, fileName=None):
        '''
        Retrieves a gene for the species
        @param geneStableID: Ensembl gene ID (ENSG000001234)
        @param fileName: if specified, a name of the file to write the sequence to (fasta format)
        @return: A Gene object
        '''
        gene = self.genome.getGeneByStableId (StableId = geneStableID)
        
        
        if (fileName != None):
            fasta_file = open(fileName, 'w');
            fasta_file.write(">%s %s\n%s", gene.StableId, gene.Status, gene.Seq);
            fasta_file.close()
            
        return gene
    
    def getGeneSequence (self, gene, fileName=None):
        '''
        Retrieves the gene sequence for ensembl stable id
        @param gene: Gene object / Ensembl gene ID (ENSG000001234)
        @param fileName: if specified, a name of the file to write the sequence to (fasta format)
        @return: string containing gene sequence
        '''
        if (not isinstance(gene, Gene)):
            gene = self.getGene (gene, fileName)
        return gene.Seq
    
    def getTranscript (self, gene, transcriptStableId, fileName=None):
        '''
        Retrieves the transcript for the gene
        @param gene: Gene object / Ensembl gene ID (ENSG000001234)
        @param transcriptStableId: Ensembl transcript stable id
        @return: A Transcript object 
        '''
        
        # if the gene object is not accessible, fetch it
        if (not isinstance(gene, Gene)):
            gene = self.genome.getGeneByStableId(gene)
        
        transcript = gene.getMember (StableId = transcriptStableId)
        
        if (fileName != None):
            fastaFile = open(fileName, 'w');
            fastaFile.write(">%s dna:%s %s\n", transcript.StableId, gene.StableId, transcript.Status, transcript.Seq)
            fastaFile.close()
            
        return transcript
    
    def getTranscriptSequence (self, geneStableId, gene, transcriptStableId, fileName=None):
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
    
    def getExonsForTranscript (self, transcript, gene, fileName = None, exon_type = "Standard"):
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
        if (exon_type == "Standard"):
            exons = list(transcript.Exons)
        elif (exon_type == "Translated"):
            exons = list(transcript.TranslatedExons)
        else:
            print("ERROR: Exon type needs to be Standard / Translated. Using Standard by default.")
            exons = list(transcript.Exons)
        
        
        #Write exons to a file. 
        #Fasta format (Ensembl like):
        #>transcript_id exon_id dna:gene_id transcript_status    
        if (fileName != None):
            fastaFile = open(fileName, 'w')
            for exon in exons:
                print exon.StableId
                fastaFile.write(">%s %s dna:%s %s\n" % (transcript.StableId, exon.StableId, gene.StableId, transcript.Status))
                fastaFile.write("%s\n\n" % exon.Seq)
            fastaFile.close()
                
        return exons
    
    def fetchExonsFromDescriptionFile (self, descrFileName):
        descrFile = open(descrFileName, 'r')
        i           = 0
        species     = ""
        proteinId   = ""
        location    = ""
        status      = ""
        geneId      = ""
        exon_database = "%s/exons/db/" % self.sessionFolder
        
        for line in descrFile.readlines():
            line = line.strip()
            #empty line
            if (len(line) == 0):
                continue
            if (i % 2 == 0):
                species = line
            elif (i % 2) == 1:
                data = line.split();
                if (len(data) == 3):
                    (proteinId, location, status) = tuple(data)
                else:
                    (proteinId, location, status, geneId, transcriptId) = tuple(data)
            
            # there is no gene, transcript and exon information 
            # in Ensembl about proteins predicted by Genscan
            if (status == "abinitio" and i % 2 == 1):
                i = i+1
                continue
            
            if (i % 2 == 1):
                engine = RemoteEnsemblSearchEngine(species, 66)
                try:
                    engine.getExonsForTranscript(transcriptId, geneId, "%s/%s.fa" % (exon_database, species), "Standard")
                except RuntimeError:
                    print "Species %s not available by PyCogent." % species
                
            i = i+1
                    
            
            
                    
                
                
                
                


if __name__ == '__main__':
    remote_ens = RemoteEnsemblSearchEngine('human', 66)
    #gene = remote_ens.getGene("ENSG00000196642", None)
    #exons = remote_ens.getExonsForTranscript("ENST00000311502", "ENSG00000196642", None, "Standard")
    #for exon in exons:
    #    print exon.Seq
    remote_ens.setSessionFolder("/home/intern/Project/workspaceBIO/sessions/parf")
    remote_ens.fetchExonsFromDescriptionFile("/home/intern/Project/workspaceBIO/sessions/parf/mutual_best_out.descr")
        
        