'''
Created on Mar 26, 2012

@author: anana
'''

import ConfigParser;
import os;

class BioMartSearchEngine(object):
    '''
    Queries the BioMart for exons
    '''
    
    def __init__ (self):
        '''
        Loads the project specific configuration
        '''
        config = ConfigParser.RawConfigParser()
        config.read("../../config.cfg")
        
        rootDir             = config.get("Project root", "project_root_folder")
        self.sessionDir     = config.get("Project root", "session_folder")
        self.sessionDir     = "%s/%s" % (rootDir, self.sessionDir)
        self.descrFileName  = config.get("Session files", "descr_output")
        
        self.perlBiomartScript  = "%s/ensembl_search/src/BioMartRemoteAccess.pl" % (rootDir)
        self.tmpXmlFile         = "Query.xml"
        
        self.templateXML = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
    <Dataset name = "%s_gene_ensembl" interface = "default" >
        <Filter name = "ensembl_transcript_id" value = "%s"/>
        <Attribute name = "gene_exon" />
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "ensembl_transcript_id" />
    </Dataset>
</Query>
        '''
        
    def parseDescriptionFile (self, descrFileName):
        descrFile = open(descrFileName, 'r')
        
        transcripts = []
        abinitioSpecies = []
        i = 0
        for line in descrFile.readlines():
            line = line.strip()
            
            if (len(line) == 0):
                continue
            if (i % 2 == 0):
                species= line
                ensemblSpeciesName = "%s%s" % (species[0].lower(), species.split('_')[1])  
            else :
                data = line.split()
                if (len(data) == 5):
                    transcriptId = data[4]
                    transcripts.append((species, ensemblSpeciesName, transcriptId))
                else :
                    abinitioSpecies.append(species)
            i = i + 1
            
        
        descrFile.close()
        return (transcripts,abinitioSpecies)
        
        
    def populateExonDatabase (self, proteinDirectory):
        '''
        Retrieves exons for the species with known genes (not abinitio)
        @param proteinDirectory: basically just the name of the protein (not the whole path)
        '''
        print "Exons for protein: %s" % proteinDirectory
        descrFileName = "%s/%s/%s" % (self.sessionDir, proteinDirectory, self.descrFileName)
        # find transcripts for all the species with known genes
        (transcripts, abinitioSpecies) = self.parseDescriptionFile(descrFileName)
        
        exonDatabase   = "%s/%s/exons/db/" % (self.sessionDir, proteinDirectory)
        
        for (species, ensemblSpeciesName, transcriptId) in transcripts :
            queryFile = open(self.tmpXmlFile, 'w')
            query = self.templateXML % (ensemblSpeciesName, transcriptId)
            queryFile.write(query)
            queryFile.close()
            exonFileName = "%s/%s.fa" % (exonDatabase, species)
            cmd = "perl %s %s > %s" % (self.perlBiomartScript, self.tmpXmlFile, exonFileName)
            os.system(cmd)
            
        return abinitioSpecies
            
            
            
if __name__ == '__main__':
    biomart = BioMartSearchEngine()
    biomart.populateExonDatabase("parf")
        
        
        
        
         