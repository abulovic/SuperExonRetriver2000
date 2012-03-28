'''
Created on Mar 26, 2012

@author: anana
'''

import ConfigParser;
import os, re;

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
        <Attribute name = "exon_chrom_start" />
        <Attribute name = "exon_chrom_end" />
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "gene_exon" />
        <Attribute name = "ensembl_exon_id" />
    </Dataset>
</Query>
        '''
        
    def parseDescriptionFile (self, descrFileName):
        descrFile = open(descrFileName, 'r')
        
        transcripts = []
        abinitioSpecies = []
        knownSpecies = []
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
                    knownSpecies.append(species)
                else :
                    abinitioSpecies.append(species)
            i = i + 1
            
        
        descrFile.close()
        return (transcripts,knownSpecies,abinitioSpecies)
        
        
    def populateExonDatabase (self, proteinDirectory):
        '''
        Retrieves exons for the species with known genes (not abinitio)
        @param proteinDirectory: basically just the name of the protein (not the whole path)
        '''
        print "Exons for protein: %s" % proteinDirectory
        descrFileName = "%s/%s/%s" % (self.sessionDir, proteinDirectory, self.descrFileName)
        # find transcripts for all the species with known genes
        (transcripts, knownSpecies, abinitioSpecies) = self.parseDescriptionFile(descrFileName)
        
        print transcripts
        
        exonDatabase   = "%s/%s/exons/db/" % (self.sessionDir, proteinDirectory)
        queryFileName = "%s/%s/%s" % (self.sessionDir, proteinDirectory, self.tmpXmlFile)
        
        for (species, ensemblSpeciesName, transcriptId) in transcripts :
            queryFile = open(queryFileName, 'w')
            query = self.templateXML % (ensemblSpeciesName, transcriptId)
            queryFile.write(query)
            queryFile.close()
            exonFileName = "%s/%s.fa" % (exonDatabase, species)
            cmd = "perl %s %s > %s" % (self.perlBiomartScript, queryFileName, exonFileName)
            print cmd
            os.system(cmd)
            
            
        for (species, ensemblSpeciesName, transcriptId) in transcripts :
            exonFileName = "%s/%s.fa" % (exonDatabase, species)
            self.reorderExons(exonFileName)
            
        return (knownSpecies, abinitioSpecies)
            
            
    def reorderExons (self, exonFileName):
        
        print "reordering exons... %s" % exonFileName
        headerPattern = re.compile(r'>(\d+)\|(\d+)\|(\w+)\|(\w+)')
        exonFile = open(exonFileName, 'r')
        exonSequence = ""
        exons = []
        headerInfoList = []
        
        for line in exonFile.readlines():
            print line
            headerInfo = re.match(headerPattern, line)
            if (headerInfo != None):
                headerInfo = headerInfo.groups()
                if exonSequence != "":
                    print headerInfoList
                    exon = AuxExon(headerInfoList[0], headerInfoList[1], headerInfoList[2], headerInfoList[3], exonSequence)
                    print "Adding exon"
                    exons.append(exon)
                    exonSequence = ""
                headerInfoList = headerInfo
            else:
                exonSequence = exonSequence + line
            
        exonFile.close()
        # if no appropriate exons found (genewise output, then do nothing
        if (len(exons) == 0):
            return
        exonFile = open(exonFileName, 'w')
        exonCount = 1    
        for exon in sorted(exons):
            exonFile.write(">%d %s|%s|%s|%s\n%s" % (exonCount, exon.beg, exon.end, exon.transId, exon.exonId, exon.seq))
            ">%d %s|%s|%s|%s\n%s" % (exonCount, exon.beg, exon.end, exon.transId, exon.exonId, exon.seq)
            exonCount = exonCount + 1
        exonFile.close()
       
       
       
class AuxExon (object):
    def __init__ (self, beg, end, transId, exonId, seq):
        self.beg = beg
        self.end = end
        self.transId = transId
        self.exonId = exonId
        self.seq = seq
    def __lt__ (self, other):
        if (other.beg > self.beg):
            return True
        else:
            return False
            
if __name__ == '__main__':
    biomart = BioMartSearchEngine()
    biomart.populateExonDatabase("ENSP00000252816")
        
        
        
        