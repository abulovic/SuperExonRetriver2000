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
        self.config = ConfigParser.RawConfigParser()
        self.config.read("../../config.cfg")
        
        rootDir             = self.config.get("Project root", "project_root_folder")
        self.sessionDir     = self.config.get("Project root", "session_folder")
        self.sessionDir     = "%s/%s" % (rootDir, self.sessionDir)
        self.descrFileName  = self.config.get("Session files", "descr_output")
        
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
        <Attribute name = "strand" />
    </Dataset>
</Query>
        '''
        
    def generate_file_name (self, masked, species, id_type, id):
        ensembldb = self.config.get("Ensembl cfg", "ensembldb")
        file_name = "%s/%s/dna/" % (ensembldb, species.lower())
        tmp_file=""
        for file in os.listdir(file_name):
            if (file != "README"):
                tmp_file = file 
                break
        m = re.findall ('(.*).dna', tmp_file)   
        if (masked == True):
            file_name = "%s/%s.dna_rm." % (file_name, m[0])
        else :
            file_name = "%s/%s.dna." % (file_name, m[0])
        if (id_type == 'chromosome'):
            file_name = "%schromosome.%s.fa" % (file_name, id)
        else :
            file_name = "%stoplevel.fa" % (file_name)
        return file_name
        
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
                    locationData = data[1].split(':')
                    transcriptLocationType = locationData[0]
                    transcriptLocationID = locationData[2]
                    strand = locationData[5]
                    transcripts.append((species, ensemblSpeciesName, transcriptId, transcriptLocationType, transcriptLocationID, strand))
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
        
        for (species, ensemblSpeciesName, transcriptId, transcriptLocationType, transcriptLocationID, strand) in transcripts :
            queryFile = open(queryFileName, 'w')
            query = self.templateXML % (ensemblSpeciesName, transcriptId)
            queryFile.write(query)
            queryFile.close()
            exonFileName = "%s/%s.fa" % (exonDatabase, species)
            cmd = "perl %s %s > %s" % (self.perlBiomartScript, queryFileName, exonFileName)
            print cmd
            while (True):
     
                os.system(cmd)
                fileSize = os.path.getsize(exonFileName)
                if (fileSize != 0):
                    break
         
            
        for (species, ensemblSpeciesName, transcriptId, transcriptLocationType, transcriptLocationID, strand) in transcripts :
            exonFileName = "%s/%s.fa" % (exonDatabase, species)
            exons = self.reorderExons(exonFileName)
            self.removeFrameshiftIntrons(species, exonFileName, exons,transcriptId, transcriptLocationType, transcriptLocationID, strand)
            
        return (knownSpecies, abinitioSpecies)
            
            
    def reorderExons (self, exonFileName):
        
        print "reordering exons... %s" % exonFileName
        headerPattern = re.compile(r'>(\d+)\|(\d+)\|(\w+)\|(\w+)\|(.*)')
        
        fileSize = os.path.getsize(exonFileName)
        if (fileSize == 0):
            os.remove(exonFileName)
            return
        
        exonFile = open(exonFileName, 'r')
        exonSequence = ""
        exons = []
        headerInfoList = []
        
        for line in exonFile.readlines():
            headerInfo = re.match(headerPattern, line)
            if (headerInfo != None):
                headerInfo = headerInfo.groups()
                if exonSequence != "":
                    exon = _AuxExon(headerInfoList[0], headerInfoList[1], headerInfoList[2], headerInfoList[3], headerInfoList[4], exonSequence)
                    exons.append(exon)
                    exonSequence = ""
                headerInfoList = headerInfo
            else:
                exonSequence = exonSequence + line
        exon = _AuxExon(headerInfoList[0], headerInfoList[1], headerInfoList[2], headerInfoList[3], headerInfoList[4], exonSequence)
        exons.append(exon) 
            
        exonFile.close()
        # if no appropriate exons found (genewise output, then do nothing
        if (len(exons) == 0):
            return
        return sorted(exons)
        '''
        exonFile = open(exonFileName, 'w')
        exonCount = 1    
        for exon in sorted(exons):
            exonFile.write(">%d %s|%s|%s|%s\n%s" % (exonCount, exon.beg, exon.end, exon.transId, exon.exonId, exon.seq))
            ">%d %s|%s|%s|%s\n%s" % (exonCount, exon.beg, exon.end, exon.transId, exon.exonId, exon.seq)
            exonCount = exonCount + 1
        exonFile.close()
        '''
        
    def removeFrameshiftIntrons (self, species, exonOutputFileName, exons, transcriptId, transcriptLocationType, transcriptLocationID, strand):
        
        print "Removing frameshift introns for species %s, transcript %s" % (species, transcriptId)
        exonFile = open(exonOutputFileName, 'w')
        dnaDB = self.generate_file_name (False, species, transcriptLocationType, transcriptLocationID)
        if (transcriptLocationType == "chromosome"):
            fastacmdSearchID = "chrom%s" % str(transcriptLocationID)
        else:
            fastacmdSearchID = transcriptLocationID
        
        cmd = "fastacmd -d %s -s %s -S %s -L %s,%s -p F -o %s"  %           (dnaDB,          # database name
                                                                             fastacmdSearchID,                # id
                                                                             strand,     # strand
                                                                             "%d",     # seq beginning
                                                                             "%d",     # seq ending
                                                                             "tmp.fa")
        
        if (exons == None):
            return
        
        newExons = []
        
        startingLocation = exons[0].beg
        endLocation = exons[0].end
        
        for exonCounter in range (1, len(exons)):
            exon = exons[exonCounter]
            if (species == "Sorex_araneus"):
                print exon.seq
            
            if (strand == "1"):

                if (abs(exon.beg - endLocation) <= 10):
                    print "MERGING!!"
                    endLocation = exon.end
                else :
                    cmdEx = cmd % (startingLocation, endLocation)
                    print cmdEx
                    os.system(cmdEx)
                    exonSeq = open("tmp.fa", 'r')
                    exonSeq.readline()
                    newExons.append(_AuxExon(startingLocation, endLocation, transcriptId, "", strand, exonSeq.read()))
                    exonSeq.close()
                    
                    startingLocation = exon.beg
                    endLocation = exon.end
                    
            else:
                if (abs(exon.end - startingLocation) <= 10):
                    print "MERGING!!"
                    startingLocation = exon.beg
                else:
                    cmdEx = cmd % (startingLocation, endLocation)
                    print cmdEx
                    os.system(cmdEx)
                    exonSeq = open("tmp.fa", 'r')
                    exonSeq.readline()
                    newExons.append(_AuxExon(startingLocation, endLocation, transcriptId, "", strand, exonSeq.read()))
                    exonSeq.close()
                    
                    startingLocation = exon.beg
                    endLocation = exon.end
                
        cmdEx = cmd % (startingLocation, endLocation)
        print cmdEx
        os.system(cmdEx)
        exonSeq = open("tmp.fa", 'r')
        exonSeq.readline()
        newExons.append(_AuxExon(startingLocation, endLocation, transcriptId, "", strand, exonSeq.read()))
        exonSeq.close()
        exonFile.close()
        print "Species %s, number of exons:: %d" % (species, len(newExons))
        
        exonFile = open(exonOutputFileName, 'w')
        exonCount = 1    
        for exon in sorted(newExons):
            exonFile.write(">%d %s|%s|%s\n%s" % (exonCount, exon.beg, exon.end, exon.transId, exon.seq))
            exonCount = exonCount + 1
        exonFile.close()
       
       
       
class _AuxExon (object):
    def __init__ (self, beg, end, transId, exonId, strand, seq):
        self.beg = int(beg)
        self.end = int(end)
        self.transId = transId
        self.exonId = exonId
        self.strand = strand
        self.seq = seq
    def __lt__ (self, other):
        if (self.strand == "1"):
            if (other.beg > self.beg):
                return True
            else:
                return False
        else:
            if (other.beg > self.beg):
                return False
            else:
                return True
            
if __name__ == '__main__':
    biomart = BioMartSearchEngine()
    #biomart.populateExonDatabase("ENSP00000252816")
    biomart.populateExonDatabase("ENSP00000311134")
    #biomart.populateExonDatabase("ENSP00000341765")
        
        
        
        