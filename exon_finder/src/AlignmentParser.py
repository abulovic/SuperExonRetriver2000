'''
Created on Mar 27, 2012

@author: Mario, anana
'''

import ConfigParser
from Exon import Exon

class AlignmentParser(object):
    '''
    Class serves as a parent class for all the alignment parsing classes,
    @see: blastAlignmentParser, @see: swAlignmentParser. 
    '''


    def __init__(self):
        '''
        Constructor. Loads the necessary configuration.
        '''
        config_file         = "../../config.cfg"
        self.config              = ConfigParser.RawConfigParser()
        self.config.read(config_file)
        
        self.readConfigurations()
        
    def readConfigurations (self):
        
        self.projectRoot            = self.config.get('Project root', 'project_root_folder')
        self.sessionsFolder         = "%s/%s" % (self.projectRoot, self.config.get('Project root', 'session_folder'))
        self.geneRegions            = self.config.get('Gene regions path', 'regions')
        self.expandedGeneRegions    = self.config.get('Gene regions path', 'expanded_regions')
        self.exons                  = self.config.get('Exon database path', 'exons_path')
        self.exonsDb                = "%s/%s" % (self.exons, self.config.get('Exon database path', 'exons_db'))
        self.blastout               = self.config.get('Blastout path', 'blastout')
        self.swout                  = self.config.get('SWout path', 'swout')
        
    def setProteinFolder (self, proteinDir):
        '''
        Merges all the session folders to the protein directory
        and generates absolute paths 
        @param proteinDir: Ensembl protein ID (most probably, according to the last configuration)
        '''
        self.geneRegionsAbs             = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.geneRegions)
        self.expandedGeneRegionsAbs     = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.expandedGeneRegions)
        self.exonsAbs                   = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.exons)
        self.exonsDbAbs                 = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.exonsDb)
        self.blastoutAbs                = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.blastout)
        self.swoutAbs                   = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.swout)
        
        
      
    def fetchDNAseq(self, speciesName):
        '''
        Fetches DNA sequence for particular species in one string
        @param speciesName: Ensembl name of species
        @return: DNA sequence in string format
        '''
        fasta_f = self.geneRegionsAbs + speciesName + '.fa'
        DNA_f   = open (fasta_f, 'r');
        DNA_f.readline();
        seq     = DNA_f.read();
        seq     = "".join(seq.split('\n'))
        return seq
    
    def initializeExonDictionary(self, number_of_exons):
        '''
        Create a dictionary of exons.
        This is an auxiliary method.
        '''
        exons = {}
        for current_exon in range(1, number_of_exons + 1):
            exons[current_exon] = Exon(current_exon, 0, "", "")
        return exons
      
    def parseOutput (self, alignmentOutputFile, numberOfExons):
        pass
    
    