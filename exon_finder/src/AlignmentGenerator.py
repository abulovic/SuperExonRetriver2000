'''
Created on Mar 28, 2012

@author: intern
'''

import ConfigParser
import sys, os, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from AlignmentParserSW import AlignmentParserSW

class AlignmentGenerator(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Loads the configurations and sets default reference species.
        '''
        config_file         = "../../config.cfg"
        self.config              = ConfigParser.RawConfigParser()
        self.config.read(config_file)
        
        self.readConfigurations()
        
        self.algorithms             = self.enum('BLASTN', 'TBLASTN', 'SW', 'GENEWISE')
        self.referenceSpeciesAll    = "Homo_sapiens"
        self.referenceSpeciesDict   = None
        
    def enum(self, *sequential, **named):
        '''
        Auxiliary class. Python does not have an enum :(
        '''
        enums = dict(zip(sequential, range(len(sequential))), **named)
        return type('Enum', (), enums)
        
    def readConfigurations (self):
        '''
        This will be transfered into a Singleton class.
        So far so ugly.
        '''
        # alignment algorithms
        self.blastn                 = self.config.get('Blast cfg', 'blastn')
        self.tblastn                = self.config.get('Blast cfg', 'tblastn')
        self.swSharp                = self.config.get('SW#', 'swSharp')
        self.wise                   = self.config.get('Wise cfg', 'wise')
        self.wise_flags             = self.config.get('Wise cfg', 'flags')
        sys.path.append(self.config.get('SW#', 'src'))
        # project session folders
        self.projectRoot            = self.config.get('Project root', 'project_root_folder')
        self.sessionsFolder         = "%s/%s" % (self.projectRoot, self.config.get('Project root', 'session_folder'))
        
        self.geneRegions            = self.config.get('Gene regions path', 'regions')
        self.expandedGeneRegions    = self.config.get('Gene regions path', 'expanded_regions')
        
        self.proteins               = self.config.get('Found proteins path', 'proteins')
        
        self.exons                  = self.config.get('Exon database path', 'exons_path')
        self.exonsDb                = "%s/%s" % (self.exons, self.config.get('Exon database path', 'exons_db'))
        self.exonsWiseout                = "%s/%s" % (self.exons, self.config.get('Exon database path', 'exons_wiseout'))
        
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
        
        self.proteinsAbs                = "%s/%s/%s"  % (self.sessionsFolder, proteinDir, self.proteins)
        
        self.exonsAbs                   = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.exons)
        self.exonsDbAbs                 = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.exonsDb)
        self.exonsWiseoutAbs            = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.exonsWiseout)
        
        self.blastoutAbs                = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.blastout)
        self.swoutAbs                   = "%s/%s/%s/" % (self.sessionsFolder, proteinDir, self.swout) 
        
    def setReferenceSpeciesForAll (self, referenceSpeciesName):
        '''
        Sets the reference species name.
        This should be invoked before any of the alignment algorithms
        because default reference species is Homo_sapiens, which 
        may or may not be good for your species.
        '''
        self.referenceSpeciesAll = referenceSpeciesName
        
    def setReferenceSpecies (self, referenceSpeciesDict):
        '''
        @param referenceSpeciesDict: Dictionary where keys are all the species (latin names),
                                     and values are their reference species (for blast alignments).
        '''
        self.referenceSpeciesDict = referenceSpeciesDict
                
                
    def runBatchAlignmentAlgorithm (self, algorithm):
        '''
        Runs a batch alignment.
        @param algoritm: has to be one of the alignmentGeneratorObject.algorithms.
                        {BLASTN, TBLASTN, SW, GENEWISE}
        '''
        {self.algorithms.BLASTN     : self.runBatchBlastn (),
         self.algorithms.TBLASTN    : self.runBatchTblastn(),
         self.algorithms.SW         : self.runBatchSW(),
         self.algorithms.GENEWISE   : self.runBatchWise()    } [algorithm] ()
         
    def runBatchBlastn (self, expanded = True):
        '''
        Runs the blastn for all the species in the gene_regions directory
        @param expanded: if True, batch test will be run on the expanded
                        gene regions. Default is True
        '''
        
        queryDirectory = ""
        if (expanded is True):
            queryDirectory = self.expandedGeneRegionsAbs
        else:
            queryDirectory = self.geneRegionsAbs
            
        for fastaFile in os.listdir(queryDirectory):
            if (fastaFile.endswith('.fa')):
                species = fastaFile.split('.')[0]
                self.runBlastnForSpecies(species, queryDirectory)
            
    def runBatchTblastn (self):
        '''
        Runs the tblastn for all the species in the gene_regions directory
        @param expanded: if True, batch test will be run on the expanded
                        gene regions. Default is True
        '''
        for fastaFile in os.listdir(self.proteinsAbs):
            if (fastaFile.endswith('.fa')):
                species = fastaFile.split('.')[0]
                self.runTblastnForSpecies(species)
        pass
    
    def runBatchSW(self, swType = "ex-dna", expanded = True):
        '''
        Runs a batch Smith-Waterman alignment for all the species in the 
        source folder. Source folder is dependant of the type of the
        batch query.
        @param type: type of query. There are two types. ex-dna means that
        the human exons will be aligned to species protein coding region,
        and ex-ex means that each of the species exons will be aligned
        with all of the human exons.
        @param expanded: if True, protein coding regions will be expanded
        Default is True.
        '''
        
        if (swType == "ex-dna"):
            if (expanded == True):
                queryDirectory = self.expandedGeneRegionsAbs
            else :
                queryDirectory = self.geneRegionsAbs
        else :
            queryDirectory = self.exonsDbAbs
            
        for fastaFile in os.listdir(queryDirectory) :
            if (fastaFile.endswith('.fa')):
                speciesName = fastaFile.split('.')[0]
                self.runSWForSpecies(speciesName, swType, expanded)
                
        else :
            queryDirectory = self.exonsDbAbs
            for fastaFile in os.listdir(queryDirectory):
                speciesName = fastaFile.split('.')[0]
                self.runSWForSpecies(speciesName, swType="ex-ex")
                

    def runBatchGenewise (self, speciesList, expanded = False):
        
        if (expanded == True):
            dnaFolder = self.expandedGeneRegionsAbs
        else:
            dnaFolder = self.geneRegionsAbs
        
        for species in speciesList :
            proteinFile = "%s/%s.fa" % (self.proteinsAbs, species)
            dnaFile = "%s/%s.fa" % (dnaFolder, species)
            wise_out = "{0}/{1}.{2}".format(self.exonsWiseoutAbs, species, self.config.get('Wise cfg', 'outfile'))
            exons_out = "%s/%s.fa" % (self.exonsDbAbs, species) 
            
            # run wise
            cmd = "{0} {1} {2} {3} > {4}".format(self.wise, proteinFile, dnaFile, self.wise_flags, wise_out)
            print cmd
            os.system(cmd)
            
            self.createExonsDescription(wise_out, dnaFile, exons_out)
    
    
    def runBlastnForSpecies (self, speciesName, queryDirectory=None):
        '''
        Runs blastn for single species
        @param speciesName: Ensembl species name, underscore delimited
        @param queryDirectory: if not given, expanded gene regions will be used
        '''
        # in case query directory not specified, use the expanded regions
        if (queryDirectory is None):
            queryDirectory = self.expandedGeneRegionsAbs
            
        blastoutFile = "%s/dna/%s.blastout" % (self.blastoutAbs, speciesName)
        querySequence = "%s/%s.fa" % (queryDirectory, speciesName)
        
        refSpec = ""
        if (self.referenceSpeciesDict is not None):
            refSpec = self.referenceSpeciesDict[speciesName]
        else:
            refSpec = self.referenceSpeciesAll
            
        exonDbFile = "%s/%s.fa" % (self.exonsDbAbs, refSpec)
        self.checkIfDatabaseExists(exonDbFile, False)
        
        if (not os.path.isfile(querySequence)):
            print querySequence, "file does not exist."
            return
        
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.blastn, exonDbFile, querySequence, blastoutFile)
        print cmd
        os.system(cmd)
         
    def runTblastnForSpecies (self, speciesName):
       
        '''
        Runs tblasn for single species. It alignes species protein to
        a database of human exons.
        @param speciesName: Ensembl species name, underscore delimited
        '''
        refSpec = ""
        if (self.referenceSpeciesDict is not None):
            refSpec = self.referenceSpeciesDict[speciesName]
        else:
            refSpec = self.referenceSpeciesAll
            
        querySequence = "%s/%s.fa" % (self.proteinsAbs, speciesName)
        exonDbFile = "%s/%s.fa" % (self.exonsDbAbs, refSpec)
        blastoutFile = "%s/protein/%s.blastout" % (self.blastoutAbs, speciesName)
        
        self.checkIfDatabaseExists(exonDbFile, False)
        
        cmd = "{0} -d {1} -i {2} -o {3}".format(self.tblastn, exonDbFile, querySequence, blastoutFile)
        print cmd
        os.system(cmd)
        

            
            
    
    def runSWForSpecies (self, speciesName, swType = "ex-dna", expanded = True):
        '''
        Runs Smith-Waterman alignment for single species.
        @param swType: there are two types of supported SW alignment.  
        First is to align species DNA to database of reference species exons.(ex-dna)
        Second is to align each of the exons of the species to database
        of reference species exons. (ex-ex)
        @param expanded: if true, expanded gene region will be used. Default is true.
        '''
        refSpec = ""
        if (self.referenceSpeciesDict is not None):
            refSpec = self.referenceSpeciesDict[speciesName]
        else:
            refSpec = self.referenceSpeciesAll
        
        exonDbFile = "%s/%s.fa" % (self.exonsDbAbs, refSpec)
        unimportant_file = "tmp.txt" # serves just to supress the sw output
            
        if (swType == "ex-dna"):
            if (expanded == True):
                querySequence = "%s/%s.fa" % (self.expandedGeneRegionsAbs, speciesName)
            else:
                querySequence = "%s/%s.fa" % (self.geneRegionsAbs, speciesName)
                
            outputFile = "%s/dna/%s.swout" % (self.swoutAbs, speciesName)
            
            self.checkIfDatabaseExists(exonDbFile, False)
            
            cmd = "{0} -i {1} -j {2} --out {3} > {4}".format(self.swSharp, querySequence, exonDbFile, outputFile, unimportant_file)
            os.system(cmd)
            os.remove(unimportant_file)
            print cmd
            #os.system(cmd)
         
        else:
            exonSpecies = "%s/%s.fa" % (self.exonsDbAbs, speciesName)
            outputDir = "%s/exon/%s" % (self.swoutAbs, speciesName)
            if (not os.path.isdir(outputDir)):
                os.makedirs(outputDir)
            species_exons = self.disectExonsAndRunSW(exonSpecies, exonDbFile, outputDir)


    def disectExonsAndRunSW(self, exonSpeciesFasta, exonDbFile, outputDir):
        '''
        Auxiliary function. Reads the species exons file and saves each exon
        to a temporary file. After that it does a SW  alignment to a database
        of exons and saves the output to a file.
        @param exonSpeciesFasta: fasta file which contains all the known (or
        predicted exons for the species protein)
        @param exonDbFile: fasta file with exons of reference species
        @param outputDir: output directory for alignment results
        '''
        tmpFile = "tmp.fa"
        
        '''
        There are two types of exon formating.
        Ensembl - >ExonID|TranscriptID
        Wiseout - >Exon exonNum
        '''   
        
        exonSpeciesFile = open(exonSpeciesFasta, 'r')
        exonSeq = ""
        exonHeader = ""
        exonCounter = 0
        for line in exonSpeciesFile.readlines():
            if (line.startswith('>')):
                # if exon Sequence available, store it in the tmp file and 
                if (exonSeq != ""):
                    exonCounter = exonCounter + 1
                    outputFile = "%s/exon%d.swout" % (outputDir, exonCounter)
                    
                    tmp = open(tmpFile, 'w')
                    tmp.write("%s\n%s\n" % (exonHeader, exonSeq))
                    tmp.close()
                    
                    cmd = "{0} -i {1} -j {2} --out {3}".format(self.swSharp, tmpFile, exonDbFile, outputFile)
                    print cmd
                    #os.system(cmd)
                    
                exonSeq = ""
                exonHeader = line.strip()
                
            else:
                exonSeq = exonSeq + line.strip()
                
        exonCounter = exonCounter + 1
        outputFile = "%s/exon%d.swout" % (outputDir, exonCounter)
        tmp = open(tmpFile, 'w')
        tmp.write("%s\n%s\n" % (exonHeader, exonSeq))
        tmp.close()
        cmd = "{0} -i {1} -j {2} --out {3}".format(self.swSharp, tmpFile, exonDbFile, outputFile)
        print cmd
        #os.system(cmd)
                
        exonSpeciesFile.close()
        
        
    def checkIfDatabaseExists(self, dbFileName, protein):
        '''
        Check if database is already existent. If no, create it using formatdb tool.
        @param dbFileName: database file
        @param protein: if True, database will be created for proteins and vice versa
        '''
        if (protein == True):
            db1 = "%s.phr" % dbFileName
            db2 = "%s.pin" % dbFileName
            db3 = "%s.psd" % dbFileName
            db4 = "%s.psi" % dbFileName
            db5 = "%s.psq" % dbFileName
            p = "T"
        else:
            db1 = "%s.nhr" % dbFileName
            db2 = "%s.nin" % dbFileName
            db3 = "%s.nsd" % dbFileName
            db4 = "%s.nsi" % dbFileName
            db5 = "%s.nsq" % dbFileName
            p = "F"
            
        # if database does not exist, create it
        if (not (os.path.isfile(db1) and os.path.isfile(db2) and os.path.isfile(db3) and os.path.isfile(db4) and os.path.isfile(db5))):
            '''
            What this mysterious chunk of code does in this:
            - copies the database file to a temporary file
            - reads the temporary file and writes the exons with preceding enumerations to the new file
            - new file has the name of the (old) database file
            - formatdb is called on the new file
            - after these ugly manipulations, the temporary file is copied back to the old database file
            - now we have a database of enumerated exons, but still uniform exon database files
            '''
            '''
            tmpFileName = "%s.tmpFile" % dbFileName
            shutil.copy(dbFileName, tmpFileName)
            
            tmpFile = open(tmpFileName,'r')
            dbFile = open(dbFileName, 'w')
            
            exonCount = 1
            for line in tmpFile.readlines():
                if (line.startswith('>')):
                    
                    dbFile.write(">%d %s" % (exonCount, line[1:]))
                    exonCount = exonCount + 1
                else :
                    dbFile.write(line)
                    
            tmpFile.close()
            dbFile.close()  
            '''
            cmd = "formatdb -i %s -p %s -o F -a F" % (dbFileName, p)
            os.system(cmd)
            #shutil.move(tmpFileName, dbFileName)
            
    def createExonsDescription(self, wise_f, dna_f, exons_out):
        # database config
        
        protein_exon_regions_out = "%s/protein_regions" % (self.exonsAbs)
        
        formatdb = self.config.get('Formatdb cfg', 'formatdb')
        formatdb_flgs = self.config.get('Formatdb cfg', 'dna_flgs')
        # Read DNA file
        input_DNA_f = open(dna_f, 'r');
        input_DNA_f.readline();
        DNA = re.sub("\n", "", input_DNA_f.read());
        
        DNA = Seq( DNA, IUPAC.unambiguous_dna )
        input_DNA_f.close()
        # Read wise2_out file
        wise_txt = open(wise_f, 'r');
        # Create exon file
        exon_outfile = open(exons_out, 'w')
        protein_exon_mapping_file = open(protein_exon_regions_out, 'w')
        # Fill the exon file
        exon_cnt = 0
        regex_exon = re.compile(r'  Exon (\d+) (\d+) phase \d+')
        regex_protein = re.compile(r'     Supporting \d+ \d+ (\d+) (\d+)')
        for line in wise_txt:
            exon_match = re.match(regex_exon, line)
            if exon_match is not None:
                exon_cnt += 1
                lower = int(exon_match.groups()[0]) - 1
                upper = int(exon_match.groups()[1])
                record = SeqRecord(Seq(str(DNA[lower:upper]),
                           IUPAC.unambiguous_dna),
                           id=str(exon_cnt), description="exon length {0}".format(upper - lower))
                SeqIO.write(record, exon_outfile, "fasta")
            protein_match = re.match(regex_protein, line)
            if protein_match is not None:
                lower = int(protein_match.groups()[0])
                upper = int(protein_match.groups()[1])
                protein_exon_mapping_file.write("exon {0} {1} {2}\n".format(str(exon_cnt), lower, upper))
        exon_outfile.close()
        protein_exon_mapping_file.close()
        # Format database
        cmd = "%s -i %s %s" % (formatdb, exons_out, formatdb_flgs)
        os.system(cmd)
        return exon_cnt
        
            
       
        
if __name__ == '__main__':
    alGen = AlignmentGenerator()
    alGen.setProteinFolder("ENSP00000311134")
    #alGen.runBatchBlastn(True)
    #alGen.runBatchTblastn()
    alGen.runBatchSW("ex-ex")


        
        