'''
Created on Mar 29, 2012

@author: intern
'''
import sys, csv, re
import ConfigParser
from operator import attrgetter
from timeit import itertools

class StatisticsGenerator(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        config_file         = "../../config.cfg"
        self.config              = ConfigParser.RawConfigParser()
        self.config.read(config_file)
        self.readConfigurations()
        self.referenceSpeciesAll    = "Homo_sapiens"
        self.referenceSpeciesDict   = None
        
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
        
        self.statistics             = self.config.get('Statistics', 'exon_finder')
        
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
        
        self.statisticsAbs              = "%s/%s/%s" % (self.sessionsFolder, proteinDir, self.statistics)
        
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
        
        
    def generate_statistics(self, exons_via_proteins, exons_via_dna, exons_SW, exonsSWEnsembl, all_species, protein_id):
        statistics_header = ["Protein_ID", 
                             "Species", 
                             "Type_of_search", 
                             "Exon_number", 
                             "Length", 
                             "Score", 
                             "Alignment_matches", 
                             "Alignment_length"]
        tblastn_statistics = []
        blastn_statistics = []
        SW_statistics = []
        SWE_statistics = []
        
        base_exon_length = self.get_exon_lenghts(len(exons_via_proteins))
        
        for species in sorted(all_species):
            if (exons_via_proteins.has_key(species)):
                tblastn_species_statistic = self.generate_statistics_based_on_search(exons_via_proteins[species], base_exon_length)
                for exon in tblastn_species_statistic:
                    exon[3] *= 3    # proteins -> dna
                    exon[4] *= 3
                    exon.insert(0, "tblastn")
                    exon.insert(0, species)
                    exon.insert(0, protein_id)
                tblastn_statistics.append(tblastn_species_statistic)
                
            if (exons_via_dna.has_key(species)):
                blastn_species_statistic = self.generate_statistics_based_on_search(exons_via_dna[species], base_exon_length)
                for exon in blastn_species_statistic:
                    exon.insert(0, "blastn")
                    exon.insert(0, species)
                    exon.insert(0, protein_id)
            
                blastn_statistics.append(blastn_species_statistic)
                
            if (exons_SW.has_key(species)):
                SW_species_statistic = self.generate_statistics_based_on_search(exons_SW[species], base_exon_length)
                for exon in SW_species_statistic:
                    exon.insert(0, "SW")
                    exon.insert(0, species)
                    exon.insert(0, protein_id)
                SW_statistics.append(SW_species_statistic)
                
            if (exonsSWEnsembl.has_key(species)):
                SWE_species_statistic = self.generate_statistics_based_on_search(exons_SW[species], base_exon_length)
                for exon in SWE_species_statistic:
                    exon.insert(0, "SW_ensembl")
                    exon.insert(0, species)
                    exon.insert(0, protein_id)
                SWE_statistics.append(SWE_species_statistic)
            #print generate_statistics_based_on_search(exons_SW[species])
        ###WRITE TO CSV FILE###
        statout = csv.writer(open("%s.csv" % self.statisticsAbs, 'wb+'), delimiter = ',')
        statout.writerow(statistics_header)
        for species in tblastn_statistics:
            for exon in species:
                statout.writerow(exon)
        for species in blastn_statistics:
            for exon in species:
                statout.writerow(exon)
        for species in SW_statistics:
            for exon in species:
                statout.writerow(exon)
        for species in SWE_statistics:
            for exon in species:
                statout.writerow(exon)
        ####
        return [statistics_header, tblastn_statistics, blastn_statistics]
    
    
    def get_exon_lenghts (self, number_of_exons):
        exons       = {}
        exon_db = "%s/db/%s.fa" % (self.exonsAbs, self.referenceSpeciesAll)
        exon_file   = open(exon_db, "r")
        
        p1 = re.compile('>(\d+) exon length (\d+)')
        p2 = re.compile('>(\d+) (\d+)\|(\d+)\|(\w+)\|(\w+)')
        for line in exon_file.readlines():
            m1 = re.match(p1, line)
            m2 = re.match(p2, line)
            if m1 is not None:
                exons[int(m1.groups()[0])] = int(m1.groups()[1])
            if m2 is not None:
                exons[int(m2.groups()[0])] = int(m2.groups()[2]) - int(m2.groups()[1])
        return exons
    
    def generate_statistics_based_on_search (self, exons, base_exon_lengths):
        statistics = []
        
        for current_exon in range(1, len(exons) + 1):
            if exons[current_exon].viability:
                statistics.append([current_exon,
                                   base_exon_lengths[current_exon],
                                   exons[current_exon].score, 
                                   exons[current_exon].no_of_matches, 
                                   exons[current_exon].alignment_length])
            else:
                statistics.append([current_exon,
                                   base_exon_lengths[current_exon],
                                   0, 
                                   0, 
                                   0])
        return statistics
    
if __name__ == '__main__':
    statGen = StatisticsGenerator()
    statGen.setProteinFolder("ENSP00000253237")
    