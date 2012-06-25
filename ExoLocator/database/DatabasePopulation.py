'''
Created on Jun 25, 2012

@author: marioot
'''
from utilities.ConfigurationReader import Singleton, ConfigurationReader
import MySQLdb
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.utilities import generate_structure
from data_analysis.containers.DataMapContainer import DataMapContainer
from utilities import FileUtilities
from utilities.DescriptionParser import DescriptionParser
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.base.EnsemblExon import EnsemblExon
from data_analysis.base.GenewiseExon import GenewiseExon
from data_analysis.analysis.AlignmentPostprocessing import annotate_spurious_alignments_batch
from data_analysis.base.Exon import Exon
from data_analysis.containers.GeneContainer import GeneContainer
from data_analysis.translation import BestExonAlignment, Runner
from data_analysis.translation.AlignmentExonPiece import AlignmentExonPiece
from data_analysis.translation.BestExonAlignmentContainer import BestExonAlignmentContainer
import data_analysis
from database.DBManager import DBManager

def populate_gene_table():
    dbm = DBManager.Instance()
    dmc = DataMapContainer.Instance()
    
    protein_id_list = FileUtilities.get_protein_list()
    species_list = FileUtilities.get_default_species_list()
    data_map_list = []
    for ref_protein_id in protein_id_list:
        for species in species_list:
            try:
                data_map = dmc.get((ref_protein_id[0], species))
                data_map_list.append(data_map)
            except KeyError:
                print "PROTEIN_ID %s ERROR" % (ref_protein_id[0])
    print data_map_list
    dbm.update_gene_table(data_map_list)

def populate_exon_table():
    dbm = DBManager.Instance()
    ec = ExonContainer.Instance()
    dmc = DataMapContainer.Instance()
    
    protein_id_list = FileUtilities.get_protein_list()
    species_list = FileUtilities.get_default_species_list()
    exon_type_list = ["ensembl", "genewise", "blastn", "tblastn", "sw_gene"]
    
    exon_list = []
    for ref_protein_id in protein_id_list:
        for species in species_list:
            for exon_type in exon_type_list:
                exon_key = (ref_protein_id[0], species, exon_type)
                try:
                    exons = ec.get(exon_key).get_ordered_exons()
                    for exon in exons:
                        if type(exon) is Exon:
                            if exon.viability:
                                exon_list.append(exon)
                        else:
                            exon_list.append(exon)
                except KeyError:
                    pass
    dbm.update_exon_table(exon_list)
    dbm.update_alignment_table(exon_list)
                   
def populate_exon_alignment_piece_table():
    dbm = DBManager.Instance()
    ec = ExonContainer.Instance()
    beac = BestExonAlignmentContainer.Instance()
    
    protein_id_list = FileUtilities.get_protein_list()
    species_list = FileUtilities.get_default_species_list()
    
    exon_aln_list = []
    for (ref_protein_id, exon_num) in protein_id_list:
        for species in species_list:
                try:
                    ref_exons = ec.get((ref_protein_id, 'Homo_sapiens', 'ensembl'))
                    for ref_exon in ref_exons.get_coding_exons():
                        best_exon_alignment = beac.get(ref_exon.exon_id, species)
                        if best_exon_alignment and best_exon_alignment.sw_gene_alignment:
                            for aln_piece in best_exon_alignment.sw_gene_alignment.alignment_pieces:
                                if aln_piece.type in ('coding', 'insertion'):
                                    exon_aln_list.append([ref_exon.exon_id, species, aln_piece])
                except KeyError, e:
                    print e
    dbm.update_exon_alignment_piece_table(exon_aln_list)             

def populate_protein_table():
    dbm = DBManager.Instance()
    pc = ProteinContainer.Instance()
    dmc = DataMapContainer.Instance()
    
    protein_id_list = FileUtilities.get_protein_list()
    species_list = FileUtilities.get_default_species_list()
    protein_list = []
    for ref_protein_id in protein_id_list:
        for species in species_list:
            try:
                protein_id = dmc.get((ref_protein_id[0], species))
                protein = pc.get(protein_id.protein_id)
                protein_list.append(protein)
            except KeyError:
                print "PROTEIN_ID %s ERROR" % (ref_protein_id[0])
    dbm.update_protein_table(protein_list)

def populate_ortholog_table():
    dbm = DBManager.Instance()
    dmc = DataMapContainer.Instance()
    
    protein_id_list = FileUtilities.get_protein_list()
    species_list = FileUtilities.get_default_species_list()
    data_map_list = []
    for ref_protein_id in protein_id_list:
        for species in species_list:
            try:
                data_map = dmc.get((ref_protein_id[0], species))
                data_map_list.append(data_map)
            except KeyError, e:
                pass
                #print e
                #print "PROTEIN_ID %s ERROR" % (ref_protein_id[0])
    print data_map_list
    dbm.update_ortholog_table(data_map_list)