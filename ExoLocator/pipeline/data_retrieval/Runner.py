'''
Created on Apr 19, 2012

@author: marioot
'''
from utilities                  import FileUtilities
from pipeline.data_retrieval    import LocalDbSearchEngine, RemoteDbSearchEngine
            
def populate_genes(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['MUTUAL_BEST'] != 'OK':
            print "ABORTING {0}: no description file!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'GENE_RETRIEVAL', 'FAILED')
            continue
        print "POPULATING {0} genes".format(protein_id)
        if LocalDbSearchEngine.populate_sequence_gene(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'GENE_RETRIEVAL', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'GENE_RETRIEVAL', 'PARTIAL') 
            
def populate_expanded_genes(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['MUTUAL_BEST'] != 'OK':
            print "ABORTING {0}: no description file!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'EXP_GENE_RETRIEVAL', 'FAILED')
            continue
        print "POPULATING {0} expanded genes".format(protein_id)
        if LocalDbSearchEngine.populate_sequence_expanded_gene(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'EXP_GENE_RETRIEVAL', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'EXP_GENE_RETRIEVAL', 'PARTIAL')           

def populate_proteins(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['MUTUAL_BEST'] != 'OK':
            print "ABORTING {0}: no description file!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'PROTEIN_RETRIEVAL', 'FAILED')
            continue
        print "POPULATING {0} proteins".format(protein_id)
        if LocalDbSearchEngine.populate_sequence_protein(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'PROTEIN_RETRIEVAL', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'PROTEIN_RETRIEVAL', 'PARTIAL') 

def populate_ensembl_exons(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['MUTUAL_BEST'] != 'OK':
            print "ABORTING {0}: no description file!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'ENSEMBL_EXON_RETRIEVAL', 'FAILED')
            continue
        print "POPULATING {0} ensembl exons".format(protein_id)
        if RemoteDbSearchEngine.populate_sequence_exon_ensembl(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'ENSEMBL_EXON_RETRIEVAL', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'ENSEMBL_EXON_RETRIEVAL', 'PARTIAL') 

def populate_genewise_exons(protein_list):
    pass

def main ():
    protein_list_raw = FileUtilities.get_protein_list()
    protein_list = []
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
        
    populate_genes(protein_list)
    populate_expanded_genes(protein_list)
    populate_proteins(protein_list)
    populate_ensembl_exons(protein_list)
    
if __name__ == '__main__':
    main()