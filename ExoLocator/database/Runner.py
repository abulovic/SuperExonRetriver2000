'''
Created on Jun 25, 2012

@author: marioot
'''
from data_analysis.translation import Runner
from database import DatabasePopulation

def main():
    #Runner.fill_all_containers(True)
    Runner.main()
    
    DatabasePopulation.populate_exon_alignment_piece_table()
    DatabasePopulation.populate_gene_table()
    DatabasePopulation.populate_protein_table()
    DatabasePopulation.populate_exon_table()  #This populates ALIGNMENT table, also!
    DatabasePopulation.populate_ortholog_table()
    
if __name__ == '__main__':
    main()