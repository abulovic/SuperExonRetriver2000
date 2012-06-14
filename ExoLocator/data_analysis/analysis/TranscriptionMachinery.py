'''
Created on Mar 22, 2012

@author: marioot
'''

# Python imports
import re, sys, math, os
from timeit                                         import itertools

# BioPython imports
from Bio.Seq                                        import Seq
from Bio.Alphabet                                   import IUPAC
from Bio                                            import SeqIO
from Bio.SeqRecord                                  import SeqRecord

# utilities imports
from utilities.DescriptionParser                    import DescriptionParser
from utilities.FileUtilities                        import check_status_file, get_species_list, write_failed_species_to_status
from utilities.Logger                               import Logger
from utilities.DirectoryCrawler                     import DirectoryCrawler

# pipeline utilities imports
from pipeline.utilities.AlignmentCommandGenerator   import AlignmentCommandGenerator

# Data analysis imports
from data_analysis.utilities.ExonUtils              import LongestCommonSubstring

from data_analysis.analysis.Exon_translation        import Exon_translation

from data_analysis.base.TranslatedProtein                   import TranslatedProtein
from data_analysis.containers.TranslatedProteinContainer    import TranslatedProteinContainer
from data_analysis.containers.ExonContainer                 import ExonContainer
from data_analysis.containers.EnsemblExonContainer          import EnsemblExonContainer
from data_analysis.containers.ProteinContainer              import ProteinContainer
from data_analysis.containers.DataMapContainer              import DataMapContainer



# to our future selves: we truly and honestly apologize. And again.
ignored_nucleotide_count = 0


def assemble_and_store_protein(protein_id, species, exons_for_transcription, target_prot_seq, assembled_protein_fasta):
    '''
    Assembles the protein from the sw_gene alignment, creates the fasta 
    in the appropriate directory (sequence/assembled_protein)
    
    @param exons_for_transcription: exons that are suitable for transcription (correct order ecc)
    @param target_prot_seq: sequence of the reference species protein
    @param assembled_protein_fasta: path to the assembled protein fasta file
    '''
    
    tpc                 = TranslatedProteinContainer.Instance()

    # actual translation
    (resulting_protein_seq, exon_translations) = transcribe_exons(exons_for_transcription, target_prot_seq)
    
    # write the protein to the appropriate fasta file
    assembled_protein_seq = SeqRecord(Seq(resulting_protein_seq), species, description = "")
    fasta_file = open(assembled_protein_fasta, "w")
    SeqIO.write(assembled_protein_seq, fasta_file, "fasta")
    fasta_file.close()
    
    translated_protein = TranslatedProtein(protein_id, species)
    tpc.add(protein_id, species, translated_protein)
    
    


def create_protein_alignment(protein_id, species):
    '''
    Generates the SW alignment of three protein sequences:
    reference species protein, the assembled protein and the ensembl species protein
    @param protein_id: referent protein id
    @param species: species (latin)
    '''
    
    sequences_for_fasta = []
    
    dc                  = DirectoryCrawler()
    pc                  = ProteinContainer.Instance()
    dmc                 = DataMapContainer.Instance()
    acg                 = AlignmentCommandGenerator()
    tpc                 = TranslatedProteinContainer.Instance()
    
    data_map            = dmc.get((protein_id, species))
    
    # get all the proteins
    ref_protein         = pc.get(protein_id)
    species_protein     = pc.get(data_map.protein_id)
    assembled_protein   = tpc.get(protein_id, species)
    
    sequences_for_fasta.append(ref_protein.get_sequence_record())
    sequences_for_fasta.append(assembled_protein.get_sequence_record())
    sequences_for_fasta.append(species_protein.get_sequence_record())
    
    msa_fasta       = "%s/%s.fa" % (dc.get_mafft_path(protein_id), species)
    msa_afa         = "%s/%s.afa" % (dc.get_mafft_path(protein_id), species)
    msa_fasta_file  = open(msa_fasta, "w")
    SeqIO.write(sequences_for_fasta, msa_fasta_file, "fasta")
    msa_fasta_file.close()
    
    mafft_cmd = acg.generate_mafft_command(msa_fasta, msa_afa)
    os.system(mafft_cmd)


def translate_alignment_exons_for_protein(protein_id, exon_number):
    '''
    Translates all the proteins for which there is SW to gene alignment
    '''
    algorithm = "sw_gene"
    
    # instantiate all the utilities
    logger              = Logger.Instance()
    dc                  = DirectoryCrawler()
    translation_logger  = logger.get_logger("translator")
    
    # instantiate all the containters
    eec                 = EnsemblExonContainer.Instance()
    ec                  = ExonContainer.Instance()
    pc                  = ProteinContainer.Instance()    

    failed_species          = []
    assembled_protein_path  = dc.get_assembled_protein_path(protein_id)

    # for all the species for which it is required to generate translated protein
    for species in get_species_list(protein_id, assembled_protein_path):
        
        # get all you need for the processing
        assembled_protein_fasta = "%s/%s.fa" % (dc.get_assembled_protein_path(protein_id), species)
        exon_key                = (protein_id, species, algorithm)
        target_prot             = pc.get(protein_id)
        target_prot_seq         = target_prot.get_sequence_record().seq
        
        try:
            exons = ec.get(exon_key)
        except KeyError:
            translation_logger.error("%s,%s,%s" % (protein_id, species, "No exons available"))
            failed_species.append(species)
            continue
        exons_for_transcription = []

        # THIS PART WILL NOT EXIST IN THE NEAR FUTURE
        last_translated_exon = False
        for al_exon in exons.get_ordered_exons():

            ref_exon     = eec.get(al_exon.ref_exon_id)
            trans_exon   = Exon_translation(ref_exon, al_exon)
            # if we've already bumped into exon with UTR on its end, all the other exons are not viable
            if last_translated_exon:
                trans_exon.viability = False
                
            if trans_exon.viability:
                (trans_exon, last_translated_exon)  = chop_off_start_utr(al_exon.ref_exon_id, trans_exon, target_prot_seq, exon_number)
                trans_exon                          = chop_off_end_utr (al_exon.ref_exon_id, trans_exon, target_prot_seq, exon_number, protein_id)
            
            exons_for_transcription.append(trans_exon)
        # up to here - this will get trashed
        
        assemble_and_store_protein (protein_id, species, exons_for_transcription, target_prot_seq, assembled_protein_fasta)
        create_protein_alignment   (protein_id, species)
                   
    write_failed_species_to_status(failed_species, assembled_protein_path)
    return failed_species


def get_exon(exons, exon_id):
    for exon in exons:
        if exon.id is exon_id:
            return exon
        
def set_exon(exons, exon):
    for exon_id in exons:
        if exon.id is exon_id:
            exon_id = exon
    return exons

def find_frame(coding_dna, target_transcript):
    
    
    # if the coding dna is 4 bases or shorter, 
    # there is no guarantee that the translation will be valid
    if len(coding_dna) <= 4:
        return [-1, "", target_transcript]
    
    printin = False
    if printin:
        print "Finding frame..."
        print coding_dna
        print target_transcript
    dna = Seq(coding_dna, IUPAC.ambiguous_dna)
    closest_match_on_transcript = len(target_transcript)
    closest_match_frame = -1
    closest_match_translation = ""
    for frame in range(3):
        translation             = dna[frame:].translate().tostring()
        location_on_transcript  = target_transcript.find(translation)
        if location_on_transcript is not -1 and location_on_transcript < closest_match_on_transcript and len(translation):
            closest_match_on_transcript = location_on_transcript
            closest_match_frame = frame
            closest_match_translation = translation
    if closest_match_on_transcript == len(target_transcript):
        raise Exception('No viable frames!')
    if printin:
        print "Left part: %s" % closest_match_translation
        print "Right part: %s" % target_transcript[closest_match_on_transcript + len(closest_match_translation):]
    return [closest_match_frame, 
            closest_match_translation, 
            target_transcript[closest_match_on_transcript + len(closest_match_translation):]]

def indentify_initial_exon_frames(exons, target_prot):
    for exon in exons:
        if exon.viability:
            frame = find_frame(exon.target.replace("-", ""), 
                               target_prot)[0]
            exon.set_frame(frame)
    return exons

def repair_broken_exons(exons, viable_exons):
    for i in range(len(viable_exons) - 1):
        if viable_exons[i + 1] == viable_exons[i] + 1:
            exon_ = get_exon(exons, viable_exons[i])
            _exon = get_exon(exons, viable_exons[i + 1])
            if exon_.t_end == exon_.length and _exon.t_start == 1 and _exon.frame != 0:
                query_nucleotides   =   exon_.query[len(exon_.query) - (3 - _exon.frame) : ]
                target_nucleotides  =   exon_.target[len(exon_.target) - (3 - _exon.frame) : ]
                
                _exon.query         =   query_nucleotides + _exon.query
                _exon.target        =   target_nucleotides + _exon.target
                
                exon_.query         =   exon_.query[ : len(exon_.query) - (3 - _exon.frame)]
                exon_.target        =   exon_.target[ : len(exon_.target) - (3 - _exon.frame)]
                
                exons               =   set_exon(exons, exon_)
                exons               =   set_exon(exons, _exon)
           
def translate_dna(coding_dna, frame):
    dna = Seq(coding_dna, IUPAC.ambiguous_dna)
    return dna[frame:].translate().tostring()

def transcribe_exons(exons, target_prot):
    query_prot      = ""
    viable_exons    = []
    exon_translations = []
    for exon in exons:
        if exon.viability:
            viable_exons.append(exon.id)
    
    indentify_initial_exon_frames(exons, 
                                  target_prot)
    viable_exons = sorted(viable_exons)
    #repair_broken_exons(exons,
    #                    viable_exons)
    for exon in exons:
        if exon.viability:
            exon_to_protein = analyse_SW_alignment(exon.query, 
                                               exon.target, 
                                               target_prot)
            query_prot += exon_to_protein
            exon_translations.append((exon.id, exon_to_protein))
            #print query_prot
    return query_prot, exon_translations

def get_to_next_valid_position(dna, position):
    while dna[position] is '-':
        position    += 1
        if position == len(dna): 
            return -1
    return position

def round_up_to_triplet(position):
    if (position + 1) % 3 == 0:
        return position + 1
    else:
        return position + 2

def count_initial_gaps(target_dna):
    cnt = 0
    while target_dna[cnt] is '-':
        cnt += 1
    return cnt
    
def analyse_query_dna(query_dna, target_dna, target_prot):
    
    global ignored_nucleotide_count
    
    query_prot  =   ""
    
    next_gap    =   query_dna.find('-')
    if next_gap is not -1:
        [frame, left_target_prot, right_target_prot] = find_frame(target_dna[:next_gap], 
                                                                  target_prot) 
        
        if frame == -1:
            # cdna is 4 bases or shorter
            ignored_nucleotide_count += len(target_dna[:next_gap])
        else:
            # calculate number of skipped AAs
            protein_x_count =  int(math.ceil(float(ignored_nucleotide_count + frame)/ 3))
            query_prot += "X"*protein_x_count
            # actually skip those AAs
            target_prot = target_prot[protein_x_count:]
            # reset the ignored base count
            ignored_nucleotide_count = 0

            query_prot          +=  translate_dna(query_dna[:next_gap], 
                                                  frame)
            next_valid_position =   get_to_next_valid_position(query_dna, 
                                                               next_gap)
            number_of_gaps      =   next_valid_position - next_gap
            untranslated_tail   =   len(query_dna[frame:next_gap])%3
            if (untranslated_tail > 0) and (untranslated_tail + number_of_gaps) % 3 is 0:
                query_prot      += "X"
                
            query_prot          +=  analyse_query_dna(query_dna[next_valid_position:], 
                                                      target_dna[next_valid_position:], 
                                                      right_target_prot)
            query_dna           =   query_dna[:next_gap]
    else:
        frame               =   find_frame(target_dna, 
                                           target_prot)[0]
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  translate_dna(query_dna, 
                                              frame)
    return query_prot
    
def analyse_SW_alignment(query_dna, target_dna, target_prot):
    query_prot  =   ""
    init_gaps   =   count_initial_gaps(target_dna)
    insertion   =   translate_dna(query_dna[0:init_gaps], 0)
    query_prot  +=  insertion
    query_dna   =   query_dna[init_gaps:]
    target_dna  =   target_dna[init_gaps:]
    next_gap    =   target_dna.find('-')
    if next_gap is not -1:
        #not the end of sequence
        [frame, left_target_prot, right_target_prot] = find_frame(target_dna[:next_gap], 
                                                                  target_prot)        
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  analyse_query_dna(query_dna[frame:next_gap], 
                                                  target_dna[frame:next_gap], 
                                                  left_target_prot)
        next_valid_position =   get_to_next_valid_position(target_dna, 
                                                           next_gap)
        number_of_gaps      =   next_valid_position - next_gap
        insertion_size      =   number_of_gaps/3
        if insertion_size > 0:
            #print "Insertion %d" % insertion_size
            #print "Range: %d - %d" % (next_gap, next_gap+3*insertion_size)
            #print query_dna[next_gap:next_gap + 3*insertion_size]
            #print query_dna
            #print target_dna
            query_prot      +=  translate_dna(query_dna[next_gap:next_gap + 3*insertion_size], 0).lower()
        query_prot          +=  analyse_SW_alignment(query_dna[next_valid_position:], 
                                                     target_dna[next_valid_position:], 
                                                     right_target_prot)
    else:
        #end of the sequence
        #print "No next gap:"
        #print query_dna
        #print target_dna
        [frame, left_target_prot, target_prot] = find_frame(target_dna, 
                                                            target_prot)        
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  analyse_query_dna(query_dna[frame:], 
                                                  target_dna[frame:], 
                                                  left_target_prot)
        #check to see if the last segment ends as a uncomplete triplet
        if len(target_dna[frame:])%3 > 0:
            query_prot      +=  "X"
    return query_prot

def isSorted(inList):
    for i in range(len(inList) - 1):
        if inList[i] > inList[i+1]: return False
    return True     

def calculate_total_score(valid_id_list, exons):
    score = 0.0
    for exon in exons:
        if exon.id in valid_id_list:
            score += exon.fitness()
    return score

def exaustive_search(exons):
    full_array          = []
    highest_score       = 0.0
    best_combination    = []
    for exon in exons:
        full_array.append(exon.id)
    for i in range(len(full_array)):  
        all_combinations                = itertools.combinations(full_array, len(full_array) - i)
        for comb in all_combinations:
            if isSorted(comb): 
                score                   = calculate_total_score(comb, exons)
                if score > highest_score: 
                    highest_score       = score
                    best_combination    = comb
        if highest_score > 0.0: return best_combination
            
def find_candidates_for_NW(exons):
    # move somewhere else!!!
    treshold_matches = 0.7
    treshold_alignment = 0.8
    #
    candidates = []
    for exon in exons:
        if not exon.viability:
            continue
        percentage_matching     =   float(exon.no_of_matches) / exon.alignment_length
        percentage_aligned      =   float(exon.alignment_length) / exon.length
        if percentage_matching > treshold_matches and percentage_aligned < treshold_alignment:
            candidates.append(exon)
    return candidates

def find_ref_exon_translation(ref_exon_id, target_protein):
    '''
    Auxiliary function, serves to check if the translation of the first
    exon alignment is located in the actual first exon translation
    '''

    eec = EnsemblExonContainer.Instance()
    ref_exon = eec.get(ref_exon_id)
    #print ref_exon.exon_id, ref_exon.sequence
    #ref_exon_seq = Seq (ref_exon.sequence.tostring(), IUPAC.unambiguous_dna)
    
    longest_translation = ""
    
    for frame in range (0,3):
        first_exon_translation = LongestCommonSubstring(ref_exon.sequence[frame:].translate(), target_protein)
        if len(first_exon_translation) > len(longest_translation):
            longest_translation = first_exon_translation

    return longest_translation


def chop_off_start_utr (ref_exon_id, exon, target_protein, number_of_exons):
    '''
    Removes the untranslated region from the aligned exon
    If the exon is also the last exon that gets translated, the last_exon marker is returned True
    If the exon is not translated at all, it is marked as not viable
    '''
    exon_dna = Seq(exon.target.replace("-", ""), IUPAC.ambiguous_dna)
    longest_translation = ""
    longest_translation_frame = -1
    last_exon = False
    
    first_exon_translation  = find_ref_exon_translation (ref_exon_id, target_protein)
    
    # if it's the last exon, we might have the end UTR, 
    # and should not mark it false if the longest
    # translation is not in the end of cdna translation
    if target_protein.endswith(first_exon_translation):
        last_exon = True
    
    for frame in range(0,3):
        translation = exon_dna[frame:].translate().tostring()
        translation = LongestCommonSubstring(translation, first_exon_translation)
        if len(translation) > len(longest_translation):
            longest_translation = translation
            longest_translation_frame = frame
    
    # check if this is truly the beginning of the protein
    # this check shouln't be so rigorous  
    cdna_translation = exon_dna[longest_translation_frame:].translate()
    
    # there is no UTR    
    if first_exon_translation.find(cdna_translation) != -1:
        return (exon, last_exon)
    
    #if exon.id != 1 and not last_exon:
    #    raise Exception ("Not first exon, and has UTR")
          
    location_on_target_prot = first_exon_translation.find(longest_translation)
    if location_on_target_prot == -1:
        exon.viability = False
        return (exon, last_exon)

    
    location_on_cdna_translation = cdna_translation.find(longest_translation)
    if (location_on_cdna_translation + len(longest_translation)) != len(cdna_translation) and exon.id != number_of_exons and not last_exon:
        exon.viability = False
        return (exon, last_exon)
    if location_on_target_prot > len(first_exon_translation)-len(longest_translation):
        exon.viability = False
        #print "Exon not viable?"
        
    else:
        # locate the translated region on the pure region translation
        location_on_translated_pure_dna = exon_dna[longest_translation_frame:].translate().tostring().find(longest_translation)
        location_on_pure_dna_start = location_on_translated_pure_dna * 3 + longest_translation_frame
        location_on_pure_dna_stop  = (location_on_translated_pure_dna + len(longest_translation)) * 3 + longest_translation_frame 
        # extract the pure coding dna
        coding_dna = exon_dna[location_on_pure_dna_start:location_on_pure_dna_stop]
        
        # create regex to match the newly found region
        regex = ".*("
        for char in coding_dna:
            regex += "%s-*" % char
        regex = regex[0:len(regex)-2] + ").*"
        pattern = re.compile(regex)
        pattern_match = re.match(pattern, exon.target)
        
        actual_cdna = pattern_match.groups()[0]
        actual_start = exon.target.find (actual_cdna)
        #print "CDNA 1ST:", actual_cdna
        #actual_stop = actual_start + len(actual_cdna)
        
        exon.target = exon.target[actual_start:]
        exon.query = exon.query[actual_start:]
        exon.t_start += actual_start
        exon.q_start += actual_start
        
        
    return (exon, last_exon)

def chop_off_end_utr (ref_exon_id, exon, target_prot_seq, number_of_exons, protein_id):
    '''
    Removes the ending untranslated region from the exon
    '''

    # first, find the actual translation    
    actual_translation = find_ref_exon_translation(ref_exon_id, target_prot_seq)
    
    # get the '-' free sequence
    exon_dna = Seq(exon.target.replace("-", ""), IUPAC.ambiguous_dna)
    # find the longest translation
    longest_translation = ""
    longest_translation_frame = -1
    
    for frame in range (0,3):
        translation = exon_dna[frame:].translate()
        translation = LongestCommonSubstring(translation, target_prot_seq)
        
        if len(translation) > len(longest_translation):
            longest_translation = translation
            longest_translation_frame = frame
            
    # if this longest_translation isn't in the actual translation, mark the exon non-viable
    cdna_translation = exon_dna[longest_translation_frame:].translate()
    
    # if there is no UTR
    if actual_translation.find(cdna_translation) != -1:
        return exon
    
    if cdna_translation.find(longest_translation) != 0:
        exon.viability = False
        return exon
    
    if actual_translation.find(longest_translation) == -1:
        print "Exon not viable?"
        exon.viability = False
        
    else:
        #if exon.id != number_of_exons:
        #    raise Exception ("Protein %s has UTR" % (protein_id))
        # chop off the final part:
        location_of_actual_translation_start = actual_translation.find(longest_translation)
        # we accidentally found some part of the protein in the untranslated region:
        if (location_of_actual_translation_start * 3 + longest_translation_frame) > len(actual_translation)*3:
            print "TROLOLOLOLO"
        coding_dna = exon_dna[:len(longest_translation)*3+longest_translation_frame]
       
        # create regex to find the position in the original sequence (infested with '-')
        regex = ".*("
        for char in coding_dna:
            regex += "%s-*" % char
        regex = regex[0:len(regex)-2] + ").*"
        pattern = re.compile(regex)
        pattern_match = re.match (pattern, exon.target)
        actual_cdna = pattern_match.groups()[0]
        
        actual_start = exon.target.find(actual_cdna)
        
        actual_stop = actual_start + len(actual_cdna)
        exon.target = actual_cdna
        exon.query = exon.query[actual_start:actual_stop]
        exon.t_start  += actual_start
        exon.t_stop = exon.t_start + (actual_stop - actual_start)
        exon.q_start += actual_start
        exon.q_stop = exon.q_start + (actual_stop - actual_start)    
      
    return exon
        
   
def translate_ensembl_exons(protein_list):
    
    ec = ExonContainer.Instance()
    pc = ProteinContainer.Instance()
    dmc = DataMapContainer.Instance()
    
    for protein_id in protein_list:
        
        if not check_status_file(protein_id):
            continue
        
        species_list = DescriptionParser().get_species(protein_id)
        for species in species_list:
            
            data_map = dmc.get((protein_id, species))
            species_protein = pc.get(data_map.protein_id)
            species_protein_seq = species_protein.get_sequence_record() 
            
            exon_key = (protein_id, species, "ensembl")
            try:
                exons = ec.get(exon_key)
            except Exception:
                continue
            (cdna, locations) = exons.get_cDNA()
            translation_len = 0
            for frame in range(0,3): 
                translated_protein = cdna[frame:].seq.translate()
                common_translation = LongestCommonSubstring(translated_protein, species_protein_seq.seq)
                if len(common_translation) > translation_len:
                    longest_common_translation = common_translation
                    translation_len = len(common_translation)
                
            if not str(longest_common_translation) == str(species_protein_seq.seq):
                print "not OK"
                print "Original:   " + species_protein_seq.seq
                print "Translated: " + longest_common_translation 
            
       
          

#Kao command line argument mu se predaje path do SW outfilea

def main(argv=None):
    if argv is None:
        argv = sys.argv
    sw_output = argv[1]
    exons_SW = parse_SW_output(sw_output)
    exons_SW = discard_FP(exons_SW)
    exons_NW = find_candidates_for_NW(exons_SW)
    
    target_prot = "MFSALKKLVGSDQAPGRDKNIPAGLQSMNQALQRRFAKGVQYNMKIVIRGDRNTGKTALWHRLQGRPFVEEYIPTQEIQVTSIHWSYKTTDDIVKVEVWDVVDKGKCKKRGDGLKMENDPQEAESEMALDAEFLDVYKNCNGVVMMFDITKQWTFNYILRELPKVPTHVPVCVLGNYRDMGEHRVILPDDVRDFIDNLDRPPGSSYFRYAESSMKNSFGLKYLHKFFNIPFLQLQRETLLRQLETNQLDMDATLEELSVQQETEDQNYGIFLEMMEARSRGHASPLAANGQSPSPGSQSPVVPAGAVSTGSSSPGTPQPAPQLPLNAAPPSSVPPVPPSEALPPPACPSAPAPRRSIISRLFGTSPATEAAPPPPEPVPAAEGPATVQSVEDFVPDDRLDRSFLEDTTPARDEKKVGAKAAQQDSDSDGEALGGNPMVAGFQDDVDLEDQPRGSPPLPAGPVPSQDITLSSEEEAEVAAPTKGPAPAPQQCSEPETKWSSIPASKPRRGTAPTRTAAPPWPGGVSVRTGPEKRSSTRPPAEMEPGKGEQASSSESDPEGPIAAQMLSFVMDDPDFESEGSDTQRRADDFPVRDDPSDVTDEDEGPAEPPPPPKLPLPAFRLKNDSDLFGLGLEEAGPKESSEEGKEGKTPSKEKKKKKKKGKEEEEKAAKKKSKHKKSKDKEEGKEERRRRQQRPPRSRERTAADELEAFLGGGAPGGRHPGGGDYEEL"
    print transcribe_exons(exons_SW, target_prot)
   
    return exons_SW
    

if __name__ == '__main__':
    sys.exit(main())



''' 
query_dna = "GTCCTCCGCAAAGGCCTCAAGGCCACATCGGGGCGCAGCTCCCA---GGACCACAGAGCCCCTCTGG------------TCGGGTGGCAC---------CAAGCCCCCCACT--GAGGGCTCCTCCCGAGGGCACGAGGACAGGAGGGACAAGCAGGAGTCCTCA---GAGAGCGACCCCGAGGGGCCCATTGCCGCCCAGATGCTGTCCTTTGTCATGGACGACCCTGACTTTGAGAGCGAC---TCAGATACTCAGCGGACAGCG"
target_dna = "GTCCTCC---CCAGCTTCGAAGCCACGGAGGGGGACAGCTCCCACGAGGACCGCAGCACCCCCCTGGCCAGGCGGTGTCTCTGTTCGCACAGGTCCGGAGAAGCGCAGCAGCACCAGGCCCCCTGCTGAG--ATGGAGCCGGGGAAGGGTGAGCAGGCCTCCTCGTCGGAGAGTGACCCCGAGGGACCCATTGCTGCACAAATGCTGTCCTTCGTCATGGATGACCCCGACTTTGAGAGCGAGGGATCAGACACACAGCGCAGGGCG"
target_prot = "MFSALKKLVGSDQAPGRDKNIPAGLQSMNQALQRRFAKGVQYNMKIVIRGDRNTGKTALWHRLQGRPFVEEYIPTQEIQVTSIHWSYKTTDDIVKVEVWDVVDKGKCKKRGDGLKMENDPQEAESEMALDAEFLDVYKNCNGVVMMFDITKQWTFNYILRELPKVPTHVPVCVLGNYRDMGEHRVILPDDVRDFIDNLDRPPGSSYFRYAESSMKNSFGLKYLHKFFNIPFLQLQRETLLRQLETNQLDMDATLEELSVQQETEDQNYGIFLEMMEARSRGHASPLAANGQSPSPGSQSPVVPAGAVSTGSSSPGTPQPAPQLPLNAAPPSSVPPVPPSEALPPPACPSAPAPRRSIISRLFGTSPATEAAPPPPEPVPAAEGPATVQSVEDFVPDDRLDRSFLEDTTPARDEKKVGAKAAQQDSDSDGEALGGNPMVAGFQDDVDLEDQPRGSPPLPAGPVPSQDITLSSEEEAEVAAPTKGPAPAPQQCSEPETKWSSIPASKPRRGTAPTRTAAPPWPGGVSVRTGPEKRSSTRPPAEMEPGKGEQASSSESDPEGPIAAQMLSFVMDDPDFESEGSDTQRRADDFPVRDDPSDVTDEDEGPAEPPPPPKLPLPAFRLKNDSDLFGLGLEEAGPKESSEEGKEGKTPSKEKKKKKKKGKEEEEKAAKKKSKHKKSKDKEEGKEERRRRQQRPPRSRERTAADELEAFLGGGAPGGRHPGGGDYEEL"

#print analyse_SW_alignment(query_dna, target_dna, target_prot)
exons_SW = parse_SW_output("/home/marioot/Desktop/working_dir/pig_test/SW/Sus_scrofa/sus_scrofa_out.txt")
exons_SW = discard_FP(exons_SW)
print transcribe_exons(exons_SW, target_prot)
'''

        
if __name__ == '__main__':   
    query_dna = "GTCCTCCGCAAAGGCCTCAAGGCCACATCGGGGCGCAGCTCCCA---GGACCACAGAGCCCCTCTGG------------TCGGGTGGCAC---------CAAGCCCCCCACT--GAGGGCTCCTCCCGAGGGCACGAGGACAGGAGGGACAAGCAGGAGTCCTCA---GAGAGCGACCCCGAGGGGCCCATTGCCGCCCAGATGCTGTCCTTTGTCATGGACGACCCTGACTTTGAGAGCGAC---TCAGATACTCAGCGGACAGCG"
    target_dna = "GTCCTCCATACCAGCTTCGAAGCCACGGAGGGGGACAGCTCCCACGAGGACCGCAGCACCCCCCTGGCCAGGCGGTGTCTCTGTTCGCACAGGTCCGGAGAAGCGCAGCAGCACCAGGCCCCCTGCTGAG--ATGGAGCCGGGGAAGGGTGAGCAGGCCTCCTCGTCGGAGAGTGACCCCGAGGGACCCATTGCTGCACAAATGCTGTCCTTCGTCATGGATGACCCCGACTTTGAGAGCGAGGGATCAGACACACAGCGCAGGGCG"
    target_prot = "MFSALKKLVGSDQAPGRDKNIPAGLQSMNQALQRRFAKGVQYNMKIVIRGDRNTGKTALWHRLQGRPFVEEYIPTQEIQVTSIHWSYKTTDDIVKVEVWDVVDKGKCKKRGDGLKMENDPQEAESEMALDAEFLDVYKNCNGVVMMFDITKQWTFNYILRELPKVPTHVPVCVLGNYRDMGEHRVILPDDVRDFIDNLDRPPGSSYFRYAESSMKNSFGLKYLHKFFNIPFLQLQRETLLRQLETNQLDMDATLEELSVQQETEDQNYGIFLEMMEARSRGHASPLAANGQSPSPGSQSPVVPAGAVSTGSSSPGTPQPAPQLPLNAAPPSSVPPVPPSEALPPPACPSAPAPRRSIISRLFGTSPATEAAPPPPEPVPAAEGPATVQSVEDFVPDDRLDRSFLEDTTPARDEKKVGAKAAQQDSDSDGEALGGNPMVAGFQDDVDLEDQPRGSPPLPAGPVPSQDITLSSEEEAEVAAPTKGPAPAPQQCSEPETKWSSIPASKPRRGTAPTRTAAPPWPGGVSVRTGPEKRSSTRPPAEMEPGKGEQASSSESDPEGPIAAQMLSFVMDDPDFESEGSDTQRRADDFPVRDDPSDVTDEDEGPAEPPPPPKLPLPAFRLKNDSDLFGLGLEEAGPKESSEEGKEGKTPSKEKKKKKKKGKEEEEKAAKKKSKHKKSKDKEEGKEERRRRQQRPPRSRERTAADELEAFLGGGAPGGRHPGGGDYEEL"
    #print analyse_SW_alignment(query_dna, target_dna, target_prot)
    ec = ExonContainer.Instance()
    eec = EnsemblExonContainer.Instance()
    exon_key = ("ENSP00000340983", "Ailuropoda_melanoleuca", "sw_gene")
    exons = ec.get(exon_key)
    exons_for_transcription = []
    for al_exon in exons.get_ordered_exons():
        al_exon = al_exon[0]
        ref_exon = eec.get(al_exon.ref_exon_id)
        trans_exon = Exon_translation(al_exon.ordinal, 
                                      ref_exon.length, 
                                      al_exon.alignment_info["query_seq"], 
                                      al_exon.alignment_info["sbjct_seq"])
        trans_exon.set_intervals(al_exon.alignment_info["query_start"], 
                                 al_exon.alignment_info["query_end"], 
                                 al_exon.alignment_info["sbjct_start"], 
                                 al_exon.alignment_info["sbjct_end"])
        trans_exon.set_identity(al_exon.alignment_info["identities"], al_exon.alignment_info["length"])
        trans_exon.set_viablity(al_exon.viability)
        exons_for_transcription.append(trans_exon)
        
    print transcribe_exons(exons_for_transcription, target_prot)
    '''
    #print analyse_SW_alignment(query_dna, target_dna, target_prot)
    exons_SW = parse_SW_output("/home/intern/Documents/sw_input/sw_output/ex1.swout")
    for ex in exons_SW:
        print ex.match, ex.length, ex.score
    print transcribe_exons(exons_SW, target_prot)
    '''