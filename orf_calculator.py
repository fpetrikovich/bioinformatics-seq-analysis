import math
from Bio import SeqIO
from constants import START_CODON, END_CODONS
from file_helper import save_file
from bio_utils import translate_sequence_to_protein, transcribe_dna, check_if_valid_orf_with_cds, get_cds_protein_from_record

codon_size = 3
orf_count = 6

def run_exercise_1(genbank_file, fasta_nuc_output_file, fasta_prot_output_file):

    nucleotide_handle = open(fasta_nuc_output_file, "w")
    proteins_handle = open(fasta_prot_output_file, "w")

    orf_status = [None] * (orf_count + 1)
    
    orf = 1

    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        print('Working with GenBank sequence record %s' % seq_record.id)

        sequence = transcribe_dna(seq_record.seq)
        reverse_sequence = sequence.reverse_complement()

        for i in range(0, 3):
            orf_status = handle_new_open_reading_frame(orf, sequence, i, False, orf_status, nucleotide_handle, proteins_handle)
            orf += 1
            orf_status = handle_new_open_reading_frame(orf, reverse_sequence, i, True, orf_status, nucleotide_handle, proteins_handle)
            orf += 1

    nucleotide_handle.close()
    proteins_handle.close()


def find_start_stop_codons(nucleotides):
    """
    Finds posible proteins by checking start and stop codons. 
    Start codon: AUG
    Stop codons: UAA, UAG y UGA
    Arguments: 
        nucleotides: string with ARNm nucleotides already transcribed
    Returns: substring of nucleotides between start and stop codons.
    """
    final_start_idx = 0
    final_end_idx = 0
    start_idx = 0
    end_idx = 3
    found_start = False
    coding_sections = []
    
    while end_idx < len(nucleotides):
        codon = nucleotides[start_idx:end_idx]

        if not found_start:
            if START_CODON == codon:
                found_start = True
                final_start_idx = start_idx
        else:
            if codon in END_CODONS:
                final_end_idx = start_idx
                coding_sections.append(nucleotides[final_start_idx:final_end_idx])
                found_start = False

        start_idx = end_idx
        end_idx = end_idx + 3

    max_length = 0
    final_nucleotides = ""
    for sec in coding_sections:
        if len(sec) > max_length:
            max_length = len(sec)
            final_nucleotides = sec

    return final_nucleotides


def handle_new_open_reading_frame(orf_id, nucleotides, start_offset, is_reverse, orf_status, nucleotide_handle, protein_handle):
    # Calculate where the reading frame starts and ends       
    start = start_offset
    end = math.floor((len(nucleotides) - start_offset)/codon_size) * codon_size + start_offset
    # Get the possible protein from analyzing the codons
    final_nucleotides = find_start_stop_codons(nucleotides[start:end])
    # translate sequence to protein for a given ORF
    protein_sequence = translate_sequence_to_protein(final_nucleotides, cds=False)
    # Save the translation to compare with the CDS
    orf_status[orf_id] = {
        "start": start,
        "end": end,
        "is_reverse": is_reverse,
        "protein": protein_sequence
    }
    # Write the ORF to the file
    protein_handle.write('>lcl|ORF%i\n%s\n' % (orf_id, protein_sequence))
    nucleotide_handle.write('>lcl|ORF%i\n%s\n' % (orf_id, nucleotides[start:end]))

    return orf_status