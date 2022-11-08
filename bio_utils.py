from Bio.Seq import Seq

def transcribe_dna(sequence):
    # Nucleotide sequences are normally read from the 5' to 3; direction
    coding_dna = sequence
    template_dna = coding_dna.reverse_complement()

    # Want to transcribe the coding strand into the corresponding mRNA
    messenger_rna = coding_dna.transcribe()
    
    # If you do want to do a true biological transcription starting with 
    # the template strand, then this becomes a two-step process:
    template_dna.reverse_complement().transcribe() 
    
    return messenger_rna

# Accepts mRNA sequence of coding DNA
def translate_sequence_to_protein(seq, table = "Standard", stop_symbol = "*", cds = True):
    return seq.translate(table=table, stop_symbol=stop_symbol, cds=cds)

# Most protein sequences are derived from translations of CoDing Sequence (CDS)
# derived from gene predictions. A CoDing Sequence (CDS) is a region of DNA or RNA 
# whose sequence determines the sequence of amino acids in a protein. It should not 
# be mixed up with an Open Reading Frame (ORF), which is a continuous stretch of DNA 
# codons that begins with a start codon and ends at a STOP codon. All CDS are ORFs, 
# but not all ORFs are CDS...
def get_cds_feature_from_sequence_record(seq_record):
    for feature in seq_record.features:
        if feature.type == 'CDS':
            return feature
    
    return None

# Check if the translated open reading frame contains the coding sequence protein
# described in the genbank file. 
def check_if_valid_orf_with_cds(seq_record, orf_protein):
    cds_feature = get_cds_feature_from_sequence_record(seq_record)
        
    gene_sequence = cds_feature.extract(seq_record.seq)
    protein_sequence = translate_sequence_to_protein(seq=gene_sequence, cds=True)

    # Halt if the translation does not match
    assert protein_sequence == cds_feature.qualifiers["translation"][0]

    # ORF must include the protein sequence of the CDS for it to be a valid ORF
    return orf_protein.count(str(protein_sequence)) == 1

def get_cds_protein_from_record(seq_record):
    cds_feature = get_cds_feature_from_sequence_record(seq_record)    
    return cds_feature.qualifiers["translation"][0]
