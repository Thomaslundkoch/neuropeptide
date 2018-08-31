# Additional file X: Program for identification of preprohormones in the transcriptomes from various cubomedusae
# The script ‘neuropeptideFinder.py’ is written in Python3.
# To apply it to the transcriptome, the relevant database should be connected to the script.



def reverse_complement(sequence):
    """input: DNA string. output: reverse complement of DNA string
    Revised by TLK 18.07.18"""
    reverse_sequence=sequence[::-1]
    reverse_complement=''
    substitutions={'A':'T','T':'A','G':'C','C':'G','N':'N'}
    for nucleotide in reverse_sequence:
        try:
            reverse_complement+=substitutions[nucleotide]
        except:
            continue
    return reverse_complement

def translate(seq):
    """Input nucleotide string. Output: translated protein"""

    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""

    for i in range(0, len(seq), 3):
        try:
            codon = seq[i:i + 3]
            protein+= table[codon]
        except:
            quit
    return protein

def find_neuropeptides_transcriptome(database, motifs, processing_site_threshold, peptide_length,identity_threshold):
    """Input: database is a nucleotide sequence containing fasta-file from transcriptome, the other inputs are to the helper functions
    Output: list of proteins that passed both criteria (motifs and peptide identity)
    The program translate the sequence in all three positive reading frames and splits them into segments without stop-codons, which are analyzed
    Written by TLK, 17.06.18"""
    title=''
    sequence=''
    for line in database:
        line=line.replace('\n','')
        if line.startswith('>'):
            reverse_complement_sequence=reverse_complement(sequence)
            for reading_frame in [translate(sequence),translate(sequence[1:]),translate(sequence[2:]),
            translate(reverse_complement_sequence),translate(reverse_complement_sequence[1:]),translate(reverse_complement_sequence[2:])]:
                for open_reading_frame in reading_frame.split('*'):
                    if len(open_reading_frame)>1500: #ignore proteins that are too long to encode prepropeptide
                        title=line
                        sequence=''
                    elif 'GRGRGRGRGR' in open_reading_frame: #ignore repetitive motifs of GR's that are somewhat common
                        title=line
                        sequence=''
                    else:
                        (search,peptide_list)=count_and_list(open_reading_frame,motifs,processing_site_threshold,peptide_length)
                        if search:
                            if peptide_identity(peptide_list,identity_threshold):
                                print(title)
                                print(open_reading_frame)
            title=line
            sequence=''
        else:
            sequence+=line

find_neuropeptides_transcriptome(database, ['GR','GKR','GKK'],7,5,80)
