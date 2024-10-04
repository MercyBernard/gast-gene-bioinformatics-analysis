from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from collections import Counter

# Provide email for NCBI's Entrez tool
Entrez.email = "mercybernard.3399@gmail.com"

# Fetch the gene sequence by accession number
accession_number = "NM_000805"  # Accession for GAST gene
handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")

# Parse and read the sequence
record = SeqIO.read(handle, "genbank")
handle.close()

print(f"ID: {record.id}")
print(f"Description: {record.description}")

for feature in record.features:
    if feature.type == "CDS":
        cds_sequence = feature.extract(record.seq)
        print(f"CDS Sequence: {cds_sequence}")
        print(f"CDS Location: {feature.location}")
        print(f"CDS Qualifiers: {feature.qualifiers}")
        break

gast_sequence = str(cds_sequence)

print("GAST coding sequence lenght:", len(gast_sequence))

nucleotide_count = Counter(gast_sequence)
print(nucleotide_count)


# Calculate GC content
def calculate_gc_content(dna_sequence):

    dna_sequence = gast_sequence

    g_count = gast_sequence.count('G')
    c_count = gast_sequence.count('C')

    gc_content = ((g_count + c_count)/len(dna_sequence)) * 100

    return gc_content

gc_percentage = calculate_gc_content(gast_sequence)

print(f"GC Content: {gc_percentage:.2f}%")


#Get the Reverse Complementary strand for the DNA sequence
def reverse_complement_strand(dna_sequence):

    dna_sequence = gast_sequence
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_comp = ''.join([complement[base] for base in dna_sequence[::-1]])

    return rev_comp

reverse_complement_sequence = reverse_complement_strand(gast_sequence)
print(f"Reverse Complement Strand: {reverse_complement_sequence}")


def transcribe_dna_to_rna(dna_sequence):

    dna_sequence = gast_sequence
    rna_sequence = dna_sequence.replace('T', 'U')
    return rna_sequence

rna_sequence = transcribe_dna_to_rna(gast_sequence)
print(f"RNA Sequence: {rna_sequence}")


def translate_rna_to_protein(rna_sequence):

    codon_dictinary = {
        'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
        'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
        'UAU':'Y', 'UAC':'Y', 'UAA':'Stop', 'UAG':'Stop',
        'UGU':'C', 'UGC':'C', 'UGA':'Stop', 'UGG':'W',
        'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
        'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
        'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
        'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
        'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
        'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
        'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    }

    rna_sequence = rna_sequence
    protein_sequence = []

    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]

        amino_acid = codon_dictinary.get(codon, '*')
        if amino_acid == 'Stop':
            break

        protein_sequence.append(amino_acid)

    return ''.join(protein_sequence)

protein_sequence = translate_rna_to_protein(rna_sequence)
print(f"Protein Sequence: {protein_sequence}")


def mutate_sequence(dna_sequence, position, new_base):

    dna_sequence = gast_sequence

    dna_sequence_list = list(dna_sequence)
    dna_sequence_list[position] = new_base
    return ''.join(dna_sequence_list)

mutate_dna_sequence = mutate_sequence(gast_sequence, 9, 'T')
print(f'Mutated Sequence: {mutate_dna_sequence}')


def simulate_primer_design(dna_sequence, reverse_complement_sequence, primer_length=20):

    dna_sequence = gast_sequence
    reverse_complement_sequence = reverse_complement_sequence

    forward_primer = dna_sequence[:primer_length]
    reverse_primer = reverse_complement_sequence[-primer_length:]
    return forward_primer, reverse_primer

forward_primer, reverse_primer = simulate_primer_design(gast_sequence, reverse_complement_sequence)
print(f"Forward Primer: {forward_primer}")
print(f"Reverse Primer: {reverse_primer}")