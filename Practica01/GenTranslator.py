codones_traduccion = { # Alanina
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    # Cisteina
    "UGU":"C", "UGC":"C",
    # Acido aspartico
    "GAU":"D", "GAC":"D",
    # Acido glutamico
    "GAA":"E", "GAG":"E",
    # Fenilalanina
    "UUU":"F", "UUC":"F",
    # Glicina
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    # Histidina
    "CAU":"H", "CAC":"H",
    # Isoleucina
    "AUA":"I", "AUU":"I", "AUC":"I",
    # Lisina
    "AAA":"K", "AAG":"K",
    # Leucina
    "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    # Metionina
    "AUG":"M",
    # Aspargina
    "AAU":"N", "AAC":"N",
    # Prolina
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    # Glutamina
    "CAA":"Q", "CAG":"Q",
    # Arginina
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
    # Serina
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S",
    # Treonina
    "ACU":"U", "ACC":"U", "ACA":"U", "ACG":"U",
    # Valina
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    # Triptofano
    "UGG":"W",
    # Tirosina
    "UAU":"Y", "UAC":"Y",
    # Stop
    "UAA":"_", "UAG":"_", "UGA":"_"}

# Función para leer la primera secuencia de un archivo fasta.
def read_first_sequence_from_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return lines[1].strip()

# Función para obtener la cadena de DNA complementaria.
def get_complementary_dna(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in dna_sequence])

# Función para obtener el RNA mensajero.
def get_mrna(dna_sequence):
    return dna_sequence.replace('T', 'U')

# Función para traducir RNA a aminoácidos.
def translate_rna_to_aminoacids(rna_sequence, codon_translation):
    amino_acids = []
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_acids.append(codon_translation.get(codon, '-'))
    return ''.join(amino_acids)

# Leer la secuencia del archivo
seq = read_first_sequence_from_fasta("./VentilasHidrotermales/fragment_2.fna")

# Generar cDNA, mRNA y aminoácidos
cDNA = get_complementary_dna(seq)
mRNA = get_mrna(seq)
amino_acids = translate_rna_to_aminoacids(mRNA, codones_traduccion)

# Guardar en archivos fasta
with open("cDNA.fasta", "w") as output:
    output.write(">cDNA_sequence\n" + cDNA)

with open("mRNA.fasta", "w") as output:
    output.write(">mRNA_sequence\n" + mRNA)

with open("aminoacidos.fasta", "w") as output:
    output.write(">amino_acids_sequence\n" + amino_acids)
