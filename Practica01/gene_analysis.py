from Bio import SeqIO
from Bio.Seq import Seq
"""
    Casas Espinosa, Axel
    Jimenez Reyes, Abraham}; 318230577
    Villarreal Maldonado, Jorge Manuel; 307312637
"""
# Leer un archivo FASTA y devolver un diccionario con IDs de secuencia y sus respectivas secuencias.
def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences
# Identificar ORFs en una secuencia de DNA que comiencen con 'ATG' y terminen con 'TAA', 'TAG' o 'TGA'.
def find_orfs(sequence, min_length=100):
    orfs = []
    sequence_length = len(sequence)
    for i in range(sequence_length):
        if sequence[i:i + 3] == 'ATG':
            for j in range(i, sequence_length, 3):
                if sequence[j:j + 3] in {'TAA', 'TAG', 'TGA'}:
                    if j + 3 - i >= min_length:
                        orfs.append(sequence[i:j + 3])
                    break
    return orfs
# Traducir un ORF a una secuencia de aminoácidos.
def translate_orf(orf):
    return Seq(orf).translate(to_stop=True)

file_path = "./VentilasHidrotermales/fragment_1.fna"
sequences = read_fasta(file_path)

for seq_id, seq in sequences.items():
    print(f"ID de la Secuencia: {seq_id}, Longitud: {len(seq)}")

orfs_in_sequences = {}
for seq_id, seq in sequences.items():
    orfs_in_sequences[seq_id] = find_orfs(seq)

for seq_id, orfs in orfs_in_sequences.items():
    print(f"ID de la Secuencia: {seq_id}, Número de ORFs: {len(orfs)}")
    translated_orfs = [translate_orf(orf) for orf in orfs]
