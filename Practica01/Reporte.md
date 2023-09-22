# 🧬 El Dogma Central de la Biología Molecular 🧪
# 🛠 Práctica 01 📑

Introducción
En este proyecto, estamos explorando la genómica de organismos que habitan en ambientes extremos, como las ventilas 
hidrotermales. Utilizando muestras recolectadas de estas áreas, hemos secuenciado fragmentos de ADN que podrían 
pertenecer a organismos desconocidos. El objetivo de este proyecto es analizar estas secuencias de ADN para identificar
posibles genes y proteínas codificadas por estos fragmentos de ADN.

## Dependencias
* Python: Utilizamos Python para escribir y ejecutar nuestros scripts de análisis de datos. Asegúrate de tener una versión
reciente de Python instalada.
* Biopython: Una herramienta esencial para el análisis bioinformático en Python. Se utiliza para leer archivos FASTA y 
traducir secuencias de nucleótidos a aminoácidos.

## Parte 1: Descripción del Proyecto

Eres el computólogo a cargo en un equipo multidisciplinario donde tus colegas son biólogos, matemáticos y terrólogos.
Resulta que tu grupo se dedica a investigar la vida que hay en ambientes extremos. Ésta vez,lograron extraer muestras 
provenientes de ventilas hidrotermales. Los biólogos se hicieron cargo de su trabajo y te otorgaron ésta carpeta con 
cuatro archivos FASTA.

## Marco teorico

Un gen es una unidad de herencia que ocupa una ubicación específica (locus) en un cromosoma. En el contexto de la 
biología molecular, un gen es una secuencia de nucleótidos en el ADN que codifica la síntesis de una cadena de 
polipéptidos o de una molécula de ácido ribonucleico (ARN) con una función conocida.

Un Marco de Lectura Abierto (ORF, por sus siglas en inglés, Open Reading Frame) es una secuencia continua de nucleótidos
que tiene la potencialidad de codificar una proteína. Un ORF comienza con un codón de inicio (generalmente ATG, que 
codifica para el aminoácido metionina) y termina con uno de los tres codones de parada (TAA, TAG o TGA), sin ningún otro
codón de parada en medio. Un ORF, por lo tanto, representa una parte de la secuencia de un gen que tiene el potencial de
codificar una proteína. Para lograr este trabajo, se siguieron los siguientes pasos:

## Análisis de Búsqueda de Genes

Para realizar el análisis de búsqueda de genes en los archivos de secuencia, utilizamos una combinación de scripts de 
Bash y comandos de línea de comandos, incluyendo `grep` y `sed`. Para contar el número de genes (o encabezados de secuencia)
presentes en cada archivo, utilizamos el comando `grep`con la opción `-c` para contar las líneas que comienzan con el 
carácter `>` (que denota el inicio de una nueva secuencia). 



```bash
grep -c "^>" *.fasta
  fragment_1.fna:1
  fragment_2.fna:1
  fragment_3.fna:1
  fragment_4.fna:1
```
### Cálculo de la longitud de las secuencias

Para calcular la longitud de las secuencias en un archivo FASTA, podemos utilizar el comando awk. 

```bash
awk '{if($0 !~ />/) print length($0)}' file_name.fna
```


### Análisis de Codones de Inicio y Parada 

Para llevar a cabo un análisis inicial y localizar los codones de inicio y parada en las secuencias, podemos utilizar el
comando grep. A continuación, presentamos una guía sobre cómo usar este comando para identificar líneas que contienen los
codones de inicio ("ATG") y los codones de parada ("TAA", "TAG", "TGA") en un archivo FASTA:

#### Identificación de Codones de Inicio y Parada

* Codón de Inicio - **"ATG"**  
Para identificar todas las líneas que contienen el codón de inicio "ATG", utiliza el siguiente comando:

```bash
grep -n 'ATG' file_name.fna
```
* Codón de Parada - **"TAA"**  
Para identificar todas las líneas que contienen el codón de parada "TAA", utiliza el siguiente comando:

```bash
grep -n 'TAA' file_name.fna
```
* Codón de Parada - "TAG"  
Para identificar todas las líneas que contienen el codón de parada "TAG", utiliza el siguiente comando:
```bash
grep -n 'TAG' file_name.fna
```
* Codón de Parada - "TGA"  
Para identificar todas las líneas que contienen el codón de parada "TGA", utiliza el siguiente comando:
```bash
grep -n 'TGA' file_name.fna
```

#### Identificación de Posibles Genes (ORFs) 

Para identificar posibles genes (ORFs) en un archivo FASTA, podemos utilizar grep para buscar marcos de lectura abiertos
(ORFs) que comienzan con el codón de inicio "ATG" y terminan con uno de los codones de parada ("TAA", "TAG" o "TGA"):

```bash
grep -o -P 'ATG(?:[ATGC]{3}){30,}?T(?:AA|AG|GA)' file_name.fna | wc -l
```
##### Explicación del Comando

`grep`: Una herramienta para buscar cadenas específicas en archivos.
* `-o`: Opción que permite imprimir solo las partes de la línea que coinciden con el patrón.
* `-P`: Opción que permite interpretar la expresión como una expresión regular de Perl, lo que habilita el uso del constructo de no captura (?: ...).
* `'ATG(?:[ATGC]{3}){30,}?T(?:AA|AG|GA)'`: Es una expresión regular que busca una cadena que comienza con "ATG", seguida de cualquier combinación de 30 o más tripletes de nucleótidos (asegurando así una longitud mínima para el ORF), y que termina con uno de los codones de parada ("TAA", "TAG" o "TGA").
* `file_name.fna`: Es el archivo donde se está realizando la búsqueda.
* `|`: Es un operador de tubería que pasa la salida del comando anterior como entrada para el siguiente comando.
* `wc -l`: Cuenta el número de líneas en la entrada, lo que en este caso indica el número de coincidencias encontradas por grep.

## Utilizando BioPython para la Identificación de ORFs

En este proyecto, también nos apoyaremos en la biblioteca BioPython, una herramienta poderosa y flexible para el análisis 
computacional de secuencias biológicas. La biblioteca BioPython facilita la escritura de scripts de Python para trabajar 
con datos biológicos.

### Instalación de BioPython

Antes de comenzar, asegúrate de tener instalada la biblioteca BioPython. Puedes instalarla utilizando el siguiente comando:

```bash
pip install biopython
```

### Lectura de Archivos FASTA
Utilizaremos BioPython para leer los archivos FASTA y extraer las secuencias que contienen. Esto se puede hacer usando 
el módulo SeqIO de BioPython.

### Identificación de ORFs
Posteriormente, escribiremos scripts para identificar ORFs en las secuencias extraídas. Esto implica buscar secuencias 
que comiencen con el codón de inicio "ATG" y terminen con un codón de parada ("TAA", "TAG" o "TGA"), con una longitud 
mínima especificada para considerarse un ORF válido.

### Codigo en Python

Este script de Python utiliza la biblioteca BioPython para analizar archivos en formato FASTA. Identifica y traduce los
marcos de lectura abiertos (ORFs) en las secuencias nucleotídicas presentes en el archivo. Los ORFs se definen como 
secuencias que comienzan con el codón "ATG" y terminan con uno de los codones de parada ("TAA", "TAG" o "TGA"), con una 
longitud mínima de 90 nucleótidos (30 tripletes). El script también ofrece detalles sobre la longitud de cada secuencia 
y el número de ORFs identificados.

```python
from Bio import SeqIO
from Bio.Seq import Seq

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

file_path = "./fasta/fragment_1.fna"
sequences = read_fasta(file_path)

for seq_id, seq in sequences.items():
    print(f"ID de la Secuencia: {seq_id}, Longitud: {len(seq)}")

orfs_in_sequences = {}
for seq_id, seq in sequences.items():
    orfs_in_sequences[seq_id] = find_orfs(seq)

for seq_id, orfs in orfs_in_sequences.items():
    print(f"ID de la Secuencia: {seq_id}, Número de ORFs: {len(orfs)}")
    translated_orfs = [translate_orf(orf) for orf in orfs]

```

## Resultados (Parte 1)

### Análisis de Fragmentos de Secuencia

#### Archivo Fragmento 1:

- Número de registros (posibles genes): 1
- Detalles del registro:
  - ID: Fragmento1
  - Longitud de la secuencia: 632,428 bases
  - Número de ORFs identificados: 3,611

#### Archivo Fragmento 2:

- Número de registros (posibles genes): 1
- Detalles del registro:
  - ID: Fragmento2
  - Longitud de la secuencia: 100,531 bases
  - Número de ORFs identificados: 527

#### Archivo Fragmento 3:

- Número de registros (posibles genes): 1
- Detalles del registro:
  - ID: Fragmento3
  - Longitud de la secuencia: 1,154,456 bases
  - Número de ORFs identificados: 5,741

#### Archivo Fragmento 3:

- Número de registros (posibles genes): 1
- Detalles del registro:
  - ID: Fragmento4
  - Longitud de la secuencia:7,115,445 bases
  - Número de ORFs identificados: 35,429

### Observaciones:

* **Fragmento 1:** Aunque su tamaño sugiere que podría ser un genoma bacteriano, el número de ORFs identificados 
(3,611) es típico para muchas bacterias. Esto refuerza la idea de que este fragmento puede ser un genoma bacteriano completo.

* **Fragmento 2:** A pesar de su menor tamaño, el número de ORFs (527) es notable. Aunque no es suficiente para ser un 
genoma bacteriano completo, podría representar un plásmido o un genoma viral.

* **Fragmento 3:** Con 5,741 ORFs, este fragmento tiene una diversidad genética considerable. Esto, junto con su longitud, 
sugiere que es muy probable que sea un genoma bacteriano completo.

* **Fragmento 4:** La enorme cantidad de ORFs (35,429) refuerza la idea de que este fragmento pertenece a un organismo
eucariota, como un hongo o incluso un protista.

## Parte 2: Identificación y Extracción de un Gen Clave en Ambientes Extremos

Tus colegas descubren que el primer gen del fragmento de DNA más corto, podría ser importante para los organismos que 
viven en ventilas hidrotermales. 

### Código Python (GenTranslator)

El programa `Gen Translator` se diseñó para analizar y procesar genes de interés en secuencias de DNA. A partir de un 
archivo `.fasta` con una secuencia genética, el programa realiza las siguientes tareas:

* **Generación del cDNA**: Calcula la secuencia complementaria de DNA a partir de la secuencia de entrada y la guarda en 
un archivo `cDNA.fasta`.
* **Transcripción a ARNm**: Convierte la secuencia de DNA en su correspondiente ARN mensajero, eliminando timinas y 
reemplazándolas por uracilos. El resultado se almacena en `mRNA.fasta`.
* **Traducción a Aminoácidos**: Traduce el ARN mensajero en una cadena de aminoácidos utilizando el código genético,
generando una secuencia proteica. Esta secuencia se guarda en `aminoacidos.fasta`.

```python
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
seq = read_first_sequence_from_fasta("file_name.fna")

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
```

### Resultado (Parte2)

A partir de la información genómica y proteómica obtenida y con base en los archivos generados (`cDNA.fasta`, `mRNA.fasta`,
`aminoacidos.fasta`), hemos deducido que el organismo en estudio es eucarionte. A continuación, se detallan las razones principales:

#### Gen Relacionado con la Formación del Núcleo Celular

Se identificó un gen que está asociado con la formación o función del núcleo celular. Es esencial destacar que solo las
células eucariontes poseen un núcleo celular bien definido. La presencia de este gen es una fuerte evidencia de la 
naturaleza eucariota del organismo.

#### cDNA y su Relevancia en Eucariontes

El archivo `cDNA.fasta` representa el ADN complementario formado a partir del ARNm. En eucariontes, el cDNA es crucial 
para estudiar la expresión génica, ya que refleja exclusivamente los genes expresados, excluyendo las regiones intrónicas. 
Los procariotas carecen de intrones, por lo que la importancia del cDNA en este contexto sugiere un origen eucarionte.



#### cDNA y su Relevancia en Eucariontes
El archivo `cDNA.fasta` representa el ADN complementario formado a partir del ARNm. En eucariontes, el cDNA es crucial 
para estudiar la expresión génica, ya que refleja exclusivamente los genes expresados, excluyendo las regiones intrónicas.
Los procariotas carecen de intrones, por lo que la importancia del cDNA en este contexto sugiere un origen eucarionte.

###### Genes en el ADN: Exones e Intrones

Los genes en el ADN están compuestos por dos tipos principales de secuencias: exones e intrones.

###### Exones

Son las secuencias de ADN que se transcriben y traducen en proteínas. Es decir, tienen la información codificada que se
utilizará para producir una proteína específica.

###### Intrones

Son las secuencias de ADN que se encuentran entre los exones, pero no se traducen en proteínas. Durante el proceso de 
formación del ARN mensajero (ARNm) en eucariotas, los intrones se transcriben inicialmente, pero luego son eliminados en
un proceso llamado "empalme" o "splicing", dejando solo los exones en el ARNm maduro.

#### Características del ARN Mensajero

El archivo `mRNA.fasta` contiene secuencias de ARNm. En eucariontes, este ARNm pasa por un proceso de maduración que 
incluye adiciones específicas y el empalme para eliminar intrones. Estas características, si se detectan en el archivo, 
indican un proceso típico de maduración del ARNm eucarionte.

####  Análisis Proteómico

El archivo `aminoacidos.fasta` presenta las proteínas traducidas a partir del ARNm. Un análisis posterior de estas 
secuencias reveló funciones específicas asociadas a eucariontes.
