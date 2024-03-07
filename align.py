# Il programma prende in input un file bam e un file gtf e conta quanti allineamenti ci sono per ogni trascritto nel file gtf prendendo come riferimento le reads del file bam.
# Successivamente stampa a schermo una tabella con una riga per ogni trascritto e una colonna per gli allineamenti.
# E' possibile specificare un parametro opzionale -q, --quality che specifica il valore di qualità minimo per gli allineamenti da considerare. Il valore di default è 0.
# Esempio di utilizzo (python): python align.py -bam file.bam -gtf file.gtf -q 30
# Esempio di utilizzo (python3): python3 align.py -bam file.bam -gtf file.gtf -q 30

import pysam
import argparse
import pandas as pd
import re

from pysam import AlignmentFile

# Tolleranza per l'allineamento degli esoni
TOLLERANCE = 1

# Funzione per controllare se un blocco è allineato su un esone
def check_block(exons, block):
    for exon in exons:
        if block[0] >= exon[0] - TOLLERANCE and block[1] <= exon[1] + TOLLERANCE:
            return True
    return False

# Funzione per controllare se un blocco è allineato su un esone e il blocco successivo è allineato sull'esone successivo (ovvero se gli introni sono allineati)
def check_next_block(exons, block, block_next):
    exons_iter = iter(exons)
    for exon in exons_iter:
        next_exon = next(exons_iter, -1)
        if next_exon != -1:
            if exon[1] - TOLLERANCE <= block[1] <= exon[1] + TOLLERANCE and next_exon[0] - TOLLERANCE <= block_next[0] <= next_exon[0] + TOLLERANCE:
                return True
    return False

# Funzione per controllare se una read è allineata sugli esoni
def check_read(exons, read):
    block_iter = iter(read.get_blocks())
    for block in block_iter:
        if not check_block(exons, block):
            return False
        next_item = next(block_iter, -1)
        if next_item != -1:
            if not check_next_block(exons, block, next_item):
                return False
    return True

def main(args):
    # Dizionario vuoto per i trascritti e i loro allineamenti
    transcripts_dict = {}

    # Indicizzazione del file bam
    pysam.index(args.bam_file)
    bam_file = AlignmentFile(args.bam_file, 'rb')
    # Lettura del file gtf
    with open(args.gtf_file, 'r') as gtf_file:
        for line in gtf_file:
            fields = line.strip().split('\t')
            if fields[0][0] == '#':
                continue
            # Estrazione delle informazioni sui trascritti
            if fields[2] == 'transcript':
                id = re.search(r'transcript_id "([^"]+)"', line).group(1)
                # Estrazione del cromosoma di riferimento
                chrom = fields[0]
                # Estrazione delle coordinate del trascritto
                start = int(fields[3])
                end = int(fields[4])
                # Aggiunta del trascritto al dizionario
                transcripts_dict[id] = {'chrom': chrom, 'start': start, 'end': end, 'alignments': 0, 'exons': list()}
            if fields[2] == 'exon':
                id = re.search(r'transcript_id "([^"]+)"', line).group(1)
                # Estrazione delle coordinate dell'esone
                start = int(fields[3])
                end = int(fields[4])
                # Aggiunta dell'esone al dizionario
                transcripts_dict[id]['exons'].append((start, end))
    gtf_file.close()
    # Scorrimento del file bam per contare gli allineamenti
    for item in transcripts_dict:
        chrom = transcripts_dict[item]['chrom']
        start = transcripts_dict[item]['start']
        end = transcripts_dict[item]['end']
        exons = transcripts_dict[item]['exons']
        for read in bam_file.fetch(chrom, start, end):
            # Se la read è allineata sugli esoni e ha una qualità maggiore o uguale a quella specificata, la aggiungo al set
            if check_read(exons, read) and read.mapping_quality >= args.quality:
                transcripts_dict[item]['alignments'] += 1

    # Chiusura del file
    bam_file.close()
    # Creazione del dataframe e stampa a schermo
    df = pd.DataFrame.from_dict(transcripts_dict, orient='index', columns=['alignments'])
    df.index.name = 'ID del trascritto'
    df.columns = ['Allineamenti']
    print(df)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Conta quanti allineamenti completi e parziali ci sono per ogni trascritto nel file gtf prendendo come riferimento il file bam.')

    # Aggiunta degli argomenti: file bam, file gtf e qualità minima
    parser.add_argument('-bam', '--bam_file', type=str, help='File bam in input', required=True)
    parser.add_argument('-gtf', '--gtf_file', type=str, help='File gtf in input', required=True)
    parser.add_argument('-q', '--quality', type=int, help='Valore di qualità minimo', required=False, default=0)
    
    # Parsing degli argomenti
    args = parser.parse_args()
    main(args)