'''
Arabidopsis genes -> MapMan4 bins (GoMapMan 2), for CKN and the PSS functional clusters.

Input (raw_data/, downloaded from the GoMapMan 2 exports page, https://gomapman2.nib.si):
    MapMan4_gene-binCode_pairs.zip    Identifier, plantName, binCode
    MapMan4_ontology.zip              binCode, binName, ...

Output: transformed_data/ath_mapman4.tsv, one row per gene:
    gene      AGI, upper case (AT2G38470)
    mapman    the gene's bins, '|'-joined (as in CKN and the PSS downloads), each "<bin code>_<full bin name>", the full name being the names of
              the bin and its parents joined with '.' (bin names contain no '.', ';' or '|'), e.g.
              26.11.3.2.1_External stimuli response.pathogen.defense mechanisms.WRKY33-dependent plant
              immunity.transcription factor *(WRKY33)
Bin 35 (not assigned / not annotated) is left out: no bins = no annotation.

    python scripts/ath_mapman4.py
'''
import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
RAW, OUT = BASE / 'raw_data', BASE / 'transformed_data'
NOT_ASSIGNED = '35'


def read_zipped_tsv(path):
    with zipfile.ZipFile(path) as z:
        (name,) = z.namelist()
        yield from csv.DictReader(io.TextIOWrapper(z.open(name), 'utf-8'), delimiter='\t')


def main():
    names = {r['binCode']: r['binName'] for r in read_zipped_tsv(RAW / 'MapMan4_ontology.zip')}

    def label(code):
        parts = code.split('.')
        return code + '_' + '.'.join(names['.'.join(parts[:i])] for i in range(1, len(parts) + 1))

    def sort_key(code):
        return [int(x) for x in code.split('.')]

    bins = defaultdict(set)
    for r in read_zipped_tsv(RAW / 'MapMan4_gene-binCode_pairs.zip'):
        code = r['binCode']
        if r['plantName'] == 'Arabidopsis' and code.split('.')[0] != NOT_ASSIGNED:
            bins[r['Identifier'].upper()].add(code)

    OUT.mkdir(exist_ok=True)
    with open(OUT / 'ath_mapman4.tsv', 'w', newline='') as f:
        out = csv.writer(f, delimiter='\t', lineterminator='\n')
        out.writerow(['gene', 'mapman'])
        for gene in sorted(bins):
            out.writerow([gene, '|'.join(label(c) for c in sorted(bins[gene], key=sort_key))])
    print(f"{len(bins)} genes, {sum(map(len, bins.values()))} gene-bin pairs -> {OUT / 'ath_mapman4.tsv'}")


if __name__ == '__main__':
    main()
