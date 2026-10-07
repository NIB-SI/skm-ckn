'''
CKN v2.0.1: CKN v2 (2023.06) in the format of the PSS network exports, with MapMan updated.
Format and MapMan only: the same nodes and edges as v2.

    python scripts/ckn_v2.0.1.py <CKN v2 folder> <output folder>

<CKN v2 folder> has AtCKN-v2-2023.06.tsv.gz (edges) and AtCKN-v2-2023.06_node-annot.tsv.gz (nodes);
the output folder gets AtCKN-v2.0.1-2026.10.tsv.gz and AtCKN-v2.0.1-2026.10_node-annot.tsv.gz.
MapMan comes from ETL/gmm_annot/transformed_data/ath_mapman4.tsv (scripts/ath_mapman4.py there).

Changes from v2:
- lists joined with '|', as in PSS (synonyms, mapman and interactionSources already were; tissue was ',');
- nodes: node_ID -> id, full_name -> description, GMM -> mapman: GoMapMan 2 (MapMan4) bins of the gene,
  "<bin code>_<full bin name>", as in PSS (v2: classic MapMan; genes only now: no bins for metabolites);
  node_type is the PSS class (PlantCoding, PlantNonCoding, PlantPseudogene, Metabolite, Complex, Process,
  ForeignAbiotic, and for the biotic nodes, by id, ForeignCoding / ForeignEntity / ForeignNonCoding as in PSS);
  the v2 node_type of genes and RNAs (the TAIR locus type, e.g. protein_coding, mirna) moves to locus_type;
  new display_label: the short name, or the id when there is none;
- edges: source, interaction, target first (extended SIF, as PSS): interaction = the influence from effect
  (act: positive-influence, inh: negative-influence, anything else: unknown-influence; effect is kept);
  isDirected 0/1 -> directed False/True (undirected edges are still listed once);
- text: HTML entities spelled out as in PSS (&beta; -> beta, &prime; -> ', &#8208; -> -, the fi/fl
  ligatures -> fi/fl), no-break spaces -> spaces; other Unicode (e.g. o with acute in names) is kept.
'''
import csv
import gzip
import html
import re
import sys
from collections import Counter
import unicodedata
from pathlib import Path

VERSION_IN, VERSION_OUT = 'v2-2023.06', 'v2.0.1-2026.10'
MAPMAN = Path(__file__).resolve().parent.parent / 'ETL/gmm_annot/transformed_data/ath_mapman4.tsv'

INFLUENCE = {'act': 'positive-influence', 'inh': 'negative-influence'}

NODE_COLUMNS = ['id', 'node_type', 'locus_type', 'species', 'TAIR', 'display_label', 'short_name', 'synonyms', 'description',
                'mapman', 'note', 'tissue']
EDGE_COLUMNS = ['source', 'interaction', 'target', 'directed', 'rank', 'effect', 'type', 'species',
                'isTFregulation', 'interactionSources']

# v2 node_type -> PSS class; the v2 types of genes and RNAs are kept as locus_type
LOCUS_TYPES = {
    'protein_coding': 'PlantCoding', 'transposable_element_gene': 'PlantCoding',
    'pseudogene': 'PlantPseudogene',
    'mirna': 'PlantNonCoding', 'antisense_long_noncoding_rna': 'PlantNonCoding', 'pre_trna': 'PlantNonCoding',
    'other_rna': 'PlantNonCoding', 'small_nuclear_rna': 'PlantNonCoding', 'small_nucleolar_rna': 'PlantNonCoding',
}
OTHER_TYPES = {'metabolite': 'Metabolite', 'complex': 'Complex', 'process': 'Process', 'abiotic': 'ForeignAbiotic'}
# biotic nodes by id, with the classes of these entities in PSS
BIOTIC = {
    **dict.fromkeys(['bacteria_flg22', 'virus_6K2', 'virus_CI', 'virus_CP', 'virus_HC-Pro', 'virus_NIa-Pro',
                     'virus_NIb', 'virus_P1', 'virus_P3', 'virus_VPg'], 'ForeignCoding'),
    **dict.fromkeys(['bacteria', 'virus_PVY'], 'ForeignEntity'),
    **dict.fromkeys(['virus_dsRNA', 'virus_vsiRNA', 'virus_me-vsiRNA'], 'ForeignNonCoding'),
}


def node_class(node):
    ''' (PSS class, locus type) of a v2 node; an unknown type or biotic node is an error '''
    node_type = node['node_type']
    if node_type in LOCUS_TYPES:
        return LOCUS_TYPES[node_type], node_type
    if node_type in OTHER_TYPES:
        return OTHER_TYPES[node_type], ''
    if node_type == 'biotic' and node['node_ID'] in BIOTIC:
        return BIOTIC[node['node_ID']], ''
    raise ValueError(f"no PSS class for {node['node_ID']} (node_type {node_type})")


GREEK_LETTERS = {unicodedata.name(chr(c)).split()[-1].lower() for c in range(0x3B1, 0x3CA)}
CHARACTERS = {'‐': '-', '‑': '-', '′': "'", '⁄': '/', ' ': ' '}
ENTITY = re.compile(r'&(?:[A-Za-z]+|#[0-9]+|#x[0-9A-Fa-f]+);')


def spelled_out(character):
    ''' Greek letters spelled out (β -> beta, Δ -> Delta), ligatures split (ﬁ -> fi), look-alikes as ASCII '''
    words = unicodedata.name(character, '').split()
    if words[:2] == ['GREEK', 'SMALL'] and words[-1].lower() in GREEK_LETTERS:
        return words[-1].lower()
    if words[:2] == ['GREEK', 'CAPITAL'] and words[-1].lower() in GREEK_LETTERS:
        return words[-1].capitalize()
    if 'LIGATURE' in words:
        return unicodedata.normalize('NFKC', character)
    return CHARACTERS.get(character, character)


def text(value):
    ''' A value with its HTML entities (and no-break spaces) spelled out as in PSS '''
    value = ENTITY.sub(lambda m: ''.join(spelled_out(c) for c in html.unescape(m.group(0))), value)
    return value.replace(' ', ' ')


def read(path):
    ''' v2 is CSV-quoted where a value has '"' (3 nodes, e.g. "U2B""-LIKE|U2BL") '''
    with gzip.open(path, 'rt', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')


def write(path, columns, rows):
    with gzip.open(path, 'wt', newline='') as f:
        # no quoting, as the PSS exports ('"' is written as it is; values have no tabs or newlines)
        out = csv.DictWriter(f, columns, delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_NONE, quotechar=None)
        out.writeheader()
        out.writerows(rows)
    print(f"{path}: {len(rows)} rows")


def main(folder_in, folder_out):
    folder_in, folder_out = Path(folder_in), Path(folder_out)
    folder_out.mkdir(parents=True, exist_ok=True)

    with open(MAPMAN, newline='') as f:
        mapman = {r['gene']: r['mapman'] for r in csv.DictReader(f, delimiter='\t')}

    nodes = []
    for n in read(folder_in / f'AtCKN-{VERSION_IN}_node-annot.tsv.gz'):
        short_name = text(n['short_name'])
        node_type, locus_type = node_class(n)
        nodes.append({
            'id': n['node_ID'], 'node_type': node_type, 'locus_type': locus_type, 'species': n['species'], 'TAIR': n['TAIR'],
            'display_label': short_name or n['node_ID'], 'short_name': short_name,
            'synonyms': text(n['synonyms']), 'description': text(n['full_name']),
            'mapman': mapman.get(n['node_ID'].upper(), ''), 'note': text(n['note']),
            'tissue': '|'.join(t.strip() for t in n['tissue'].split(',') if t.strip()),
        })
    write(folder_out / f'AtCKN-{VERSION_OUT}_node-annot.tsv.gz', NODE_COLUMNS, nodes)
    for (node_type, locus_type), count in sorted(Counter((n['node_type'], n['locus_type']) for n in nodes).items()):
        print(f'  {node_type:17} {locus_type or "-":29} {count}')

    edges = []
    for e in read(folder_in / f'AtCKN-{VERSION_IN}.tsv.gz'):
        edges.append({
            'source': e['source'], 'interaction': INFLUENCE.get(e['effect'], 'unknown-influence'),
            'target': e['target'], 'directed': {'1': True, '0': False}[e['isDirected']], 'rank': e['rank'],
            'effect': e['effect'], 'type': e['type'], 'species': e['species'],
            'isTFregulation': e['isTFregulation'], 'interactionSources': e['interactionSources'],
        })
    write(folder_out / f'AtCKN-{VERSION_OUT}.tsv.gz', EDGE_COLUMNS, edges)


if __name__ == '__main__':
    main(*sys.argv[1:3])
