'''
CKN v2.1: CKN v2 (2023.06) in the format of the PSS network exports, with MapMan updated and its PSS part
brought up to date with PSS.

    python scripts/ckn_v2.1.py <CKN v2 folder> <PSS gene network edges> <PSS gene network nodes> \
        <deprecated PSS reactions> <output folder>

<CKN v2 folder> has AtCKN-v2-2023.06.tsv.gz (edges) and AtCKN-v2-2023.06_node-annot.tsv.gz (nodes);
the PSS gene network is pss-export's 'gene-network' for public Arabidopsis reactions (edges.tsv, nodes.tsv),
the deprecated reactions a list of ids, one per line (the DeprecatedReaction nodes of the same database); both
in ETL/pss/raw_data/. The output folder gets AtCKN-v2.1-2026.10.tsv.gz and AtCKN-v2.1-2026.10_node-annot.tsv.gz.
MapMan comes from ETL/gmm_annot/transformed_data/ath_mapman4.tsv (scripts/ath_mapman4.py there).

The PSS part (edges with skm:rx sources, rank 0; see update_pss):
- nodes have their PSS names (v2 renamed some, e.g. virus_VPg, Calcium: PSS_NAMES), and non-gene nodes their
  PSS annotation;
- reactions no longer in the PSS gene network (deprecated, or now curated for other species only) are removed,
  as are those whose v2 edges PSS no longer has (e.g. complexes it no longer makes) or has with the opposite
  sign (opposite_signs): their skm:rx source is dropped, and an edge left without sources goes;
- the reactions the CKN does not have (new, and those just removed for changed edges) are added from the gene
  network (PSS_TYPES), merged into an existing edge between the same nodes (merge_pss_edge);
- the other v2 PSS edges are kept as they are (the 2023 conversion);
- an edge left with other sources only gets their rank (source_rank);
- nodes left without edges are removed.

Changes from v2 (format):
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

VERSION_IN, VERSION_OUT = 'v2-2023.06', 'v2.1-2026.10'
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


def ckn_v2(folder_in, mapman):
    ''' The v2 nodes and edges in the new format '''
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
    edges = []
    for e in read(folder_in / f'AtCKN-{VERSION_IN}.tsv.gz'):
        edges.append({
            'source': e['source'], 'interaction': INFLUENCE.get(e['effect'], 'unknown-influence'),
            'target': e['target'], 'directed': {'1': True, '0': False}[e['isDirected']], 'rank': e['rank'],
            'effect': e['effect'], 'type': e['type'], 'species': e['species'],
            'isTFregulation': e['isTFregulation'], 'interactionSources': e['interactionSources'],
        })
    return nodes, edges


##############################################################
# The PSS part
##############################################################

# v2 names of PSS nodes -> their PSS names (found through the same reactions in the PSS gene network)
PSS_NAMES = {
    'virus_6K2': '6K2', 'virus_CI': 'CI', 'virus_CP': 'CP', 'virus_HC-Pro': 'HC-Pro', 'virus_NIa-Pro': 'NIa-Pro',
    'virus_NIb': 'NIb', 'virus_P1': 'P1', 'virus_P3': 'P3', 'virus_VPg': 'VPg', 'virus_PVY': 'potyvirus',
    'virus_dsRNA': 'viral dsRNA', 'virus_vsiRNA': 'vsiRNA', 'virus_me-vsiRNA': 'me-vsiRNA',
    'bacteria': 'trichous-bacteria', 'bacteria_flg22': 'flg22',
    'abiotic_drought': 'Drought', 'abiotic_heat': 'Heat', 'abiotic_waterlogging': 'Waterlogging',
    'ath_anthocyanins': 'Anthocyanin-accumulation', 'ath_co2-deficit': 'CO2-deficiency', 'ath_growth': 'Growth',
    'ath_o2-deficit': 'O2-deficiency', 'ath_phenolics': 'Phenolics production', 'ath_stomata': 'Stomatal closure',
    'ath_thermotolerance': 'Thermotolerance', 'ath_trichomes': 'Trichome-initiation',
    'Calcium': 'Ca2+', 'Copper': 'Cu2+', 'electron': 'e-', 'Abscisic-aldehyde': 'Abscisic aldehyde',
    'IAOx-N-oxide': 'IAOx N-oxide', 'apo-beta-carotenal': "9-cis-10'-apo-beta-carotenal",
    'ent-Kaurenoic-acid': 'ent-Kaurenoic acid', 'p-Coumaric-acid': 'p-Coumaric acid',
    'R-gene': 'TNL[fc00004]', 'GPAphid2': 'GPAphid2[fc00002]', 'Rx': 'Rx[fc00196]', 'iaaM': 'iaaM[fc00253]',
    'BAK1|FLS2|bacteria_flg22': 'BAK1|FLS2|flg22', 'BIK1|PEPR|PEP': 'BIK1|PEPR1|PEP1', 'PEPR|PEP': 'PEPR1|PEP1',
    'CDK|TGA2,5,6': 'CDK|TGA', 'NPR1|TGA2,5,6': 'NPR1|TGA', 'SCL14|TGA2,5,6': 'SCL14|TGA',
    'CML|Calcium': 'CML|Ca2+', 'CML|virus_HC-Pro': 'CML|HC-Pro', 'GA|GID': 'GA|GID1', 'DELLA|GA|GID': 'DELLA|GA|GID1',
    'EDS1|MPK3,6|PAD4': 'EDS1|MPK3|PAD4', 'NPR1|WRKY6,18': 'NPR1|WRKY18', 'R-gene|virus_PVY': 'R-gene|potyvirus',
    'RISC|MIR390': 'RISC - miR390', 'RISC|virus_vsiRNA': 'RISC - vsiRNA',
}

# PSS reaction type -> CKN type, the others 'other' (as ETL/pss/scripts/pss2ckn_20260720.ipynb)
PSS_TYPES = {
    'transcriptional/translational activation': 'transcription factor regulation',
    'transcriptional/translational repression': 'transcription factor regulation',
    'binding/oligomerisation': 'binding', 'dissociation': 'binding',
    'protein activation': 'post-translational modification', 'protein deactivation': 'post-translational modification',
    'degradation/secretion': 'post-translational modification',
}
TF_REGULATION = 'transcription factor regulation'
PTM = 'post-translational modification'
REPRESSION = 'transcriptional/translational repression'
EFFECTS = {'positive-influence': 'act', 'negative-influence': 'inh'}
PSS_SOURCE = re.compile(r'^skm:(rx\d+)$')
GENE = re.compile(r'^(AT[1-5CM]G\d{5})(\.\d+)?$', re.IGNORECASE)


def node_id(name):
    ''' A gene without its version (AT5G66900.1 -> AT5G66900) '''
    m = GENE.match(name)
    return m.group(1).upper() if m else name


def sources(edge):
    return [s for s in edge['interactionSources'].split('|') if s]


def reactions(edge):
    return {m.group(1) for s in sources(edge) if (m := PSS_SOURCE.match(s))}


def pairs(edge):
    ''' The node pairs of an edge: an undirected edge is listed once, for both directions '''
    a, b = edge['source'], edge['target']
    return {(a, b)} if edge['directed'] else {(a, b), (b, a)}


def combined_effect(a, b):
    if a == b or b in ('', 'unk'):
        return a
    if a in ('', 'unk'):
        return b
    return 'act/inh'


def set_effect(edge, effect):
    edge['effect'] = effect
    edge['interaction'] = INFLUENCE.get(effect, 'unknown-influence')


def merge_pss_edge(edge, new, log):
    ''' A PSS edge (new) between the nodes of an existing CKN edge: its source is added and the rank set to 0.
    Type and effect as in ETL/pss/scripts/pss2ckn_20260720.ipynb (cells 39, 40):
    - the same type: the effects combined (act and inh: act/inh);
    - CKN binding and a functional PSS type, or CKN rank 4: PSS's type and effect;
    - PTM on both sides: PSS's effect where CKN's is neither act nor inh;
    - PSS binding and a functional CKN type, CKN small RNA interactions and PSS repression: CKN's kept;
    - anything else: CKN's kept, listed as a conflict. '''
    old = dict(edge)
    if edge['type'] == new['type']:
        if edge['type'] == PTM and edge['effect'] not in ('act', 'inh'):
            set_effect(edge, new['effect'])
        else:
            set_effect(edge, combined_effect(edge['effect'], new['effect']))
    elif (edge['type'] == 'binding' and new['type'] != 'binding') or edge['rank'] == '4':
        edge['type'], edge['isTFregulation'] = new['type'], new['isTFregulation']
        set_effect(edge, new['effect'])
        if edge['rank'] == '4':
            edge['directed'] = new['directed']
    elif not (new['type'] == 'binding' or
              (edge['type'] == 'small RNA interactions' and new['reaction_type'] == REPRESSION)):
        log['conflicts'].append((old, new))
    edge['interactionSources'] = '|'.join(sources(edge) + [new['interactionSources']])
    edge['rank'] = '0'


def pss_edges(rows, wanted):
    ''' The gene network rows of the wanted reactions as CKN edges: one per node pair and direction (undirected:
    one per pair), with the reactions of the pair joined '''
    edges = {}
    for r in rows:
        if r['reaction_id'] not in wanted:
            continue
        a, b = node_id(r['source']), node_id(r['target'])
        directed = r['directed'] == 'True'
        key = (a, b, True) if directed else (*sorted((a, b)), False)
        effect = EFFECTS.get(r['interaction'], 'unk')
        source = f"skm:{r['reaction_id']}"
        if key in edges:
            e = edges[key]
            if source not in sources(e):
                e['interactionSources'] += '|' + source
            set_effect(e, combined_effect(e['effect'], effect))
            continue
        ckn_type = PSS_TYPES.get(r['reaction_type'], 'other')
        edges[key] = {
            'source': key[0], 'interaction': r['interaction'], 'target': key[1], 'directed': directed, 'rank': '0',
            'effect': effect, 'type': ckn_type, 'species': 'ath',
            'isTFregulation': '1' if ckn_type == TF_REGULATION else '0', 'interactionSources': source,
            'reaction_type': r['reaction_type'],
        }
    return list(edges.values())


def node_species(node, pss_node):
    ''' As in v2: metabolites none, foreign nodes 'foreign', complexes 'ath' or 'ath/foreign', the others 'ath' '''
    node_type = node['node_type']
    if node_type in ('Metabolite', 'MetaboliteFamily'):
        return ''
    if node_type.startswith('Foreign'):
        return 'foreign'
    if node_type == 'Complex':
        components = [pss_node.get(c, {}).get('node_type', '') for c in node['components'].split('|') if c]
        return 'ath/foreign' if any(t.startswith('Foreign') for t in components) else 'ath'
    return 'ath'


def opposite_signs(edges, gn_edges):
    ''' Reactions with an edge (with no other PSS source) that is an inhibition where the gene network's isn't, or
    the other way round, e.g. rx00016, an inhibition in 2023, now an activation. (v2's binding partners are a
    positive influence where the gene network's are unknown: not a change, the 2023 conversion.) '''
    signs = {}
    for r in gn_edges:
        pair = frozenset((node_id(r['source']), node_id(r['target'])))
        signs.setdefault((r['reaction_id'], pair), set()).add(r['interaction'])
    changed = set()
    for e in edges:
        rx = reactions(e)
        if len(rx) != 1:
            continue
        rx = rx.pop()
        gn = signs.get((rx, frozenset((e['source'], e['target']))))
        negative = e['interaction'] == 'negative-influence'
        if gn and (negative != ('negative-influence' in gn)) and (negative or gn == {'negative-influence'}):
            changed.add(rx)
    return changed


def source_rank(edge, edges):
    ''' The rank of an edge that has lost its PSS sources (rank 0 is PSS's): the rank its other sources have
    most often on the other edges, the best of them (e.g. 10.1093/molbev/msv058: 1, 10.1093/nar/gkp818: 3) '''
    ranks = {s: Counter() for s in sources(edge)}
    for e in edges:
        if e['rank'] != '0':
            for s in sources(e):
                if s in ranks:
                    ranks[s][e['rank']] += 1
    best = [c.most_common(1)[0][0] for c in ranks.values() if c]
    return min(best) if best else edge['rank']


def update_pss(nodes, edges, gn_edges, gn_nodes, deprecated):
    ''' The PSS part of CKN as in the PSS gene network (see the module); returns nodes, edges and a log '''
    log = {'renamed': Counter(), 'stripped': [], 'removed': [], 'merged': [], 'conflicts': []}

    # 1. PSS names
    for e in edges:
        for end in ('source', 'target'):
            if e[end] in PSS_NAMES:
                log['renamed'][e[end]] += 1
                e[end] = PSS_NAMES[e[end]]
    for n in nodes:
        n['id'] = PSS_NAMES.get(n['id'], n['id'])

    # 2. reactions to remove: no longer in the gene network, or with v2 edges the gene network doesn't have
    # (node pairs without direction: v2 has some binding edges undirected that the gene network has directed)
    gn_pairs = {}
    for r in gn_edges:
        gn_pairs.setdefault(r['reaction_id'], set()).add(frozenset((node_id(r['source']), node_id(r['target']))))
    ckn_pairs = {}
    for e in edges:
        for rx in reactions(e):
            ckn_pairs.setdefault(rx, set()).add(frozenset((e['source'], e['target'])))
    gone = {rx for rx in ckn_pairs if rx not in gn_pairs}
    changed = {rx for rx in ckn_pairs if rx in gn_pairs and not ckn_pairs[rx] <= gn_pairs[rx]}
    changed |= opposite_signs(edges, gn_edges) - gone
    log['gone'] = {rx: 'deprecated' if rx in deprecated else 'not public Arabidopsis' for rx in sorted(gone)}
    log['changed'] = sorted(changed)
    remove = gone | changed

    kept = []
    for e in edges:
        drop = reactions(e) & remove
        if not drop:
            kept.append(e)
            continue
        rest = [s for s in sources(e) if not ((m := PSS_SOURCE.match(s)) and m.group(1) in drop)]
        if not rest:
            log['removed'].append(e)
            continue
        log['stripped'].append((dict(e), sorted(drop)))
        e['interactionSources'] = '|'.join(rest)
        kept.append(e)
    edges = kept

    # 3. the reactions CKN doesn't have, from the gene network
    have = {rx for e in edges for rx in reactions(e)}
    log['added_reactions'] = sorted(set(gn_pairs) - have)
    index = {}
    for i, e in enumerate(edges):
        for pair in pairs(e):
            index.setdefault(pair, i)
    added = []
    for new in pss_edges(gn_edges, set(log['added_reactions'])):
        i = index.get((new['source'], new['target']))
        if i is None and not new['directed']:
            i = index.get((new['target'], new['source']))
        if i is None:
            new.pop('reaction_type')
            added.append(new)
            continue
        log['merged'].append((dict(edges[i]), new['interactionSources']))
        merge_pss_edge(edges[i], new, log)
    edges += added
    log['added'] = added
    log['no_pss_source'] = [e for e in edges if e['rank'] == '0' and not reactions(e)]
    for e in log['no_pss_source']:
        e['rank'] = source_rank(e, edges)

    # 4. nodes: non-gene PSS nodes with their PSS annotation (new ones added); those without edges removed
    pss_node = {node_id(n['id']): n for n in gn_nodes}
    used = {x for e in edges for x in (e['source'], e['target'])}
    by_id = {n['id']: n for n in nodes}
    for x in sorted(used):
        p = pss_node.get(x)
        if p is None or GENE.match(x):
            continue
        n = by_id.get(x)
        if n is None:
            n = {'id': x, 'locus_type': '', 'TAIR': '', 'mapman': '', 'note': '', 'tissue': 'not assigned'}
            nodes.append(n)
            by_id[x] = n
        n.update({'node_type': p['node_type'], 'locus_type': '', 'species': node_species(p, pss_node),
                  'display_label': p['display_label'] or x, 'short_name': p['short_name'],
                  'synonyms': p['synonyms'], 'description': p['description']})
    log['orphans'] = [n['id'] for n in nodes if n['id'] not in used]
    nodes = [n for n in nodes if n['id'] in used]
    return nodes, edges, log


def edge_key(edge):
    ''' (source, target, directed); an undirected edge's nodes sorted '''
    a, b = edge['source'], edge['target']
    return (a, b, True) if edge['directed'] else (*sorted((a, b)), False)


def check(nodes, edges, gn_edges, deprecated):
    ''' Every edge's nodes in node-annot, no deprecated reaction, the same reactions as the gene network, no
    edge twice (two edges between the same nodes are fine if one is directed and the other not, or in opposite
    directions: e.g. binding and transcription factor regulation) '''
    ids = Counter(n['id'] for n in nodes)
    assert all(c == 1 for c in ids.values()), f'nodes twice: {[i for i, c in ids.items() if c > 1][:10]}'
    missing = {x for e in edges for x in (e['source'], e['target'])} - set(ids)
    assert not missing, f'nodes without annotation: {sorted(missing)[:10]}'
    in_ckn = {rx for e in edges for rx in reactions(e)}
    in_gn = {r['reaction_id'] for r in gn_edges}
    assert not in_ckn & deprecated, f'deprecated reactions: {sorted(in_ckn & deprecated)}'
    assert in_ckn == in_gn, f'only in CKN: {sorted(in_ckn - in_gn)}; only in PSS: {sorted(in_gn - in_ckn)}'
    twice = [k for k, c in Counter(edge_key(e) for e in edges).items() if c > 1]
    assert not twice, f'edges twice: {twice[:5]}'
    pss = {x for e in edges if reactions(e) for x in (e['source'], e['target'])}
    print(f"checks passed: {len(in_ckn)} PSS reactions, {sum(1 for e in edges if reactions(e))} PSS edges, "
          f"{len(pss)} nodes in PSS edges")


def report(log):
    print(f"renamed to PSS names: {len(log['renamed'])} nodes ({sum(log['renamed'].values())} edge ends)")
    print(f"reactions removed: {len(log['gone'])}: " + ', '.join(f'{rx} ({why})' for rx, why in log['gone'].items()))
    print(f"reactions replaced (v2 edges PSS no longer has): {len(log['changed'])}: {', '.join(log['changed'])}")
    print(f"edges removed: {len(log['removed'])}; edges kept without a removed reaction's source: "
          f"{len(log['stripped'])}")
    for e, drop in log['stripped']:
        print(f"  {e['source']} -> {e['target']} ({e['type']}, rank {e['rank']}): -{','.join(drop)}: "
              f"{e['interactionSources']}")
    print(f"reactions added: {len(log['added_reactions'])}; edges added: {len(log['added'])}; merged into "
          f"existing edges: {len(log['merged'])}")
    print(f"type conflicts (CKN's kept): {len(log['conflicts'])}")
    for old, new in log['conflicts']:
        print(f"  {old['source']} -> {old['target']}: CKN {old['type']}/{old['effect']} (rank {old['rank']}), "
              f"PSS {new['type']}/{new['effect']} {new['interactionSources']}")
    print(f"rank 0 edges left without a PSS source (ranked by their other sources): {len(log['no_pss_source'])}")
    for e in log['no_pss_source']:
        print(f"  {e['source']} -> {e['target']} ({e['type']}/{e['effect']}): rank {e['rank']}: "
              f"{e['interactionSources']}")
    print(f"nodes removed (no edges left): {len(log['orphans'])}: {', '.join(log['orphans'])}")


def read_tsv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE))


def main(folder_in, gn_edges_path, gn_nodes_path, deprecated_path, folder_out):
    folder_in, folder_out = Path(folder_in), Path(folder_out)
    folder_out.mkdir(parents=True, exist_ok=True)

    with open(MAPMAN, newline='') as f:
        mapman = {r['gene']: r['mapman'] for r in csv.DictReader(f, delimiter='\t')}
    nodes, edges = ckn_v2(folder_in, mapman)

    gn_edges, gn_nodes = read_tsv(gn_edges_path), read_tsv(gn_nodes_path)
    deprecated = set(Path(deprecated_path).read_text().split())
    nodes, edges, log = update_pss(nodes, edges, gn_edges, gn_nodes, deprecated)
    report(log)
    check(nodes, edges, gn_edges, deprecated)

    write(folder_out / f'AtCKN-{VERSION_OUT}_node-annot.tsv.gz', NODE_COLUMNS, nodes)
    for (node_type, locus_type), count in sorted(Counter((n['node_type'], n['locus_type']) for n in nodes).items()):
        print(f'  {node_type:17} {locus_type or "-":29} {count}')
    write(folder_out / f'AtCKN-{VERSION_OUT}.tsv.gz', EDGE_COLUMNS, edges)


if __name__ == '__main__':
    main(*sys.argv[1:6])
