'''
GFF + params --> file with peak targets for mapping
'''

import re
from collections import defaultdict

from pybedtools import BedTool
from pybedtools.featurefuncs import five_prime


#
def is_valid_bedtools(a):
    '''
    Test if a loaded file is a valid BedTool object

    Notes
    =====
    Does not work for BedTool from iterable

    https://github.com/daler/pybedtools/issues/252
    '''

    try:
        ft = a.file_type
        if ft == 'empty':
            print("empty file!")
        else:
            print(f'input is {ft} - OK!')
    except IndexError:
        print('bad input!')


def extract_attribute_gff(x, attr='ID'):
    '''
    Extract a single attribute from GFF attributes column
    '''
    m = re.match(fr".*{attr}=(.*?);", x)
    if m:
        return m.groups()[0]
    print(x)


def generate_regulatory_intervals(bed, upstream=0, downstream=1):
    '''
    Generate regulatory intervals from a gene bed file
    '''
    for gene in bed:
        yield five_prime(gene, upstream=upstream, downstream=downstream)

def update_feature_attrs(feature):
    '''
    Make feature attribute the gene ID and the position of the "TSS"
    '''

    if feature.strand == "-":
        target = feature.end
    else:
        target = feature.start
    attrs = f"Note={target};"

    feature_id = feature.attrs['ID']
    if feature_id:
        # feature.attrs = f"ID={feature_id};Parent={feature_id.split(':')[0]};"
        # feature.attrs = f"ID={feature_id.split(':')[0]};"
        attrs += f"ID={feature_id.split('.')[0]};"

    feature.attrs  = attrs
    return feature


def prepare_peak_targets(gff_file, upstream=2500, downstream=2500, output='peak_targets.tsv', gff_attr="ID"):
    '''
    Prepare peak targets for mapping
    '''


    # Load GFF
    gff = BedTool(gff_file)
    is_valid_bedtools(gff)

    # proteins

    # Extract protein-coding genes

    # Extract gene IDs
    gff = gff.each(update_feature_attrs).saveas()



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Prepare peak targets for mapping peaks to genes')
    parser.add_argument('gff', help='GFF file')
    parser.add_argument('upstream', help='Number upstream nucleotides', default=2500)
    parser.add_argument('downstream', help='Number downstream nucleotides', default=2500)
    parser.add_argument('gff_attr', help='GFF attribute to use as gene ID', default='ID')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()

    prepare_peak_targets(args.gff, args.upstream, args.downstream, args.output, args.gff_attr)
    print('Done')