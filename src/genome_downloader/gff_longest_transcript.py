from pybedtools import BedTool
from pybedtools.featurefuncs import gff2bed
import re
import tempfile
import argparse
import pandas as pd
import numpy as np


def filter_transcript(feature):
    return not re.search("_", feature.fields[0]) and feature.fields[2]=="transcript"


def find_longest_transcript(input, output, clip_start=0, clip_end=0, clip_strand_specific=False, output_strand_specific=False):
    tmp_file = tempfile.NamedTemporaryFile(delete=False)

    # Load GFF file
    annotations = BedTool(input)

    # Select only transcripts and convert to BED format
    transcripts = annotations.filter(filter_transcript).\
        each(gff2bed, name_field="gene_id").sort().\
        saveas(tmp_file.name).\
        to_dataframe().\
        assign(length=lambda x: x.end - x.start + 1)

    # Select longest transcript per gene
    transcripts_longest = transcripts.loc[transcripts.reset_index().groupby(['name'])['length'].idxmax()].\
        drop(['length'], axis=1)

    # Clip at the end or the beginning of the gene (can be strand specific)
    pos_strand = (transcripts_longest["strand"] == "+").values | np.invert(clip_strand_specific)
    if clip_start > 0:
        transcripts_longest.loc[pos_strand,"start"] = transcripts_longest.loc[pos_strand,"start"] + clip_start
        transcripts_longest.loc[~pos_strand,"end"] = transcripts_longest.loc[~pos_strand,"end"] - clip_end
    if clip_end > 0:
        transcripts_longest.loc[pos_strand,"end"] = transcripts_longest.loc[pos_strand,"end"] - clip_end
        transcripts_longest.loc[~pos_strand,"start"] = transcripts_longest.loc[~pos_strand,"start"] + clip_start

    # Notify about genes with negative length
    transcripts_negative_length = transcripts_longest.query("start >= end").name.values
    if len(transcripts_negative_length) > 0:
        transcripts_longest = transcripts_longest.query("end > start")
        print("Removing transcripts with negative length from the output file: {}".format(", ".join(transcripts_negative_length)))

    # make a copy of genes and reverse strand
    transcripts_longest_reversed = transcripts_longest.copy()
    transcripts_longest_reversed["strand"] = ["+" if s=="-" else "-" for s in transcripts_longest_reversed["strand"].values]
    transcripts_longest_reversed["name"] = transcripts_longest_reversed["name"] + "_rev"
    transcripts_longest_stack = pd.concat([transcripts_longest, transcripts_longest_reversed]).\
        sort_values(["chrom", "start", "name", "strand"]).\
        reset_index(drop=True)

    # Save final data frame to file
    if not output_strand_specific:
        transcripts_longest_stack.to_csv(output, index=False, sep="\t", header=False)
        pass
    else:
        for s, s_name in {'+': 'pos', '-': 'neg'}.items():
            transcripts_longest_strand = transcripts_longest_stack.query("strand == @s")
            transcripts_longest_strand.to_csv("{}_{}".format(output, s_name), index=False, sep="\t", header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find longest transcript in GFF file write it to separate file')
    parser.add_argument('input', help='Input GFF file')
    parser.add_argument('output', help='Output bed file with longest transcripts')
    parser.add_argument('--output-strand', dest="output_strand_specific", action='store_true', help='If set outputs two files with .neg and .pos suffix')
    parser.add_argument('--clip-start', dest="clip_start", default=0, type=int, help='Clip N base pairs at the beginning of the longest transcript')
    parser.add_argument('--clip-end', dest="clip_end", default=0, type=int, help='Clip N base pairs at the end of the longest transcript')
    parser.add_argument('--clip-strand', dest="clip_strand_specific", action='store_true', help='Is clipping done strand specifically? If set then negative strand is clipped at reverse ends')

    args = parser.parse_args()
    find_longest_transcript(args.input, args.output, clip_start=args.clip_start, clip_end=args.clip_end, clip_strand_specific=args.clip_strand_specific, output_strand_specific=args.output_strand_specific)
