#!/usr/bin/env python3
"""Convert an anvi'o-exported minimal GFF3 (+ matching assembly FASTA) into a real,
multi-record GenBank (.gbff) file, with each CDS's locus_tag set to its GFF3 ID
(== the FAA header == ampcombi's CDS_id) so ampcombi's parse_gbks.py can match
predicted AMPs back to a contig via locus_tag.

Usage: gff3_to_gbk.py --gff3 sample.gff3 --fasta sample.fa --output sample.gbff --sample-id sample
"""

import argparse

from BCBio import GFF
from Bio import SeqIO


def convert(gff3_path, fasta_path, output_path, sample_id):
    base_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    records = []
    for record in GFF.parse(gff3_path, base_dict=base_dict):
        record.annotations["molecule_type"] = "DNA"
        for feature in record.features:
            if feature.type == "CDS" and "ID" in feature.qualifiers:
                feature.qualifiers["locus_tag"] = feature.qualifiers["ID"]
        records.append(record)

    with open(output_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")

    n_cds = sum(1 for r in records for f in r.features if f.type == "CDS")
    print("{}: wrote {} contig record(s), {} CDS feature(s) -> {}".format(
        sample_id, len(records), n_cds, output_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()
    convert(args.gff3, args.fasta, args.output, args.sample_id)
