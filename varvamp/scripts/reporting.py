#!/usr/bin/env python
"""
            INFO
------------------------------
Contains functions for writing the raw data and visualization.

           LICENCE
-------------------------------
todo

"""
# BUILT-INS
import os

# varVAMP
from varvamp.scripts import primers


def write_fasta(dir, seq_id, seq):
    """
    write fasta files
    """
    name = seq_id + ".fasta"
    out = os.path.join(dir, name)
    with open(out, 'w') as o:
        print(f">{seq_id}\n{seq}", file=o)


def write_alignment(dir, alignment):
    """
    write alignment to file
    """
    name = "alignment_cleaned.fasta"
    out = os.path.join(dir, name)
    with open(out, "w") as o:
        for seq in alignment:
            print(f">{seq[0]}\n{seq[1]}", file=o)


def conserved_to_bed(conserved_regions, dir):
    """
    write conserved regions as bed file
    """
    counter = 0
    outfile = dir+"conserved_regions.bed"
    with open(outfile, 'w') as o:
        for region in conserved_regions:
            print(
                "ambiguous_consensus",
                region[0],
                region[1],
                "region_"+str(counter),
                sep="\t",
                file=o
            )
            counter += 1


def primers_to_bed(outfile, primer_name, primer_prop, direction):
    """
    write primers as bed file
    """
    with open(outfile, 'a') as o:
        print(
            "ambiguous_consensus",
            primer_prop[1],
            primer_prop[2],
            primer_name,
            primer_prop[3],
            direction,
            sep="\t",
            file=o
        )


def write_all_primers(dir, left_primer_candidates, right_primer_candidates):
    """
    write all primers that varVAMP designed as bed file
    """
    outfile = dir + "all_primers.bed"

    for direction, primer_candidates in [("+", left_primer_candidates), ("-", right_primer_candidates)]:
        for primer in primer_candidates:
            primers_to_bed(outfile, primer, primer_candidates[primer], direction)


def write_scheme_to_files(dir, amplicon_scheme, amplicons, ambiguous_consensus, left_primer_candidates, right_primer_candidates):
    """
    write all relevant bed files and a tsv file with all primer stats
    """
    # ini
    tsv_file = os.path.join(dir,"primers.tsv")
    primer_bed_file =  os.path.join(dir, "primers.bed")
    amplicon_bed_file =  os.path.join(dir, "amplicons.bed")
    tabular_file =  os.path.join(dir, "primer_to_amplicon_assignment.tabular")
    # counter for new amplicon name
    counter = 0

    # open files to write
    with open(tsv_file, "w") as tsv, open(amplicon_bed_file, "w") as bed, open(tabular_file, "w") as tabular:
        # write header for primer tsv
        print(
            "amlicon_name\tprimer_name\tpool\tseq\tsize\tgc\ttemp\tscore",
            file=tsv
        )
        for idx, final_amp in enumerate(amplicon_scheme):
            # assign pools
            if idx % 2:
                pool = 1
            else:
                pool = 0
            # give a new amplicon name
            new_name = "amplicon_" + str(counter)
            counter += 1
            # write amplicon bed
            print(
                "ambiguous_consensus",
                amplicons[final_amp][0],  # start
                amplicons[final_amp][1],  # stop
                new_name,
                pool,
                sep = "\t",
                file = bed
            )
            # get primers
            left_name = amplicons[final_amp][2]
            left_primer = left_primer_candidates[left_name]
            left = (left_name, left_primer)
            right_name = amplicons[final_amp][3]
            right_primer = right_primer_candidates[right_name]
            right = (right_name, right_primer)
            # write primer assignments tabular file
            print(
                left_name,
                right_name,
                sep = "\t",
                file = tabular
            )
            # write primer tsv and primer bed
            for direction, primer in [("+", left), ("-",right)]:
                # get the ambiguous seq and rev complement for RIGHT primers
                seq = ambiguous_consensus[primer[1][1]:primer[1][2]]
                if direction == "-":
                    seq = primers.rev_complement(seq)
                # calc primer parameters
                gc = round(primers.calc_gc(primer[1][0]), 1)
                temp = round(primers.calc_temp(primer[1][0]), 1)
                score = round(primer[1][3], 1)
                # write tsv file
                print(
                    new_name,
                    primer[0],
                    pool,
                    seq,
                    len(primer[1][0]),
                    gc,
                    temp,
                    score,
                    sep = "\t",
                    file = tsv
                )
                # write primer bed file
                primers_to_bed(primer_bed_file, primer[0], primer[1], direction)
