"""
data writing and visualization.
"""
# BUILT-INS
import os
import math

#LIBS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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


def get_primers_from_amp_dic(amp, amplicons, left_primer_candidates, right_primer_candidates):
    """
    get primers from amplicon dictionary
    """

    left_name = amplicons[amp][2]
    left_primer = left_primer_candidates[left_name]
    right_name = amplicons[amp][3]
    right_primer = right_primer_candidates[right_name]

    return (left_name, left_primer),(right_name, right_primer)


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
        for idx, amp in enumerate(amplicon_scheme):
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
                amplicons[amp][0],  # start
                amplicons[amp][1],  # stop
                new_name,
                pool,
                sep = "\t",
                file = bed
            )
            # get primers
            left, right = get_primers_from_amp_dic(
                amp,
                amplicons,
                left_primer_candidates,
                right_primer_candidates
            )
            # write primer assignments tabular file
            print(
                left[0],
                right[0],
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


def entropy(pos, states):
    """calculate the entropy on the basis of a string and a list of unique_chars"""
    max_ent = -1/(math.log(1/float(states), 10))
    # only a rough normalization factor, not needed, but gives more
    # beautiful plots
    unique_chars = list(set(pos))
    ent = 0.0
    if len(pos) < 2:
        return ent
    for char in unique_chars:
    # calculate the entropy at the particular position
        freq = pos.count(char)
        if freq > 0:
            freq = float(freq)/float(len(pos))
            ent += freq*math.log(freq, 10)
    if ent == 0:
        return ent
    else:
        return -ent*max_ent
        # max_ent is the normalization


def alignment_entropy(alignment_cleaned):
    """
    calculate the entropy for every position in an alignment.
    return pandas df.
    """
    position = list()
    entropys = list()
    # iterate over alignment positions and the sequences
    for nuc_pos in range(0, len(alignment_cleaned[0][1])):
        pos = []
        entropy_temp = []
        for seq_number in range(0, len(alignment_cleaned)):
            pos.append(alignment_cleaned[seq_number][1][nuc_pos])
        entropys.append(entropy(pos, 4))
        position.append(nuc_pos)
    # create df
    entropy_df = pd.DataFrame()
    entropy_df["position"] = position
    entropy_df["entropy"] = entropys
    entropy_df["average"] = entropy_df["entropy"].rolling(10, center=True).mean()

    return entropy_df


def get_primer_list(amplicon_scheme, amplicons, left_primer_candidates, right_primer_candidates):
    """
    get primers as list
    """
    primer_list = []
    for amp in amplicon_scheme:
        left, right =get_primers_from_amp_dic(
            amp,
            amplicons,
            left_primer_candidates,
            right_primer_candidates
        )
        primer_list.append(left)
        primer_list.append(right)

    return primer_list


def varvamp_plot(dir, threshold, alignment_cleaned, conserved_regions, amplicon_scheme, amplicons, left_primer_candidates, right_primer_candidates):
    """
    creates overview plot for the amplicon design
    and per base coverage plots
    """
    # first: create the overview plot
    # - create pdf name
    name = "amplicon_plot.pdf"
    out = os.path.join(dir, name)
    # - create entropy df
    entropy_df = alignment_entropy(alignment_cleaned)
    # - get primer list
    primers = get_primer_list(
        amplicon_scheme,
        amplicons,
        left_primer_candidates,
        right_primer_candidates
    )
    # - ini figure
    fig, ax = plt.subplots(
        2,
        1,
        figsize=[20, 6],
        squeeze=True,
        sharex=True,
        gridspec_kw={'height_ratios': [6, 1]}
    )
    fig.subplots_adjust(hspace=0)
    # - entropy plot
    ax[0].fill_between(
        entropy_df["position"],
        entropy_df["entropy"],
        color="gainsboro",
        label="entropy"
    )
    ax[0].plot(
        entropy_df["position"],
        entropy_df["average"],
        color="black",
        label="average entropy",
        linewidth=0.5
    )
    ax[0].set_ylim((0, 1))
    ax[0].set_ylabel("entropy")
    ax[0].set_title("final amplicon design")
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    # - conserved regions plot
    for region in conserved_regions:
        ax[1].hlines(
            [1],
            region[0],
            region[1],
            linewidth=15,
            color="darkorange"
        )
    # - conserved legend
    ax[1].hlines(
        [1],
        conserved_regions[0][1],
        conserved_regions[0][1],
        label="possible primer regions",
        linewidth=15, color="darkorange"
    )
    # - amplicon plot pool 1
    for i in range(0, len(amplicon_scheme), 2):
        amp = amplicons[amplicon_scheme[i]]
        ax[1].hlines(
            0.8,
            amp[0],
            amp[1],
            linewidth=5
        )
    # - amplicon plot pool 2
    for i in range(1, len(amplicon_scheme), 2):
        amp = amplicons[amplicon_scheme[i]]
        ax[1].hlines(
            0.7,
            amp[0],
            amp[1],
            linewidth=5
        )
    # - amplicon legend
    amp = amplicons[amplicon_scheme[0]]
    ax[1].hlines(
            0.8,
            amp[0],
            amp[1],
            linewidth=5,
            label="amplicons"
        )
    # - primer plots
    for i in range(1, len(primers), 4):
        ax[1].hlines(0.8, primers[i][1][1], primers[i][1][2], linewidth=5, color="red")
    for i in range(2, len(primers), 4):
        ax[1].hlines(0.7, primers[i][1][1], primers[i][1][2], linewidth=5, color="red")
    for i in range(3, len(primers), 4):
        ax[1].hlines(0.7, primers[i][1][1], primers[i][1][2], linewidth=5, color="red")
    for i in range(4, len(primers), 4):
        ax[1].hlines(0.8, primers[i][1][1], primers[i][1][2], linewidth=5, color="red")
    # - primer legend
    ax[1].hlines(0.8, primers[0][1][1], primers[0][1][2], linewidth=5, color="red", label="primers")
    # - finalize
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].axes.get_yaxis().set_visible(False)
    ax[1].set_xlabel("genome position")
    ax[1].set_ylim((0.6, 1))
    fig.legend()
    # - save fig
    fig.savefig(out)

    # second plot: per base coverage
    # - ini name
    name = "per_base_mismatches.pdf"
    out = os.path.join(dir, name)
    # - ini multi pdf
    with PdfPages(out) as pdf:
        # - always print 4 primers to one page
        for i in range(0, len(primers), 4):
            # - ini figure
            primers_temp = primers[i:i+4]
            fig, ax = plt.subplots(len(primers_temp), figsize=(12, len(primers_temp)*4), squeeze=True)
            fig.suptitle("Per base mismatches", fontsize=18)
            fig.tight_layout(rect=[0.05, 0.05, 1, 0.98])
            fig.subplots_adjust(hspace=0.5)
            # - plotting
            for idx, primer in enumerate(primers_temp):
                x = [pos+primer[1][1] for pos in range(0, len(primer[1][4]))]
                ax[idx].bar(x, primer[1][4],
                            color='lightgrey', edgecolor='black')
                ax[idx].set_title(primer[0], loc = "left")
                ax[idx].xaxis.set_ticks(np.arange(primer[1][1], primer[1][1]+len(x), 1))
                ax[idx].xaxis.set_ticklabels(x, rotation = 45)
                ax[idx].set_ylabel(ylabel="% of sequences")
                ax[idx].set_xlabel("position")
                ax[idx].set_ylim(0, 1-threshold)
            # - to pdf
            pdf.savefig(fig, bbox_inches='tight')
