"""
data writing and visualization.
"""
# BUILT-INS
import os
import math
import itertools

# LIBS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# varVAMP
from varvamp.scripts import primers
from varvamp.scripts import config


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


def write_conserved_to_bed(conserved_regions, dir):
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


def write_primers_to_bed(outfile, primer_name, primer_properties, direction):
    """
    write primers as bed file
    """
    with open(outfile, 'a') as o:
        print(
            "ambiguous_consensus",
            primer_properties[1],  # start
            primer_properties[2],  # stop
            primer_name,
            round(primer_properties[3], 1),  # score
            direction,
            sep="\t",
            file=o
        )


def write_all_primers(dir, all_primers):
    """
    write all primers that varVAMP designed as bed file
    """
    outfile = dir + "all_primers.bed"

    for direction in all_primers:
        for primer in all_primers[direction]:
            write_primers_to_bed(outfile, primer, all_primers[direction][primer], direction)


def get_permutations(seq):
    """
    get all permutations of an ambiguous sequence. needed to
    correctly report the gc and the temperature.
    """
    groups = itertools.groupby(seq, lambda char: char not in config.ambig_nucs)
    splits = []
    for b, group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(config.ambig_nucs[nuc])
    return[''.join(p) for p in itertools.product(*splits)]


def write_scheme_to_files(dir, amplicon_scheme, ambiguous_consensus):
    """
    write all relevant bed files and a tsv file with all primer stats
    """
    # ini
    tsv_file = os.path.join(dir, "primers.tsv")
    primer_bed_file = os.path.join(dir, "primers.bed")
    amplicon_bed_file = os.path.join(dir, "amplicons.bed")
    tabular_file = os.path.join(dir, "primer_to_amplicon_assignment.tabular")

    counter = 0

    # open files to write
    with open(tsv_file, "w") as tsv, open(amplicon_bed_file, "w") as bed, open(tabular_file, "w") as tabular:
        # write header for primer tsv
        print(
            "amlicon_name\tprimer_name\tpool\tstart\tstop\tseq\tsize\tgc_best\ttemp_best\tmean_gc\tmean_temp\tscore",
            file=tsv
        )
        counter = 0

        for pool in amplicon_scheme:
            for amp in amplicon_scheme[pool]:
                # give a new amplicon name
                new_name = f"amplicon_{str(counter)}"
                counter += 1
                # get left and right primers and their names
                primer_names = list(amplicon_scheme[pool][amp].keys())
                left = (primer_names[0], amplicon_scheme[pool][amp][primer_names[0]])
                right = (primer_names[1], amplicon_scheme[pool][amp][primer_names[1]])

                # write amplicon bed
                print("ambiguous_consensus", left[1][1], right[1][2], new_name, pool, sep="\t", file=bed)
                # write primer assignments tabular file
                print(left[0], right[0], sep="\t", file=tabular)

                # write primer tsv and primer bed
                for direction, primer in [("+", left), ("-", right)]:
                    seq = ambiguous_consensus[primer[1][1]:primer[1][2]]
                    if direction == "-":
                        seq = primers.rev_complement(seq)
                    # calc primer parameters for all permutations
                    gc = 0
                    temp = 0
                    permutations = get_permutations(seq)
                    for permutation in permutations:
                        gc += primers.calc_gc(permutation)
                        temp += primers.calc_temp(permutation)
                    # write tsv file
                    print(
                        new_name,
                        primer[0],
                        pool,
                        primer[1][1],
                        primer[1][2],
                        seq,
                        len(primer[1][0]),
                        round(primers.calc_gc(primer[1][0]), 1),
                        round(primers.calc_temp(primer[1][0]), 1),
                        round(gc/len(permutations), 1),
                        round(temp/len(permutations), 1),
                        round(primer[1][3], 1),
                        sep="\t",
                        file=tsv
                    )
                    # write primer bed file
                    write_primers_to_bed(primer_bed_file, primer[0], primer[1], direction)


def write_dimers(dir, not_solved):
    """
    write dimers for which no replacement was found to file
    """
    tsv_file = os.path.join(dir, "unsolvable_primer_dimers.tsv")
    print(
        "pool\tprimer_name_1\tprimer_name_2\tdimer melting temp",
        file=tsv_file
    )
    for dimers in not_solved:
        print(
            dimers[0][0],
            dimers[0][2],
            dimers[1][2],
            round(primers.calc_dimer(dimers[0][3][0], dimers[1][3][0]).tm, 1),
            sep="\t",
            file=tsv_file
        )


def entropy(pos, states):
    """calculate the entropy on the basis of a string and a list of unique_chars"""
    max_ent = -1/(math.log(1/float(states), 10))
    # only a rough normalization factor, not needed, but gives more
    # beautiful plots
    unique_chars = list(set(pos))
    ent = 0.0
    if len(pos) < 2:
        return ent
    # calculate the entropy at the particular position
    for char in unique_chars:
        freq = pos.count(char)
        if freq > 0:
            freq = float(freq)/float(len(pos))
            ent += freq*math.log(freq, 50)
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


def varvamp_plot(dir, threshold, alignment_cleaned, conserved_regions, all_primers, amplicon_scheme):
    """
    creates overview plot for the amplicon design
    and per base coverage plots
    """

    amplicon_primers = []
    # first plot: overview
    # - create pdf name
    name = "amplicon_plot.pdf"
    out = os.path.join(dir, name)
    # - create entropy df
    entropy_df = alignment_entropy(alignment_cleaned)

    # - ini figure
    fig, ax = plt.subplots(2, 1, figsize=[22, 6], squeeze=True, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    fig.subplots_adjust(hspace=0)
    # - entropy plot
    ax[0].fill_between(entropy_df["position"], entropy_df["entropy"], color="gainsboro", label="entropy")
    ax[0].plot(entropy_df["position"], entropy_df["average"], color="black", label="average entropy", linewidth=0.5)
    ax[0].set_ylim((0, 1))
    ax[0].set_xlim(0, max(entropy_df["position"]))
    ax[0].set_ylabel("alignment entropy")
    ax[0].set_title("final amplicon design")
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)

    # - conserved regions plot
    for region in conserved_regions:
        ax[1].hlines([1], region[0], region[1], linewidth=15, color="darkorange")
    # - conserved legend
    ax[1].hlines([1], conserved_regions[0][1], conserved_regions[0][1], label="possible primer regions", linewidth=5, color="darkorange")

    # - all primer plot
    for direction in all_primers:
        if direction == "-":
            primer_position = 0.85
            primer_color = "darkgrey"
            primer_label = "all right primers"
        elif direction == "+":
            primer_position = 0.8
            primer_color = "dimgrey"
            primer_label = "all left primers"
        for primer in all_primers[direction]:
            ax[1].hlines(primer_position, all_primers[direction][primer][1], all_primers[direction][primer][2], linewidth=5, color=primer_color)
    # - legend
        ax[1].hlines(primer_position, all_primers[direction][primer][1], all_primers[direction][primer][2], linewidth=5, color=primer_color, label=primer_label)

    # - amplicon, text and primer plot
    counter = 0
    for pool in amplicon_scheme:
        for amp in amplicon_scheme[pool]:
            if pool == 0:
                position_amp = 0.7
                position_text = 0.6
            elif pool == 1:
                position_amp = 0.6
                position_text = 0.65
            primer_names = list(amplicon_scheme[pool][amp].keys())
            left = amplicon_scheme[pool][amp][primer_names[0]]
            right = amplicon_scheme[pool][amp][primer_names[1]]
            # amplicons
            ax[1].hlines(position_amp, left[1], right[2], linewidth=5)
            # text
            ax[1].text(right[2] - (right[2]-left[1])/2, position_text, str(counter), fontsize=8)
            # primers
            ax[1].hlines(position_amp, left[1], left[2], linewidth=5, color="red")
            ax[1].hlines(position_amp, right[1], right[2], linewidth=5, color="red")

            counter += 1
            # remember primers and names as they are needed for the last plot
            amplicon_primers.append((primer_names[0], left))
            amplicon_primers.append((primer_names[1], right))

    # - legends
    ax[1].hlines(position_amp, left[1]+config.PRIMER_SIZES[1], right[2]-config.PRIMER_SIZES[1], linewidth=5, label="amplicons")
    ax[1].hlines(position_amp, left[1], left[2], linewidth=5, color="red", label="primers")

    # - finalize
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].axes.get_yaxis().set_visible(False)
    ax[1].set_xlabel("genome position")
    ax[1].set_ylim((0.5, 1))
    fig.legend(loc=(0.83, 0.7))
    # - save fig
    fig.savefig(out, bbox_inches='tight')

    # second plot: per base primer mismatches
    # - ini name
    name = "per_base_mismatches.pdf"
    out = os.path.join(dir, name)
    # - ini multi pdf
    with PdfPages(out) as pdf:
        # - always print 4 primers to one page
        for i in range(0, len(amplicon_primers), 4):
            # - ini figure
            primers_temp = amplicon_primers[i:i+4]
            fig, ax = plt.subplots(len(primers_temp), figsize=(12, len(primers_temp)*4), squeeze=True)
            fig.suptitle("Per base mismatches", fontsize=18)
            fig.tight_layout(rect=[0.05, 0.05, 1, 0.98])
            fig.subplots_adjust(hspace=0.5)
            # - plotting
            for idx, primer in enumerate(primers_temp):
                x = [pos+primer[1][1] for pos in range(0, len(primer[1][4]))]
                ax[idx].bar(x, primer[1][4], color='lightgrey', edgecolor='black')
                ax[idx].set_title(primer[0], loc="left")
                ax[idx].xaxis.set_ticks(np.arange(primer[1][1], primer[1][1]+len(x), 1))
                ax[idx].xaxis.set_ticklabels(x, rotation=45)
                ax[idx].set_ylabel(ylabel="% of sequences")
                ax[idx].set_xlabel("position")
                ax[idx].set_ylim(0, 1-threshold)
            # - to pdf
            pdf.savefig(fig, bbox_inches='tight')
