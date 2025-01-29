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
from varvamp.scripts import logging


def write_fasta(path, file_name, seq_id, seq):
    """
    write fasta files
    """
    name = f"{file_name}.fasta"
    out = os.path.join(path, name)
    with open(out, 'w') as o:
        print(f">{seq_id}\n{seq}", file=o)


def write_alignment(path, alignment):
    """
    write alignment to file
    """
    name = "alignment_cleaned.fasta"
    out = os.path.join(path, name)
    with open(out, "w") as o:
        for seq in alignment:
            print(f">{seq[0]}\n{seq[1]}", file=o)


def write_regions_to_bed(primer_regions, scheme_name, path, mode=None):
    """
    write primer regions as bed file
    """

    if mode == "probe":
        outfile = f"{path}probe_regions.bed"
    else:
        outfile = f"{path}primer_regions.bed"

    with open(outfile, 'w') as o:
        for counter, region in enumerate(primer_regions):
            print(
                f"{scheme_name}_ambiguous_consensus",
                region[0],
                region[1],
                "REGION_"+str(counter),
                sep="\t",
                file=o
            )


def write_primers_to_bed(outfile, scheme_name, primer_name, primer_properties, numeric_value, direction, sequence=None):
    """
    write primers as bed file
    """
    with open(outfile, 'a') as o:
        # write header for primer bed
        if os.path.getsize(outfile) == 0 and sequence is not None:
            print("#chrom\tchromStart\tchromEnd\tprimer-name\tpool\tstrand\tprimer-sequence", file=o)
        data = [f"{scheme_name}_ambiguous_consensus",
            primer_properties[1],  # start
            primer_properties[2],  # stop
            primer_name,
            numeric_value,  # can be pool or score
            direction]
        if sequence is not None:
            data.append(sequence)
        print(
            *data,
            sep="\t",
            file=o
        )


def write_all_primers(path, scheme_name, all_primers):
    """
    write all primers that varVAMP designed as bed file
    """
    outfile = f"{path}all_primers.bed"

    for direction in all_primers:
        for primer in all_primers[direction]:
            write_primers_to_bed(outfile, scheme_name, primer, all_primers[direction][primer], round(all_primers[direction][primer][3], 2), direction)


def get_permutations(seq):
    """
    get all permutations of an ambiguous sequence. needed to
    correctly report the gc and the temperature.
    """
    groups = itertools.groupby(seq, lambda char: char not in config.AMBIG_NUCS)
    splits = []
    for b, group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(config.AMBIG_NUCS[nuc])
    return[''.join(p) for p in itertools.product(*splits)]


def calc_mean_stats(permutations):
    """
    calculate mean gc and temp over all permutations
    """
    gc = 0
    temp = 0

    for permutation in permutations:
        gc += primers.calc_gc(permutation)
        temp += primers.calc_temp(permutation)

    return round(gc/len(permutations), 1), round(temp/len(permutations), 1)


def write_qpcr_to_files(path, final_schemes, ambiguous_consensus, scheme_name, log_file):
    """
    write all relevant bed files and tsv file for the qPCR design
    """

    tsv_file = os.path.join(path, "qpcr_design.tsv")
    tsv_file_2 = os.path.join(path, "qpcr_primers.tsv")
    primer_bed_file = os.path.join(path, "primers.bed")
    amplicon_bed_file = os.path.join(path, "amplicons.bed")
    primer_fasta_file = os.path.join(path, "oligos.fasta")

    with open(tsv_file, "w") as tsv, open(tsv_file_2, "w") as tsv2, open(amplicon_bed_file, "w") as bed, open(primer_fasta_file, "w") as fasta:
        print(
            "qpcr_scheme\toligo_type\tstart\tstop\tseq\tsize\tgc_best\ttemp_best\tmean_gc\tmean_temp\tpenalty\toff_target_amplicons",
            file=tsv2
        )
        print(
            "qpcr_scheme\toff_target_amplicons\tpenalty\tdeltaG\tlength\tstart\tstop\tseq",
            file=tsv
        )
        for n, amp in enumerate(final_schemes):
            amp_name = f"{scheme_name}_{n}"
            # write bed amplicon file
            print(
                f"{scheme_name}_ambiguous_consensus",
                amp["LEFT"][1],
                amp["RIGHT"][2],
                amp_name,
                round(amp["penalty"], 1),
                ".",
                sep="\t",
                file=bed
            )
            # write tsv
            amplicon_start = amp["LEFT"][1]
            amplicon_stop = amp["RIGHT"][2]
            if "off_targets" in amp:
                if amp["off_targets"]:
                    amplicon_has_off_target = "Yes"
                    write_BLAST_warning(amp_name, log_file)
                else:
                    amplicon_has_off_target = "No"
            else:
                amplicon_has_off_target = "n.d."
            amplicon_seq = ambiguous_consensus[amplicon_start:amplicon_stop]
            print(
                amp_name,
                amplicon_has_off_target,
                round(amp["penalty"], 1),
                amp["deltaG"],
                len(amplicon_seq),
                amplicon_start + 1,
                amplicon_stop,
                amplicon_seq,
                sep="\t",
                file=tsv
            )
            # write tsv2
            for oligo_type in ["LEFT", "PROBE", "RIGHT"]:
                seq = ambiguous_consensus[amp[oligo_type][1]:amp[oligo_type][2]]
                if oligo_type == "RIGHT" or (oligo_type == "PROBE" and amp["PROBE"][5] == "-"):
                    seq = primers.rev_complement(seq)
                    direction = "-"
                else:
                    direction = "+"

                permutations = get_permutations(seq)
                gc, temp = calc_mean_stats(permutations)
                primer_name = f"{amp_name}_{oligo_type}"

                print(
                    amp_name,
                    oligo_type,
                    amp[oligo_type][1] + 1,
                    amp[oligo_type][2],
                    seq.upper(),
                    len(seq),
                    round(primers.calc_gc(amp[oligo_type][0]), 1),
                    round(primers.calc_temp(amp[oligo_type][0]), 1),
                    gc,
                    temp,
                    round(amp[oligo_type][3], 1),
                    amplicon_has_off_target,
                    sep="\t",
                    file=tsv2
                )
                # write primer bed file
                write_primers_to_bed(
                    primer_bed_file,
                    scheme_name,
                    primer_name,
                    amp[oligo_type],
                    round(amp[oligo_type][3], 2),
                    direction,
                    seq.upper()
                )
                # write fasta
                print(f">{primer_name}\n{seq.upper()}", file=fasta)


def write_scheme_to_files(path, amplicon_scheme, ambiguous_consensus, scheme_name, mode, log_file):
    """
    write all relevant bed files and a tsv file with all primer stats
    """
    tsv_file = os.path.join(path, "primers.tsv")
    primer_bed_file = os.path.join(path, "primers.bed")
    amplicon_bed_file = os.path.join(path, "amplicons.bed")
    tabular_file = os.path.join(path, "primer_to_amplicon_assignment.tabular")

    # open files to write
    with open(tsv_file, "w") as tsv, open(amplicon_bed_file, "w") as bed, open(tabular_file, "w") as tabular:
        # write header for primer tsv
        print(
            "amlicon_name\tamplicon_length\tprimer_name\tprimer_name_all_primers\tpool\tstart\tstop\tseq\tsize\tgc_best\ttemp_best\tmean_gc\tmean_temp\tpenalty\toff_target_amplicons",
            file=tsv
        )
        amplicon_bed_records = []
        primer_bed_records = []
        primer_assignment_records = []
        pools = {amp.get("pool", 0) for amp in amplicon_scheme}
        for pool in pools:
            if mode == "single":
                primer_fasta_file = os.path.join(path, "primers.fasta")
            else:
                primer_fasta_file = os.path.join(path, f"primers_pool_{pool+1}.fasta")
            with open(primer_fasta_file, "w") as primer_fasta:
                for counter, amp in enumerate(amplicon_scheme[pool::len(pools)]):
                    # give a new amplicon name
                    amplicon_index = counter*len(pools) + pool
                    amp_name = f"{scheme_name}_{amplicon_index}"
                    # get left and right primers and their names
                    amp_length = amp["RIGHT"][2] - amp["LEFT"][1]
                    if "off_targets" in amp:
                        if amp["off_targets"]:
                            amplicon_has_off_target = "Yes"
                            write_BLAST_warning(amp_name, log_file)
                        else:
                            amplicon_has_off_target = "No"
                    else:
                        amplicon_has_off_target = "n.d."
                    # write amplicon bed
                    if mode == "tiled":
                        bed_score = pool+1
                    elif mode == "single":
                        bed_score = round(amp["LEFT"][3] + amp["RIGHT"][3], 1)
                    amplicon_bed_records.append(
                        (
                            amp["LEFT"][1],
                            amp["RIGHT"][2],
                            amp_name,
                            bed_score
                        )
                    )
                    primer_assignment_records.append(
                        (
                            # will need amplicon_index for sorting
                            amplicon_index,
                            (f"{amp_name}_LEFT", f"{amp_name}_RIGHT")
                        )
                    )
                    # write primer tsv and primer bed
                    for direction, primer in [("+", amp["LEFT"]), ("-", amp["RIGHT"])]:
                        seq = ambiguous_consensus[primer[1]:primer[2]]
                        if direction == "-":
                            seq = primers.rev_complement(seq)
                            primer_name = f"{amp_name}_RIGHT"
                        else:
                            primer_name = f"{amp_name}_LEFT"
                        # write primers to fasta pool file
                        print(f">{primer_name}\n{seq.upper()}", file=primer_fasta)
                        # calc primer parameters for all permutations
                        permutations = get_permutations(seq)
                        gc, temp = calc_mean_stats(permutations)
                        # write tsv file
                        print(
                            amp_name,
                            amp_length,
                            primer_name,
                            primer[-1],
                            pool+1,
                            primer[1] + 1,
                            primer[2],
                            seq.upper(),
                            len(primer[0]),
                            round(primers.calc_gc(primer[0]), 1),
                            round(primers.calc_temp(primer[0]), 1),
                            gc,
                            temp,
                            round(primer[3], 1),
                            amplicon_has_off_target,
                            sep="\t",
                            file=tsv
                        )
                        primer_bed_records.append(
                            (
                                # will need amplicon_index for sorting
                                amplicon_index,
                                (primer_name, primer, pool+1, direction, seq.upper())
                            )
                        )
        # write amplicon bed with amplicons sorted by start position
        for record in sorted(amplicon_bed_records, key=lambda x: x[0]):
            print(
                f"{scheme_name}_ambiguous_consensus",
                *record,
                ".",
                sep="\t",
                file=bed
            )
        # use sorting by amplicon index for primer assignment file
        for record in sorted(primer_assignment_records):
            print(
                *record[1],
                sep="\t",
                file=tabular
            )
        # same for primer bed
        for record in sorted(primer_bed_records):
            write_primers_to_bed(
                primer_bed_file,
                scheme_name,
                *record[1]
            )


def write_dimers(path, primer_dimers):
    """
    write dimers for which no replacement was found to file
    """
    tsv_file = os.path.join(path, "unsolvable_primer_dimers.tsv")
    with open(tsv_file, "w") as tsv:
        print(
            "pool\tprimer_name_1\tprimer_name_2\tdimer melting temp",
            file=tsv
        )
        for pool, primer1, primer2 in primer_dimers:
            print(
                pool+1,
                primer1[1],
                primer2[1],
                round(primers.calc_dimer(primer1[2][0], primer2[2][0]).tm, 1),
                sep="\t",
                file=tsv
            )

def entropy(chars, states):
    """
    input is a list of characters or numbers.
    calculate relative shannon's entropy. relative values are
    achieved by using the number of possible states as the base
    """
    ent = 0
    n_chars = len(chars)
    # only one char is in the list
    if n_chars <= 1:
        return ent
    # calculate the number of unique chars and their counts
    value, counts = np.unique(chars, return_counts=True)
    probs = counts/n_chars
    if np.count_nonzero(probs) <= 1:
        return ent

    for prob in probs:
        ent -= prob*math.log(prob, states)

    return ent

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


def entropy_subplot(ax, alignment_cleaned, scheme_name):
    """
    creates the entropy subplot
    """
    # - create entropy df
    entropy_df = alignment_entropy(alignment_cleaned)

    ax[0].fill_between(entropy_df["position"], entropy_df["entropy"], color="gainsboro", label="entropy")
    ax[0].plot(entropy_df["position"], entropy_df["average"], color="black", label="average entropy", linewidth=0.5)
    ax[0].set_ylim((0, 1))
    ax[0].set_xlim(0, max(entropy_df["position"]))
    ax[0].set_ylabel("normalized Shannon's entropy")
    ax[0].set_title(f"{scheme_name} amplicon design")
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)


def region_subplot(ax, primer_regions, location=0.95, color="darkorange", description="possible primer regions"):
    """
    creates the primer regions subplot
    """
    for region in primer_regions:
        ax[1].hlines(location, region[0], region[1], linewidth=5, color=color)
    # legend
    ax[1].hlines(location, primer_regions[0][1], primer_regions[0][1], label=description, linewidth=5, color=color)


def all_primer_subplot(ax, all_primers):
    """
    creates the all primer subplot
    """
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
        # legend
        ax[1].hlines(primer_position, all_primers[direction][primer][1], all_primers[direction][primer][2], linewidth=5, color=primer_color, label=primer_label)


def amplicon_subplot(ax, amplicon_scheme):
    """
    creates the amplicon subplot
    """
    for counter, amp in enumerate(amplicon_scheme):
        pool = amp.get("pool", 0)
        if pool == 0:
            position_amp = 0.7
            position_text = 0.6
        elif pool == 1:
            position_amp = 0.6
            position_text = 0.65
        left = amp["LEFT"]
        right = amp["RIGHT"]
        # amplicons
        ax[1].hlines(position_amp, left[1], right[2], linewidth=5)
        # text
        ax[1].text(right[2] - (right[2]-left[1])/2, position_text, str(counter), fontsize=8)
        # primers
        ax[1].hlines(position_amp, left[1], left[2], linewidth=5, color="red")
        ax[1].hlines(position_amp, right[1], right[2], linewidth=5, color="red")
    # legends
    ax[1].hlines(position_amp, left[1]+config.PRIMER_SIZES[1], right[2]-config.PRIMER_SIZES[1], linewidth=5, label="amplicons")
    ax[1].hlines(position_amp, left[1], left[2], linewidth=5, color="red", label="primers")


def qpcr_subplot(ax, amplicon_scheme):
    """
    creates the qpcr subplot
    """
    for counter, amp in enumerate(amplicon_scheme):
        left = amp["LEFT"]
        right = amp["RIGHT"]
        probe = amp["PROBE"]
        # amplicons
        ax[1].hlines(0.8, left[1], right[2], linewidth=5)
        # text
        ax[1].text(right[2] - (right[2]-left[1])/2, 0.65, str(counter), fontsize=8)
        # primers
        ax[1].hlines(0.8, left[1], left[2], linewidth=5, color="red")
        ax[1].hlines(0.8, right[1], right[2], linewidth=5, color="red")
        # probe
        ax[1].hlines(0.75, probe[1], probe[2], linewidth=5, color="darkgrey")

    # legends
    ax[1].hlines(0.8, left[1]+config.PRIMER_SIZES[1], right[2]-config.PRIMER_SIZES[1], linewidth=5, label="amplicons")
    ax[1].hlines(0.8, left[1], left[2], linewidth=5, color="red", label="primers")
    ax[1].hlines(0.75, probe[1], probe[2], linewidth=5, color="darkgrey", label="probe")


def varvamp_plot(path, alignment_cleaned, primer_regions, scheme_name, all_primers=None, amplicon_scheme=None, probe_regions=None):
    """
    creates overview plot for the amplicon design
    and per base coverage plots
    """
    # first plot: overview
    # create pdf name
    name = "amplicon_plot.pdf"
    out = os.path.join(path, name)
    # ini figure
    fig, ax = plt.subplots(2, 1, figsize=[22, 6], squeeze=True, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    fig.subplots_adjust(hspace=0)
    # entropy plot
    entropy_subplot(ax, alignment_cleaned, scheme_name)
    # primer regions plot
    region_subplot(ax, primer_regions)
    # probe region plot for probes
    if probe_regions is not None and amplicon_scheme is not None:
        region_subplot(ax, probe_regions, 0.9, color="dimgrey", description="possible probe regions")
        qpcr_subplot(ax, amplicon_scheme)
    # all primer plot
    elif all_primers is not None and amplicon_scheme is not None:
        all_primer_subplot(ax, all_primers)
        # amplicon, text and primer plot
        amplicon_subplot(ax, amplicon_scheme)
    # finalize
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].axes.get_yaxis().set_visible(False)
    ax[1].set_xlabel("alignment position")
    ax[1].set_ylim((0.5, 1))
    fig.legend(loc=(0.83, 0.7))

    # save fig
    fig.savefig(out, bbox_inches='tight')
    plt.close()


def get_primers_for_plot(amplicon_scheme, scheme_name, mode):
    """
    get the primers for per base pair plot (single, tiled)
    """
    amplicon_primers = []

    if mode == "SINGLE/TILED":
        oligo_types = ["LEFT", "RIGHT"]
    else:
        oligo_types = ["PROBE", "LEFT", "RIGHT"]

    for counter, amp in enumerate(amplicon_scheme):
        for oligo_type in oligo_types:
            primer_name = f"{scheme_name}_{counter}_{oligo_type}"
            amplicon_primers.append((primer_name, amp[oligo_type]))

    return amplicon_primers


def per_base_mismatch_plot(path, amplicon_scheme, threshold, scheme_name, mode="SINGLE/TILED"):
    """
    per base pair mismatch multiplot
    """
    out = os.path.join(path, "per_base_mismatches.pdf")

    amplicon_primers = get_primers_for_plot(amplicon_scheme, scheme_name, mode)
    # ini multi pdf
    with PdfPages(out) as pdf:
        # always print 4 primers to one page
        for i in range(0, len(amplicon_primers), 4):
            # ini figure
            primers_temp = amplicon_primers[i:i+4]
            fig, ax = plt.subplots(len(primers_temp), figsize=(12, len(primers_temp)*4), squeeze=True)
            # edge case if primer_temp has the length 1
            if len(primers_temp) == 1:
                ax = [ax]
            fig.suptitle("Per base mismatches", fontsize=18)
            fig.tight_layout(rect=[0.05, 0.05, 1, 0.98])
            fig.subplots_adjust(hspace=0.5)
            # plotting
            for idx, primer in enumerate(primers_temp):
                x = [pos+primer[1][1] for pos in range(0, len(primer[1][4]))]
                ax[idx].bar(x, primer[1][4], color='lightgrey', edgecolor='black')
                ax[idx].set_title(primer[0], loc="left")
                ax[idx].xaxis.set_ticks(np.arange(primer[1][1], primer[1][1]+len(x), 1))
                ax[idx].xaxis.set_ticklabels(x, rotation=45)
                ax[idx].set_ylabel(ylabel="freq of sequences")
                ax[idx].set_xlabel("position")
                # set yaxis lims reasonably
                if threshold < 1:
                    ax[idx].set_ylim(0, 1-threshold)
                else:
                    ax[idx].set_ylim(0, 0.001)
            # - to pdf
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()


def write_BLAST_warning(amplicon_name, log_file):
    """
    for each primer pair that has potential unspecific amplicons
    write warnings to file.
    """
    logging.raise_error(
        f"{amplicon_name} could produce off-targets. No better amplicon in this area was found.",
        log_file,
        exit=False,
    )
