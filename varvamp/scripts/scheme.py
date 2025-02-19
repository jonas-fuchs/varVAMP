"""
amplicon search
"""

# BUILT-INS
import heapq
import math

# varVAMP
from varvamp.scripts import config, primers


def construct_graph(nodes, init_graph):
    """
    This method makes sure that the graph is symmetrical, but sets the costs
    for nodes in the reverse direction to infinity to make sure dijkstra
    never goes to an amplicon that is in the wrong direction.
    """
    graph = {}
    for node in nodes:
        graph[node] = {}

    graph.update(init_graph)
    for node, neighbors in graph.items():
        for neighbor in neighbors.keys():
            if graph[neighbor].get(node, False) is False:
                graph[neighbor][node] = (float("infinity"), 0)
    return graph


class Graph(object):
    """
    a graph class
    """

    def __init__(self, nodes, init_graph):
        self.nodes = nodes
        self.graph = construct_graph(nodes, init_graph)

    def get_nodes(self):
        """
        Returns the nodes of the graph.
        """
        return self.nodes

    def get_neighbors(self, node):
        """
        Returns the neighbors of a node.
        """
        neighbors = []
        for out_node in self.nodes:
            if self.graph[node].get(out_node, False) is not False:
                neighbors.append(out_node)
        return neighbors

    def value(self, node1, node2):
        """
        Returns the value of an edge between two nodes.
        """
        return self.graph[node1][node2]

def find_amplicons(all_primers, opt_len, max_len):
    """
    finds all possible amplicons, creates a dictionary
    """
    amplicon_number = 0
    amplicons = []

    for left_name in all_primers["+"]:
        left_primer = all_primers["+"][left_name]
        for right_name in all_primers["-"]:
            right_primer = all_primers["-"][right_name]
            amplicon_length = right_primer[2] - left_primer[1]
            if not opt_len <= amplicon_length <= max_len:
                continue
            if primers.calc_dimer(left_primer[0], right_primer[0]).tm > config.PRIMER_MAX_DIMER_TMP:
                continue
            # calculate length dependend amplicon costs as the cumulative primer
            # penalty multiplied by the e^(fold length of the optimal length).
            amplicon_costs = (right_primer[3] + left_primer[3])*math.exp(amplicon_length/opt_len)
            amplicon_name = "AMPLICON_"+str(amplicon_number)
            amplicons.append(
                {
                    "id": amplicon_name,
                    "penalty": amplicon_costs,
                    "length": amplicon_length,
                    "LEFT": left_primer + [left_name],
                    "RIGHT": right_primer + [right_name],
                }
            )
            amplicon_number += 1
    return amplicons


def has_qualifying_overlap(current_amplicon, next_amplicon, min_overlap):
    """
    check if two amplicons overlap sufficiently to connect them in the graph
    """
    # connect amplicons if they sufficiently overlap because:
    # ... the start of next amplicon lies in the second half of the prior amplicon
    if next_amplicon["LEFT"][1] < current_amplicon["LEFT"][1] + current_amplicon["length"] / 2:
        return False
    # ... the stop of the left primer of the next amplicon does not lie in the minimum amplicon insert
    if next_amplicon["LEFT"][2] > current_amplicon["RIGHT"][1] - min_overlap:
        return False
    # ... half of the next amplicon does not overlap with the previous amplicon --> enough space for a
    # further amplicon that lies in the second half next amplicon and cannot overlap with a primer of the
    # current amplicon
    if next_amplicon["RIGHT"][2] <= current_amplicon["RIGHT"][2] + next_amplicon["length"] / 2:
        return False

    return True

def create_amplicon_graph(amplicons, min_overlap):
    """
    creates the amplicon graph.
    """
    # ini graph and vertices
    amplicon_graph = {}
    nodes = []

    for current_amplicon in amplicons:
        # remember all vertices
        amplicon_id = current_amplicon["id"]
        nodes.append(amplicon_id)
        for next_amplicon in amplicons:
            if not has_qualifying_overlap(current_amplicon, next_amplicon, min_overlap):
                continue
            # --> write to graph
            if amplicon_id not in amplicon_graph:
                amplicon_graph[amplicon_id] = {
                    next_amplicon["id"]: (
                        next_amplicon.get("off_targets", False),
                        next_amplicon["penalty"]
                    )
                }
            else:
                amplicon_graph[amplicon_id][next_amplicon["id"]] = (
                    next_amplicon.get("off_targets", False),
                    next_amplicon["penalty"]
                )

    # return a graph object
    return Graph(nodes, amplicon_graph)


def dijkstra_algorithm(graph, start_node):
    """
    implementation of the dijkstra algorithm
    """

    previous_nodes = {}
    shortest_path = {node: (float('infinity'), 0) for node in graph.get_nodes()}
    shortest_path[start_node] = (0, 0)

    nodes_to_test = [((0, 0), start_node)]

    while nodes_to_test:
        current_distance, current_node = heapq.heappop(nodes_to_test)
        if current_distance > shortest_path[current_node]:
            continue
        for neighbor in graph.get_neighbors(current_node):
            off_targets, base_penalty = graph.value(current_node, neighbor)
            distance = (
                current_distance[0] + off_targets,
                current_distance[1] + base_penalty
            )
            # Only consider this new path if it's a better path
            if not distance < shortest_path[neighbor]:
                continue
            shortest_path[neighbor] = distance
            previous_nodes[neighbor] = current_node
            heapq.heappush(nodes_to_test, (distance, neighbor))

    return previous_nodes, shortest_path


def get_end_node(previous_nodes, shortest_path, amplicons):
    """
    get the target node with the lowest penalty costs out of all
    nodes that have the same maximum end position
    """
    stop_nucleotide, possible_end_nodes = 0, {}

    for node in previous_nodes.keys():
        # check if a node has a larger stop -> empty dict and set new
        # best stop nucleotide
        amplicon_stop = amplicons[node]["RIGHT"][2]
        if amplicon_stop > stop_nucleotide:
            possible_end_nodes = {node: shortest_path[node]}
            stop_nucleotide = amplicon_stop
        # if nodes have the same stop nucleotide, add to dictionary
        elif amplicon_stop == stop_nucleotide:
            possible_end_nodes[node] = shortest_path[node]

    # return the end node with the lowest penalty costs
    return min(possible_end_nodes.items(), key=lambda x: x[1])


def get_min_path(previous_nodes, start_node, target_node):
    """
    get the min path from the start to stop node from the
    previosuly calculated shortest path
    """
    path = []
    node = target_node

    while node != start_node:
        path.append(node)
        node = previous_nodes[node]
    # Add the start node manually
    path.append(start_node)

    # return the inverse list
    return path[::-1]


def create_scheme(amplicon_path, amplicons_by_id):
    """
    creates the final tiled-amplicon scheme
    """
    amplicon_scheme = []

    for pool in (0, 1):
        for amp_id in amplicon_path[pool::2]:
            amplicons_by_id[amp_id]["pool"] = pool
            amplicon_scheme.append(amplicons_by_id[amp_id])

    return amplicon_scheme


def find_best_covering_scheme(amplicons, amplicon_graph):
    """
    this brute forces the amplicon scheme search until the largest
    coverage with the minimal costs is achieved.
    """
    # ini
    best_coverage = 0
    max_stop = max(amplicons, key=lambda x: x["RIGHT"][2])["RIGHT"][2]
    lowest_costs = (float('infinity'),)
    # a dict for fast access to amplicons by their ID
    amps_by_id = {amp["id"]: amp for amp in amplicons}
    for amplicon in amplicons:
        # if the currently best coverage + start nucleotide of the currently tested amplicon
        # is smaller than the maximal stop nucleotide there might be a better amplicon
        # scheme that covers more of the genome
        if amplicon["LEFT"][1] + best_coverage <= max_stop:
            previous_nodes, shortest_path = dijkstra_algorithm(amplicon_graph, amplicon["id"])
            # only continue if there are previous_nodes
            if previous_nodes:
                target_node, costs = get_end_node(previous_nodes, shortest_path, amps_by_id)
                coverage = amps_by_id[target_node]["RIGHT"][2] - amplicon["LEFT"][1]
                # if the new coverage is larger, go for the larger coverage
                if coverage > best_coverage:
                    best_start_node = amplicon["id"]
                    best_target_node = target_node
                    best_previous_nodes = previous_nodes
                    lowest_costs = costs
                    best_coverage = coverage
                # if the coverages are identical, go for the lowest costs
                elif coverage == best_coverage:
                    if costs < lowest_costs:
                        best_start_node = amplicon["id"]
                        best_target_node = target_node
                        best_previous_nodes = previous_nodes
                        lowest_costs = costs
                        best_coverage = coverage
            else:
                # check if the single amplicon has the largest coverage so far
                coverage = amplicon["length"]
                if coverage > best_coverage:
                    best_start_node = amplicon["id"]
                    best_previous_nodes = previous_nodes
                    lowest_costs = (
                        amplicon.get("off_targets", False), amplicon["penalty"])
                    best_coverage = coverage
        # no need to check more, the best covering amplicon scheme was found and
        # has the minimal costs compared to the schemes with the same coverage
        else:
            break

    if best_previous_nodes:
        amplicon_path = get_min_path(best_previous_nodes, best_start_node, best_target_node)
    else:
        # if no previous nodes are found but the single amplicon results in the largest
        # coverage - return as the best scheme
        amplicon_path = [best_start_node]
    return best_coverage, create_scheme(amplicon_path, amps_by_id)


def test_scheme_for_dimers(amplicon_scheme):
    """
    test the lowest-cost scheme for primer dimers
    """

    primer_dimers = []
    pools = {amp["pool"] for amp in amplicon_scheme}
    for pool in pools:
        # test the primer dimers only within the respective pools
        tested_primers = []
        for amp_index, amp in enumerate(amplicon_scheme):
            if amp["pool"] != pool:
                continue
            for primer in ["LEFT", "RIGHT"]:
                # remember where the currrent primer was in the scheme
                # store the amplicon's index in the scheme, the current primer's original name and its details
                current_primer = (amp_index, amp[primer][-1], amp[primer][:-1])
                current_seq = current_primer[2][0]
                for tested in tested_primers:
                    tested_seq = tested[2][0]
                    if primers.calc_dimer(current_seq, tested_seq).tm <= config.PRIMER_MAX_DIMER_TMP:
                        continue
                    primer_dimers.append((current_primer, tested))
                # and remember all tested primers
                tested_primers.append(current_primer)

    return primer_dimers


def get_overlapping_primers(dimer, left_primer_candidates, right_primer_candidates):
    """
    get overlapping primers of a primer dimer pair that have been previously
    excluded. returns list of list with possible primers for each primer
    in primer dimer.
    """

    overlapping_primers = []
    # test each primer in dimer
    for amp_index, primer_name, primer in dimer:
        overlapping_primers_temp = []
        thirds_len = int((primer[2] - primer[1]) / 3)
        # get the middle third of the primer (here are the previously excluded primers)
        overlap_set = set(range(primer[1] + thirds_len, primer[2] - thirds_len))
        # check in which list to look for them
        if "RIGHT" in primer_name:
            primers_to_test = right_primer_candidates
        else:
            primers_to_test = left_primer_candidates
        # and check this list for all primers that overlap
        for potential_new in primers_to_test:
            primer_positions = list(range(potential_new[1], potential_new[2]))
            if not any(x in primer_positions for x in overlap_set):
                continue
            overlapping_primers_temp.append((amp_index, primer_name, potential_new))

        overlapping_primers.append(overlapping_primers_temp)

    return overlapping_primers


def test_overlaps_for_dimers(overlapping_primers):
    """
    test the overlapping primers for dimers. return new primers.
    """
    for first_overlap in overlapping_primers[0]:
        for second_overlap in overlapping_primers[1]:
            # return the first match. primers are sorted by penalty.
            # first pair that makes it has the lowest penalty
            if primers.calc_dimer(first_overlap[2][0], second_overlap[2][0]).tm <= config.PRIMER_MAX_DIMER_TMP:
                return [first_overlap, second_overlap]


def check_and_solve_heterodimers(amplicon_scheme, left_primer_candidates, right_primer_candidates, all_primers):
    """
    check scheme for heterodimers, try to find
    new primers that overlap and replace the existing ones.
    this can lead to new primer dimers. therefore the scheme
    is checked a second time. if there are still primer dimers
    present the non-solvable dimers are returned
    """

    primer_dimers = test_scheme_for_dimers(amplicon_scheme)

    if primer_dimers:
        print(f"varVAMP found {len(primer_dimers)} dimer pairs in scheme ... trying to find replacements")
    else:
        return []

    for dimer in primer_dimers:
        # get overlapping primers that have not been considered
        overlapping_primers = get_overlapping_primers(dimer, left_primer_candidates, right_primer_candidates)
        # test all possible primers against each other for dimers
        new_primers = test_overlaps_for_dimers(overlapping_primers)
        # now change these primers in the scheme
        if new_primers:
            for amp_index, primer_name, primer in new_primers:
                # overwrite in final scheme
                # ATTENTION: doesn't update the amplicon penalty currently
                # This is ok only because that value isn't used after and doesn't get reported anywhere.
                if "LEFT" in primer_name:
                    strand = "+"
                    amplicon_scheme[amp_index]["LEFT"] = primer + [primer_name]
                else:
                    strand = "-"
                    amplicon_scheme[amp_index]["RIGHT"] = primer + [primer_name]
                amplicon_scheme[amp_index]["length"] = amplicon_scheme[amp_index]["RIGHT"][2] - amplicon_scheme[amp_index]["LEFT"][1]
                # and in all primers
                all_primers[strand][primer_name] = primer
    # get remaining dimers in the revised scheme and add pool identifier for reporting
    primer_dimers = [
        (amplicon_scheme[primer1[0]]["pool"], primer1, primer2)
        for primer1, primer2 in test_scheme_for_dimers(amplicon_scheme)
    ]

    return primer_dimers


def find_single_amplicons(amplicons, n):
    """
    find non-overlapping amplicons with low penalties
    from all found amplicons. only for the SINGLE mode.
    """
    # sort amplicons
    sorted_amplicons = sorted(amplicons, key=lambda x: (x.get("off_targets", False), x["penalty"]))
    to_retain = []
    retained_ranges = []
    # find lowest non-overlapping
    for amp in sorted_amplicons:
        overlaps_retained = False
        amp_range = range(amp["LEFT"][1], amp["RIGHT"][2])
        for r in retained_ranges:
            if amp_range.start < r.stop and r.start < amp_range.stop:
                overlaps_retained = True
                break
        if not overlaps_retained:
            retained_ranges.append(amp_range)
            to_retain.append(amp)
            if len(to_retain) == n:
                break

    return to_retain
