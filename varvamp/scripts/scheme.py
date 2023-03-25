"""
amplicon search
"""

# BUILT-INS
import heapq

# varVAMP
from varvamp.scripts import config, primers


class Graph(object):
    """
    a graph class
    """

    def __init__(self, nodes, init_graph):
        self.nodes = nodes
        self.graph = self.construct_graph(nodes, init_graph)

    def construct_graph(self, nodes, init_graph):
        """
        This method makes sure that the graph is symmetrical, but sets the score
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
                    graph[neighbor][node] = float("infinity")

        return graph

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
    amplicon_dict = {}

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
            # score multiplied by the fold length of the optimal length.
            amplicon_costs = (right_primer[3] + left_primer[3])*(amplicon_length/opt_len)
            amplicon_name = "amplicon_"+str(amplicon_number)
            amplicon_dict[amplicon_name] = [
                left_primer[1],  # start
                right_primer[2],  # stop
                left_name,  # name left primer
                right_name,  # name right primer
                amplicon_length,  # amplicon length
                amplicon_costs  # costs
            ]
            amplicon_number += 1

    return amplicon_dict


def create_amplicon_graph(amplicons, min_overlap):
    """
    creates the amplicon graph.
    """
    # ini graph and vertices
    amplicon_graph = {}
    nodes = []

    # add the maximum len of a primer to ensure that possible amplicon starts
    # before the min overlap
    min_overlap = min_overlap + config.PRIMER_SIZES[2]

    for current in amplicons:
        # remember all vertices
        nodes.append(current)
        current_amplicon = amplicons[current]
        start = current_amplicon[0] + current_amplicon[4]/2
        stop = current_amplicon[1] - min_overlap
        for next in amplicons:
            next_amplicon = amplicons[next]
            # check if the next amplicon lies within the start/stop range of
            # the current amplicon and if its non-overlapping part is large
            # enough to ensure space for a primer and the min overlap of the
            # following amplicon.
            if not all((start <= next_amplicon[0] <= stop, next_amplicon[1] > current_amplicon[1] + next_amplicon[4]/2)):
                continue
            if current not in amplicon_graph:
                amplicon_graph[current] = {next: next_amplicon[5]}
            else:
                amplicon_graph[current][next] = next_amplicon[5]

    # return a graph object
    return Graph(nodes, amplicon_graph)


def dijkstra_algorithm(graph, start_node):
    """
    implementation of the dijkstra algorithm
    """

    previous_nodes = {}
    shortest_path = {node: float('infinity') for node in graph.get_nodes()}
    shortest_path[start_node] = 0

    nodes_to_test = [(0, start_node)]

    while nodes_to_test:
        current_distance, current_node = heapq.heappop(nodes_to_test)
        if current_distance > shortest_path[current_node]:
            continue
        for neighbor in graph.get_neighbors(current_node):
            distance = current_distance + graph.value(current_node, neighbor)
            # Only consider this new path if it's a better path
            if not distance < shortest_path[neighbor]:
                continue
            shortest_path[neighbor] = distance
            previous_nodes[neighbor] = current_node
            heapq.heappush(nodes_to_test, (distance, neighbor))

    return previous_nodes, shortest_path


def get_end_node(previous_nodes, shortest_path, amplicons):
    """
    get the target node with the lowest score out of all
    nodes that have the same maximum end position
    """
    stop_nucleotide = 0

    for node in previous_nodes.keys():
        # check if an node has a larger stop -> empty dict and set new
        # best stop nucleotide
        if amplicons[node][1] > stop_nucleotide:
            possible_end_nodes = {}
            possible_end_nodes[node] = shortest_path[node]
            stop_nucleotide = amplicons[node][1]
        # if nodes have the same stop nucleotide, add to dictionary
        elif amplicons[node][1] == stop_nucleotide:
            possible_end_nodes[node] = shortest_path[node]

    # return the end node with the lowest score
    return min(possible_end_nodes.items(), key=lambda x: x[1])


def get_min_path(previous_nodes, shortest_path, start_node, target_node):
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


def create_scheme_dic(amplicon_scheme, amplicons, all_primers):
    """
    creates the final scheme dictionary
    """

    scheme_dictionary = {
        0: {},
        1: {}
    }

    for pool in (0, 1):
        for amp in amplicon_scheme[pool::2]:
            scheme_dictionary[pool][amp] = {}
            primers = [amplicons[amp][2], amplicons[amp][3]]
            scheme_dictionary[pool][amp][primers[0]] = all_primers["+"][primers[0]]
            scheme_dictionary[pool][amp][primers[1]] = all_primers["-"][primers[1]]

    return scheme_dictionary


def find_best_covering_scheme(amplicons, amplicon_graph, all_primers):
    """
    this brute forces the amplicon scheme search until the largest
    coverage with the minimal costs is achieved.
    """
    # ini
    coverage = 0
    best_coverage = 0
    max_stop = max(amplicons.items(), key=lambda x: x[1])[1][1]
    best_score = float('infinity')

    for start_node in amplicons:
        # if the currently best coverage + start nucleotide of the currently tested amplicon
        # is smaller than the maximal stop nucleotide there might be a better amplicon
        # scheme that covers more of the genome
        if amplicons[start_node][0] + best_coverage <= max_stop:
            previous_nodes, shortest_path = dijkstra_algorithm(amplicon_graph, start_node)
            # only continue if there are previous_nodes
            if previous_nodes:
                target_node, score = get_end_node(previous_nodes, shortest_path, amplicons)
                coverage = amplicons[target_node][1] - amplicons[start_node][0]
                # if the new coverage is larger, go for the larger coverage
                if coverage > best_coverage:
                    best_start_node = start_node
                    best_target_node = target_node
                    best_previous_nodes = previous_nodes
                    best_shortest_path = shortest_path
                    best_score = score
                    best_coverage = coverage
                # if the coverages are identical, go for the lowest costs
                elif coverage == best_coverage:
                    if score < best_score:
                        best_start_node = start_node
                        best_target_node = target_node
                        best_previous_nodes = previous_nodes
                        best_shortest_path = shortest_path
                        best_score = score
                        best_coverage = coverage
            else:
                # check if the single amplicon has the largest coverage so far
                coverage = amplicons[start_node][1] - amplicons[start_node][0]
                if coverage > best_coverage:
                    best_start_node = start_node
                    best_previous_nodes = previous_nodes
                    best_shortest_path = shortest_path
                    best_score = amplicons[start_node][5]
                    best_coverage = coverage
        # no need to check more, the best covering amplicon scheme was found and
        # has the minimal score compared to the schemes with the same coverage
        else:
            break

    if best_previous_nodes:
        amplicon_scheme = get_min_path(best_previous_nodes, best_shortest_path, best_start_node, best_target_node)
    else:
        # if no previous nodes are found but the single amplicon results in the largest
        # coverage - return as the best scheme
        amplicon_scheme = [best_start_node]

    return best_coverage, create_scheme_dic(amplicon_scheme, amplicons, all_primers)


def test_scheme_for_dimers(amplicon_scheme):
    """
    test the best scoring scheme for primer dimers
    """

    primer_dimers = []

    for pool in amplicon_scheme:
        # test the primer dimers only within the respective pools
        tested_primers = []
        for amp in amplicon_scheme[pool]:
            for primer in amplicon_scheme[pool][amp]:
                # remember where the currrent primer was in the scheme
                current_primer = (pool, amp, primer, amplicon_scheme[pool][amp][primer])
                current_seq = current_primer[3][0]
                for tested in tested_primers:
                    tested_seq = tested[3][0]
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
    for primer in dimer:
        overlapping_primers_temp = []
        # check in which list to look for them
        overlap_range = range(primer[3][1], primer[3][2]+1)
        overlap_set = set(overlap_range)
        if "RIGHT" in primer[2]:
            primers_to_test = right_primer_candidates
        else:
            primers_to_test = left_primer_candidates
        # and check this list for all primers that overlap
        for potential_new in primers_to_test:
            primer_positions = list(range(potential_new[1], potential_new[2]+1))
            if not any(x in primer_positions for x in overlap_set):
                continue
            overlapping_primers_temp.append((primer[0], primer[1], primer[2], potential_new))

        overlapping_primers.append(overlapping_primers_temp)

    return overlapping_primers


def test_overlaps_for_dimers(overlapping_primers):
    """
    test the overlapping primers for dimers. return new primers.
    """
    for first_overlap in overlapping_primers[0]:
        for second_overlap in overlapping_primers[1]:
            # return the first match. primers are sorted by score.
            # first pair that makes it has the lowest score
            if primers.calc_dimer(first_overlap[3][0], second_overlap[3][0]).tm <= config.PRIMER_MAX_DIMER_TMP:
                return [first_overlap, second_overlap]


def check_and_solve_heterodimers(amplicon_scheme, left_primer_candidates, right_primer_candidates, all_primers):
    """
    check scheme for heterodimers, try to find
    new primers that overlap and replace the existing ones.
    this can lead to new primer dimers. therefore the
    process is repeated until no primer dimers are found
    in the updated scheme or all found primer dimers have
    no replacements.
    """
    not_solvable = []

    primer_dimers = test_scheme_for_dimers(amplicon_scheme)

    while primer_dimers:
        for dimer in primer_dimers:
            # skip the primer dimers that have not been previously solved
            if dimer in not_solvable:
                continue
            overlapping_primers = get_overlapping_primers(dimer, left_primer_candidates, right_primer_candidates)
            # test all possible primers against each other for dimers
            new_primers = test_overlaps_for_dimers(overlapping_primers)
            # now change these primers in the scheme
            if new_primers:
                for new in new_primers:
                    # overwrite in final scheme
                    amplicon_scheme[new[0]][new[1]][new[2]] = new[3]
                    # and in all primers
                    if "LEFT" in new[2]:
                        strand = "+"
                    else:
                        strand = "-"
                    all_primers[strand][new[2]] = new[3]
            # or remember the dimers for which varvamp did not find a replacement.
            else:
                not_solvable.append(dimer)
        # none of the primer dimers of this iteration could be solved
        if all(x in not_solvable for x in primer_dimers):
            break
        # some could be solved. lets check the updated scheme again.
        else:
            primer_dimers = test_scheme_for_dimers(amplicon_scheme)

    return not_solvable
