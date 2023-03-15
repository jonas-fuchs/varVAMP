"""
amplicon search
"""

# BUILT-INS
import sys
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
        This method makes sure that the graph is symmetrical.
        """
        graph = {}
        for node in nodes:
            graph[node] = {}

        graph.update(init_graph)

        for node, edges in graph.items():
            for adjacent_node, value in edges.items():
                if graph[adjacent_node].get(node, False) == False:
                    graph[adjacent_node][node] = value

        return graph

    def get_nodes(self):
        """
        Returns the nodes of the graph.
        """
        return self.nodes

    def get_outgoing_edges(self, node):
        """
        Returns the neighbors of a node.
        """
        connections = []
        for out_node in self.nodes:
            if self.graph[node].get(out_node, False) != False:
                connections.append(out_node)
        return connections

    def value(self, node1, node2):
        """
        Returns the value of an edge between two nodes.
        """
        return self.graph[node1][node2]


def find_amplicons(left_primer_candidates, right_primer_candidates, opt_len, max_len):
    """
    finds all possible amplicons, creates a dictionary
    """
    amplicon_number = 0
    amplicon_dict = {}

    for left in left_primer_candidates:
        left_primer = left_primer_candidates[left]
        for right in right_primer_candidates:
            right_primer = right_primer_candidates[right]
            amplicon_length = right_primer[2] - left_primer[1]
            if opt_len <= amplicon_length <= max_len:
                if primers.calc_dimer(right_primer[0], left_primer[0]).tm <= config.MAX_DIMER_TMP:
                    # calculate length dependend amplicon costs as the cumulative primer
                    # score multiplied by the fold length of the optimal length.
                    amplicon_costs = (right_primer[3] + left_primer[3])*(amplicon_length/opt_len)
                    amplicon_name = "amplicon_"+str(amplicon_number)
                    amplicon_dict[amplicon_name] = [
                        left_primer[1],  # start
                        right_primer[2],  # stop
                        left,  # name left primer
                        right,  # name right primer
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
        amplicon = amplicons[current]
        start_overlap_pos = amplicon[0]+config.PRIMER_SIZES[2]
        stop_overlap_pos = amplicon[1] - min_overlap
        for next in amplicons:
            possible_next = amplicons[next]
            # check if the next amplicon lies within the start/stop range of
            # the current amplicon and if its non-overlapping part is large
            # enough to ensure space for a primer and the min overlap of the
            # following amplicon.
            if all((start_overlap_pos <= possible_next[0] <= stop_overlap_pos,
                    possible_next[1] > amplicon[1] + min_overlap
                    )):
                if current not in amplicon_graph:
                    amplicon_graph[current] = {next: possible_next[5]}
                else:
                    amplicon_graph[current][next] = possible_next[5]

    # return a graph object
    return Graph(nodes, amplicon_graph)


def dijkstra_algorithm(graph, start_node):
    """
    implementation of the dijkstra algorithm
    """

    previous_nodes = {}
    shortest_path = {vertex: float('infinity') for vertex in graph.get_nodes()}
    shortest_path[start_node] = 0

    pq = [(0, start_node)]

    while len(pq) > 0:
        current_distance, current_node = heapq.heappop(pq)

        if current_distance > shortest_path[current_node]:
            continue

        for neighbor in graph.get_outgoing_edges(current_node):
            distance = current_distance + graph.value(current_node, neighbor)

            # Only consider this new path if it's a better path
            if distance < shortest_path[neighbor]:
                shortest_path[neighbor] = distance
                previous_nodes[neighbor] = current_node
                heapq.heappush(pq, (distance, neighbor))

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


def find_best_covering_scheme(amplicons, amplicon_graph):
    """
    this brute forces the amplicon scheme search until the largest
    coverage with the minimal costs is achieved.
    """
    # ini
    coverage = 0
    best_coverage = 0
    max_stop = max(amplicons.items(), key=lambda x: x[1])[1][1]
    best_score = sys.maxsize

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
        final_amplicon_scheme = get_min_path(best_previous_nodes, best_shortest_path, best_start_node, best_target_node)
    else:
        # if no previous nodes are found but the single amplicon results in the largest
        # coverage - return as the best scheme
        final_amplicon_scheme = [best_start_node]

    return best_coverage, final_amplicon_scheme
