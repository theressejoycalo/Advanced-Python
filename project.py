"""
author: Theresse Joy Calo
studentnumber: 2158606
"""
import sys
import networkx as nx
import pandas as pd
import json
import matplotlib.pyplot as plt


data = sys.argv[1]
k_value = int(data[-5])
x = int(data[4])

#   Subtask 1
#   This should be done with the function read_csv(name) that returns a dataframe with the data.


def read_csv(name: str):
    df = pd.read_csv(name, names=['segment', 'position', 'a', 'c', 'g', 't'])
    return df


dataframe = read_csv(data)

# Subtask 2
# This should be done with the function clean_data(df) where df is the dataframe with
# the data. This function returns the dataframe after cleaning all its errors.


def clean_data(df):
    df = _missing_position(df)
    df = _duplicated_position(df)
    df = _wrong_position(df)
    df = _duplicated_segments(df)
    return df


# Missing position in a segment. If the maximum position of some segment is m, but not all
# the records of the positions from 1 ... m exist. In this case, you have to ignore the whole
# segment.


def _missing_position(df):
    segment_list = list(set(df.segment))
    for i in segment_list:
        r_index = []
        position = []
        for index, row in df.iterrows():
            if row['segment'] == i:
                position.append(row['position'])
                r_index.append(index)
        m = max(position)
        n = set([item for item in range(1, m+1)])
        if set(position) != n:
            df = df.drop(r_index)
    return df


# Duplicated position in a segment. If there are multiple records of the same position with
# the same A, C, G, and T values, then remove the duplicates and keep only one copy of the
# record. However, if not all of the records agree on these values, then it is impossible to
# know which of them is the correct one. Hence, in this case, ignore the whole segment.


def _duplicated_position(df):
    df = df.drop_duplicates(subset=['segment', 'position', 'a', 'c', 'g', 't'], keep='first')
    segment_list = list(set(df.segment))
    for i in segment_list:
        positions = []
        r_index = []
        for index, row in df[df.segment == i].iterrows():
            positions.append(row.position)
            r_index.append(index)
        if len(positions) != len(set(positions)):
            df = df.drop(r_index)
    return df


# Wrong position. If at some position all the values of A, C, G, and T are 0s or more than one
# of them have the value 1. Clearly, this is a wrong record. In this case, you have to ignore
# the whole segment.

def _wrong_position(df):
    segment_list = list(set(df.segment))
    for i in segment_list:
        position = []
        r_index = []
        for index, row in df[df.segment == i].iterrows():
            position.append(row.position)
            r_index.append(index)
        for j in position:
            total = 0
            for index, row in df[(df.segment == i) & (df.position == j)].iterrows():
                total += row.a
                total += row.c
                total += row.g
                total += row.t
                if total != 1:
                    df = df.drop(r_index)
    return df

# Duplicated segments. If the sequences of two segments are identical (even if the segments
# have different numbers), then you need to remove one of them.


def _duplicated_segments(df):
    segment_list = list(set(df.segment))
    unique_segments = []
    for i in segment_list:
        position = list(df.position[df.segment == i])
        r_index = []
        segment = []
        for j in position:
            post = []
            for index, row in df[(df.segment == i) & (df.position == j)].iterrows():
                post.append(row.position)
                post.append(row.a)
                post.append(row.c)
                post.append(row.g)
                post.append(row.t)
                r_index.append(index)
                segment.append(post)
        if segment not in unique_segments:
            unique_segments.append(segment)
        else:
            df = df.drop(r_index)
    return df


dataframe_clean = clean_data(dataframe)

# Subtask 3
# Generate JSON sequences from the dataframe.
# This should be done with the function generate_sequences(df) that given the cleaned
# dataframe it returns a JSON with the sequences of the segments.


def generate_sequence(df):
    segment_list = list(set(df.segment))
    sequence_list = []
    for i in segment_list:
        sequence = ''
        position = []
        for index, row in df[df.segment == i].iterrows():
            position.append(row.position)
        position = list(df.position[df.segment == i])
        for j in position:
            for index, row in df[(df.segment == i) & (df.position == j)].iterrows():
                if row.a == 1:
                    sequence += 'a'
                elif row.c == 1:
                    sequence += 'c'
                elif row.g == 1:
                    sequence += 'g'
                elif row.t == 1:
                    sequence += 't'
        sequence_list.append(sequence)
    json_data = json.dumps(sequence_list)
    return json_data


json_object = generate_sequence(dataframe_clean)

# Subtask 4
# Construct de Bruijn graph.
# This should be done with the function construct_graph(json_data, k) where json_data
# is the one constructed from the previous step, and k is the length of k-mers. You should
# get the k from the name of the input csv file as mentioned before.


def construct_graph(json_data: json, k: int):
    g = nx.MultiDiGraph()
    k_mer = []
    json_read = json.loads(json_data)
    for item in json_read:
        length = len(item)
        i = [u for u in range(0, (length - k)+1)]
        for ii in i:
            k_mer.append(item[ii:ii+k])
    for record in k_mer:
        left_k_mer = record[0:k-1]
        right_k_mer = record[k-(k-1):k]
        if left_k_mer not in g:
            g.add_node(left_k_mer)
        if right_k_mer not in g:
            g.add_node(right_k_mer)
        g.add_edge(left_k_mer, right_k_mer)
    return g


bruijn_graph = construct_graph(json_object, k_value)

# Subtask 5
# Plot the de Bruijn graph.
# This should be done with function plot_graph(graph, filename) where graph is the
# one constructed in the previous step, while filename is the name of the output file. You
# should get the name of the file from the name of the input file as mentioned before.
# filename = DNA_[x].png


def plot_graph(graph, filename: str) -> None:
    nx.draw_networkx(graph, with_labels=True)
    plt.savefig(filename)
    return None


plot_graph(bruijn_graph, f'DNA_{x}')

# Subtask 6
# Check whether the de Bruijn graph can be sequenced.
# This should be done with the function is_valid_graph(graph) that returns a boolean
# value representing whether it satisfies the conditions of the Euler's path existence or not.


def is_valid_graph(graph) -> bool:
    """ To check whether the existence of a Euler's path conditions are met"""
    if (_condition_1(graph) and _condition_2(graph)) is True:
        return True
    else:
        return False


# Condition 1
# The graph must be connected.
# In case the graph is not Eulerian, but it satisfies
# the conditions of existence of an Euler's path, then we can make the graph Eulerian by adding
# an edge from the node that has extra incoming edge to the one that has the extra outgoing edge.

def _condition_1(graph) -> bool:
    """ Check that graph is connected."""
    nodes = list(graph.nodes)
    visited = []
    new_graph = None
    for i in nodes:
        for j in graph[i]:
            if j not in visited:
                visited.append(j)
    if set(nodes) == set(visited):
        return True
    elif nodes != list(set(visited)):
        new_graph = graph
        difference = set(nodes) - set(visited)
        new_graph.add_edge(nodes[-1], list(difference)[0])
    new_nodes = list(new_graph.nodes)
    new_visited = []
    for i in new_nodes:
        for j in new_graph[i]:
            if j not in new_visited:
                new_visited.append(j)
    if set(new_graph.nodes) == (set(new_visited)):
        return True
    else:
        return False


# Condition 2
# The number of nodes in the graph that have its InDegree different from its OutDegree is
# either 0 or 2.
# In case there are none of these nodes, then the whole graph is an Euler cycle and the graph
# is called Eulerian. Thus, the Euler's path in this case begins and ends at the same node.


def _condition_2(graph) -> bool:
    """ Check the Degree differences of the graph"""
    in_degree = set(graph.in_degree)
    out_degree = set(graph.in_degree)
    difference = len(in_degree - out_degree)
    if len(in_degree) == 0 and len(out_degree) == 0:
        return True
    elif difference == 0:
        return True
    elif difference == 2:
        if _condition_3(graph) is True:
            return True
        else:
            return False
    else:
        return False


# Condition 3
# Moreover, in case the number of these nodes is two, then one of them must have its
# InDegree - OutDegree = 1, while the other must have its OutDegree - InDegree = 1.
# Clearly, the one that has one extra outgoing edge will be the first edge to be traversed,
# while the one that the extra incoming edge will be the last edge.


def _condition_3(graph) -> bool:
    """ If nodes degree difference is equal to 2."""
    in_degree = set(graph.in_degree)
    out_degree = set(graph.in_degree)
    nodes_difference = list(in_degree - out_degree)
    node_1 = nodes_difference[0][0]
    node_2 = nodes_difference[1][0]
    if (graph.in_degree(node_1) - graph.out_degree(node_1) == 1) and \
            (graph.out_degree(node_2) - graph.in_degree(node_2) == 1):
        return True
    elif (graph.in_degree(node_2) - graph.out_degree(node_2) == 1) and \
            (graph.out_degree(node_1) - graph.in_degree(node_1) == 1):
        return True
    else:
        return False


# Subtask 7
# Construct DNA sequence.
# This should be done with the function construct_dna_sequence(graph) that returns a
# string representing the DNA sequence. Of course, this step should only be applied to valid graphs.
# In other words, the function should only be applied to a graph if the previous
# function returns True on that graph.
# Moreover, this function should print to the console the Euler's path found to get the DNA sequence.


def construct_dna_sequence(graph) -> str:
    dna_sequence = ''
    temp_sequence = []
    e_path_sequence = []
    index = len(list(graph.nodes())[0])
    if is_valid_graph(graph) is True:
        for (u, v) in graph.edges():
            if dna_sequence == '':
                dna_sequence += u
                dna_sequence += v[-1]
                e_path_sequence.append((u, v))
            else:
                if dna_sequence[-index:] == u:
                    dna_sequence += v[-1]
                    e_path_sequence.append((u, v))
                else:
                    temp_sequence.append((u, v))
        for (u, v) in temp_sequence:
            if dna_sequence[-index:] == u:
                dna_sequence += v[-1]
                temp_sequence.remove((u, v))
                e_path_sequence.append((u, v))
            elif u == v:
                u_list = [uu for (uu, vv) in e_path_sequence]
                u_location = u_list.index(u)
                v_list = [vv for (uu, vv) in e_path_sequence]
                v_location = v_list.index(v)
                if u_location - v_location == 1:
                    e_path_sequence.insert(u_location, (u, v))
                    dna_sequence = dna_sequence[:index+u_location]+v[-1]+dna_sequence[index+u_location:]
        print(e_path_sequence)
        return dna_sequence
    else:
        error_message = "DNA sequence can not be constructed."
        print(error_message)
        return error_message


dna_sequence_string = construct_dna_sequence(bruijn_graph)

# Subtask 8
# Save DNA sequence or write the error message.
# This should be done with the function save_output(s, filename) where s is the string
# of the DNA sequence or the error message that the DNA sequence can not be constructed,
# while filename is the name of the output file. You should get the name of the file from
# the name of the input file as mentioned before
# A txt file with the name DNA_[x].txt that contains:
# (1) the DNA sequence of the given data in case the de Bruijn graph satisfies the conditions
# of having an Euler's path; or
# (2) the message "DNA sequence can not be constructed.", otherwise.


def save_output(s: str, filename: str) -> None:
    """Save DNA sequence into a txt file or write the error message."""
    with open(filename, 'w') as file:
        file.write(s)
    return None


save_output(dna_sequence_string, f'DNA_{x}.txt')
