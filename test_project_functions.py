"""
author: Theresse Joy Calo
studentnumber: 2158606
"""

from pytest import mark
from project_functions import (
    clean_data,
    generate_sequence,
    construct_graph,
    is_valid_graph,
    construct_dna_sequence,
)
import networkx as nx
import pandas as pd
import json

# Subtask 2
# This should be done with the function clean_data(df) where df is the dataframe with
# the data. This function returns the dataframe after cleaning all its errors.


@mark.parametrize(
    "given, expected",
    [
        (
            {'segment': [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6],
             'position': [1, 2, 4, 5, 1, 2, 3, 4, 4, 1, 2, 3, 4, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
             'a': [0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
             'c': [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
             'g': [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
             't': [1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0]},
            {'segment': [2, 2, 2, 2],
             'position': [1, 2, 3, 4],
             'a': [0, 0, 0, 1],
             'c': [0, 0, 1, 0],
             'g': [0, 1, 0, 0],
             't': [1, 0, 0, 0]},
        ),
    ],
)
def test_clean_data(given: str, expected: dict) -> None:
    clean = (clean_data(pd.DataFrame(given))).reset_index(drop=True)
    assert clean.equals(pd.DataFrame(expected))


# Subtask 3
# Generate JSON sequences from the dataframe.
# This should be done with the function generate_sequences(df) that given the cleaned
# dataframe it returns a JSON with the sequences of the segments.

@mark.parametrize(
    "given, expected",
    [
        (
                {'segment': [2, 2, 2, 2],
                 'position': [1, 2, 3, 4],
                 'a': [0, 0, 0, 1],
                 'c': [0, 0, 1, 0],
                 'g': [0, 1, 0, 0],
                 't': [1, 0, 0, 0]},
                ['tgca'],

        ),
    ],
)
def test_generate_sequences(given, expected) -> None:
    assert generate_sequence(pd.DataFrame(given)) == json.dumps(expected)


# Subtask 4
# Construct de Bruijn graph.
# This should be done with the function construct_graph(json_data, k) where json_data
# is the one constructed from the previous step, and k is the length of k-mers. You should
# get the k from the name of the input csv file as mentioned before.


def test_construct_graph() -> None:
    k = 3
    given = ['tgca', 'gtta', 'cgta']
    expected = [('tg', 'gc'), ('gc', 'ca'), ('gt', 'tt'), ('gt', 'ta'), ('tt', 'ta'), ('cg', 'gt')]
    graph_given = construct_graph(json.dumps(given), k)
    graph_expected = nx.MultiDiGraph(expected)
    assert graph_given.nodes == graph_expected.nodes
    assert graph_given.edges == graph_expected.edges

# Subtask 6
# Check whether the de Bruijn graph can be sequenced.
# This should be done with the function is_valid_graph(graph) that returns a boolean
# value representing whether it satisfies the conditions of the Euler's path existence or not.


@mark.parametrize(
    "given, expected",
    [
        ([('tg', 'gc'), ('gc', 'ca')], True),
        ([('tg', 'gc'), ('gc', 'ca'), ('at', 'ga')], False)
    ],
)
def test_is_valid_graph(given, expected) -> None:
    assert is_valid_graph(nx.MultiDiGraph(given)) == expected


# Subtask 7
# Construct DNA sequence.
# This should be done with the function construct_dna_sequence(graph) that returns a
# string representing the DNA sequence. Of course, this step should only be applied to valid graphs.
# In other words, the function should only be applied to a graph if the previous
# function returns True on that graph.
# Moreover, this function should print to the console the Euler's path found to get the DNA sequence.


@mark.parametrize(
    "given, expected",
    [
        ([('tg', 'gc'), ('gc', 'ca'), ('ca', 'at'), ('at', 'gt')],
         'tgcatt'),
        ([('tg', 'gc'), ('gc', 'ca'), ('at', 'ga')], "DNA sequence can not be constructed.")
    ],
)
def test_construct_dna_sequence(given, expected) -> None:
    assert construct_dna_sequence(nx.MultiDiGraph(given)) == expected
