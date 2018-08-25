"""
Tests for Issue #6, provided by @rlorigro
"""
import pytest
import poagraph
import seqgraphalignment


def generate_poa_graph(sequences):
    """
    Initialize graph and align all sequences
    :param sequences: sequences to align
    :return: graph: the completed POA graph resulting from the given sequences
    """
    init_sequence = sequences[0]
    init_label = "0"

    graph = poagraph.POAGraph(init_sequence, init_label)

    for i in range(1, len(sequences)):
        sequence = sequences[i]
        alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph,
                                                        fastMethod=False,
                                                        globalAlign=True,
                                                        matchscore=1,
                                                        mismatchscore=-1,
                                                        gapscore=-2)

        graph.incorporateSeqAlignment(alignment, sequence, str(i))

    return graph


def sequences_and_test(sequences, test_sequence):
    graph = generate_poa_graph(sequences)
    alignment = seqgraphalignment.SeqGraphAlignment(test_sequence, graph,
                                                    fastMethod=False,
                                                    globalAlign=True,
                                                    matchscore=1,
                                                    mismatchscore=-1,
                                                    gapscore=-2)

    graph.incorporateSeqAlignment(alignment, test_sequence, "test")
    alignments = graph.generateAlignmentStrings()

    result = alignments[-2][1].replace("-", "")
    return graph, result


def test_order_of_alignment_case1():
    sequences = ['ATATTGTGTAAGGCACAATTAACA',
                 'ATATTGCAAGGCACAATTCAACA',
                 'ATATTGCAAGGCACACAACA',
                 'ATGTGCAAGAGCACATAACA']
    test_sequence = "ATATTGCAAGGCACACTAACA"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case2():
    sequences = ["TTA", "TGC"]
    test_sequence = "TTGC"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case3():
    sequences = ["TAGTGAAAGAGGAAAAGAA"]
    test_sequence = "GCCCAGAAATTCCAGACCAGC"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case4():
    sequences = ["CTACTTGGGAGGCTGAGGTGG", "CCACTTGAGTTGAGG", "CTACTTGGGAAGCTAGAGGTGG"]
    test_sequence = "CTACTTGGGAGGCTGAGGTGG"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence
