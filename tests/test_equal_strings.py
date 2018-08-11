"""
Test for Issue #3, reported by @rlorigro
"""
import pytest
import poagraph
import seqgraphalignment

boolean = [True, False]
nseqs = [2, 3, 4]


@pytest.mark.parametrize('fast', boolean)
@pytest.mark.parametrize('glob', boolean)
@pytest.mark.parametrize('n_sequences', nseqs)
def test_equal_strings(fast, glob, n_sequences):
    sequence = "TATACCGGCG"
    sequences = [sequence]*n_sequences

    graph = poagraph.POAGraph(sequences[0], "0")
    for i in range(1, len(sequences)):
        alignment = seqgraphalignment.SeqGraphAlignment(sequences[i], graph,
                                                        fastMethod=fast,
                                                        globalAlign=glob,
                                                        matchscore=1,
                                                        mismatchscore=-1,
                                                        gapscore=-2)

        graph.incorporateSeqAlignment(alignment, sequence, str(i))

    alignments = graph.generateAlignmentStrings()
    matches = [alignments[0][1] == alignstr for _, alignstr in alignments]
    assert all(matches)
