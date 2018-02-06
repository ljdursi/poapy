#!/usr/bin/env python
from __future__ import print_function
try:
    from builtins import range
    from builtins import object
except ImportError:
    pass
import numpy


class SeqGraphAlignment(object):
    __matchscore = 1
    __mismatchscore = -1
    __gap = -2

    def __init__(self, sequence, graph, fastMethod=True, globalAlign=False,
                 matchscore=__matchscore, mismatchscore=__mismatchscore,
                 gapscore=__gap, *args, **kwargs):
        self._mismatchscore = mismatchscore
        self._matchscore = matchscore
        self._gap = gapscore
        self.sequence    = sequence
        self.graph       = graph
        self.stringidxs  = None
        self.nodeidxs    = None
        self.globalAlign = globalAlign
        if fastMethod:
            matches = self.alignStringToGraphFast(*args, **kwargs)
        else:
            matches = self.alignStringToGraphSimple(*args, **kwargs)
        self.stringidxs, self.nodeidxs = matches

    def alignmentStrings(self):
        return ("".join([self.sequence[i] if i is not None else "-" for i in self.stringidxs]),
                "".join([self.graph.nodedict[j].base if j is not None else "-" for j in self.nodeidxs]))

    def matchscore(self, c1, c2):
        if c1 == c2:
            return self._matchscore
        else:
            return self._mismatchscore

    def matchscoreVec(self, c, v):
        res = numpy.where(v == c, self._matchscore, self._mismatchscore)
        return res

    def alignStringToGraphSimple(self):
        """Align string to graph, following same approach as smith waterman
        example"""
        if not type(self.sequence) == str:
            raise TypeError("Invalid Type")

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            pbase = node.base

            for j, sbase in enumerate(self.sequence):
                # add all candidates to a list, pick the best
                candidates = [(scores[i+1, j] + self._gap, i+1, j, "INS")]
                for predIndex in self.prevIndices(node, nodeIDtoIndex):
                    candidates += [(scores[predIndex+1, j+1] + self._gap, predIndex+1, j+1, "DEL")]
                    candidates += [(scores[predIndex+1, j] + self.matchscore(sbase, pbase), predIndex+1, j, "MATCH")]

                scores[i+1, j+1], backGrphIdx[i+1, j+1], backStrIdx[i+1, j+1], movetype = max(candidates)

                if not self.globalAlign and scores[i+1, j+1] < 0:
                    scores[i+1, j+1] = 0.
                    backGrphIdx[i+1, j+1] = -1
                    backStrIdx[i+1, j+1] = -1

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)

    def alignStringToGraphFast(self):
        """Align string to graph - using numpy to vectorize across the string
        at each iteration."""
        if not type(self.sequence) == str:
            raise TypeError("Invalid Type")

        l2 = len(self.sequence)
        seqvec = numpy.array(list(self.sequence))

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()
        inserted = numpy.zeros((l2), dtype=numpy.bool)

        # having the inner loop as a function improves performance
        # can use Cython, etc on this for significant further improvements
        # can't vectorize this since there's a loop-carried dependency
        #  along the string
        def insertions(i, l2, scores, inserted):
            inserted[:] = False
            for j in range(l2):
                insscore = scores[i+1, j] + self._gap
                if insscore >= scores[i+1, j+1]:
                    scores[i+1, j+1] = insscore
                    inserted[j] = True

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            gbase = node.base
            predecessors = self.prevIndices(node, nodeIDtoIndex)

            # calculate all best deletions, matches in one go over all
            # predecessors.

            # First calculate for the first predecessor, over all string posns:
            deletescore = scores[predecessors[0]+1, 1:] + self._gap
            bestdelete = numpy.zeros((l2), dtype=numpy.int)+predecessors[0]+1

            matchpoints = self.matchscoreVec(gbase, seqvec)
            matchscore = scores[predecessors[0]+1, 0:-1] + matchpoints
            bestmatch = numpy.zeros((l2), dtype=numpy.int)+predecessors[0]+1

            # then, the remaining
            for predecessor in predecessors[1:]:
                newdeletescore = scores[predecessor+1, 1:] + self._gap
                bestdelete     = numpy.where(newdeletescore > deletescore, predecessor+1, bestdelete)
                deletescore    = numpy.maximum(newdeletescore, deletescore)

                gbase = self.graph.nodeIdxToBase(predecessor)
                matchpoints = self.matchscoreVec(gbase, seqvec)
                newmatchscore = scores[predecessor+1, 0:-1] + matchpoints
                bestmatch     = numpy.where(newmatchscore > matchscore, predecessor+1, bestmatch)
                matchscore    = numpy.maximum(newmatchscore, matchscore)

            # choose best options available of match, delete
            deleted       = deletescore >= matchscore
            backGrphIdx[i+1, 1:] = numpy.where(deleted, bestdelete, bestmatch)
            backStrIdx [i+1, 1:] = numpy.where(deleted, numpy.arange(1, l2+1), numpy.arange(0, l2))
            scores[i+1, 1:] = numpy.where(deleted, deletescore, matchscore)

            # insertions: updated in place, don't depend on predecessors
            insertions(i, l2, scores, inserted)
            backGrphIdx[i+1, 1:] = numpy.where(inserted, i+1, backGrphIdx[i+1, 1:])
            backStrIdx[i+1, 1:] = numpy.where(inserted, numpy.arange(l2), backStrIdx[i+1, 1:])

            # if we're doing local alignment, don't let bad global alignment
            # drag us negative
            if not self.globalAlign:
                backGrphIdx[i+1, :] = numpy.where(scores[i+1, :] > 0, backGrphIdx[i+1, :], -1)
                backStrIdx [i+1, :] = numpy.where(scores[i+1, :] > 0, backStrIdx[i+1, :], -1)
                scores[i+1, :]      = numpy.maximum(scores[i+1, :], 0)

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)

    def prevIndices(self, node, nodeIDtoIndex):
        """Return a list of the previous dynamic programming table indices
           corresponding to predecessors of the current node."""
        prev = []
        for predID in list(node.inEdges.keys()):
            prev.append(nodeIDtoIndex[predID])
        # if no predecessors, point to just before the graph
        if len(prev) == 0:
            prev = [-1]
        return prev

    def initializeDynamicProgrammingData(self):
        """Initalize the dynamic programming tables:
            - set up scores array
            - set up backtracking array
            - create index to Node ID table and vice versa"""
        l1 = self.graph.nNodes
        l2 = len(self.sequence)

        nodeIDtoIndex = {}
        nodeIndexToID = {-1: None}
        # generate a dict of (nodeID) -> (index into nodelist (and thus matrix))
        ni = self.graph.nodeiterator()
        for (index, node) in enumerate(ni()):
            nodeIDtoIndex[node.ID] = index
            nodeIndexToID[index] = node.ID

        # Dynamic Programming data structures; scores matrix and backtracking
        # matrix
        scores = numpy.zeros((l1+1, l2+1), dtype=numpy.int)

        # initialize insertion score
        # if global align, penalty for starting at head != 0
        if self.globalAlign:
            scores[0, :] = numpy.arange(l2+1)*self._gap

            ni = self.graph.nodeiterator()
            for (index, node) in enumerate(ni()):
                prevIdxs = self.prevIndices(node, nodeIDtoIndex)
                best = scores[prevIdxs[0]+1, 0]
                for prevIdx in prevIdxs:
                    best = max(best, scores[prevIdx+1, 0])
                scores[index+1, 0] = best + self._gap

        # backtracking matrices
        backStrIdx = numpy.zeros((l1+1, l2+1), dtype=numpy.int)
        backGrphIdx = numpy.zeros((l1+1, l2+1), dtype=numpy.int)

        return nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx

    def backtrack(self, scores, backStrIdx, backGrphIdx, nodeIndexToID):
        """Backtrack through the scores and backtrack arrays.
           Return a list of sequence indices and node IDs (not indices, which
           depend on ordering)."""
        besti, bestj = scores.shape
        besti -= 1
        bestj -= 1
        if not self.globalAlign:
            besti, bestj = numpy.argwhere(scores == numpy.amax(scores))[-1]
        else:
            # still have to find best final index to start from
            terminalIndices = []
            ni = self.graph.nodeiterator()
            for (index, node) in enumerate(ni()):
                if node.outDegree == 0:
                    terminalIndices.append(index)
            besti = terminalIndices[0]
            bestscore = scores[besti, bestj]
            for i in terminalIndices[1:]:
                score = scores[i, bestj]
                if score > bestscore:
                    bestscore, besti = score, i

        matches = []
        strindexes = []
        while (self.globalAlign or scores[besti, bestj] > 0) and not(besti == 0 and bestj == 0):
            nexti, nextj = backGrphIdx[besti, bestj], backStrIdx[besti, bestj]
            curstridx, curnodeidx = bestj-1, nodeIndexToID[besti-1]

            strindexes.insert(0, curstridx if nextj != bestj else None)
            matches.insert   (0, curnodeidx if nexti != besti else None)

            besti, bestj = nexti, nextj

        return strindexes, matches
