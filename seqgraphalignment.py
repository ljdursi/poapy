#!/usr/bin/env python
from __future__ import print_function
try:
    from builtins import range
    from builtins import object
except ImportError:
    pass
import numpy
import collections


class SeqGraphAlignment(object):
    __matchscore = 4
    __mismatchscore = -2
    __opengapscore = -4
    __extendgapscore = -2

    def __init__(self, sequence, graph, fastMethod=True, globalAlign=False,
                 matchscore=__matchscore, mismatchscore=__mismatchscore,
                 opengapscore=__opengapscore, extendgapscore=__extendgapscore, *args, **kwargs):
        self.__mismatchscore = mismatchscore
        self.__matchscore = matchscore
        self.__opengapscore = opengapscore
        self.__extendgapscore = extendgapscore
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
            return self.__matchscore
        else:
            return self.__mismatchscore

    def matchscoreVec(self, c, v):
        res = numpy.where(v == c, self.__matchscore, self.__mismatchscore)
        return res

    def alignStringToGraphSimple(self):
        """Align string to graph, following same approach as smith waterman
        example"""
        if not type(self.sequence) == str:
            raise TypeError("Invalid Type")

        nodeIDtoIndex, nodeIndexToID, m_scores, ins_scores, del_scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            pbase = node.base

            for j, sbase in enumerate(self.sequence):
                # add all candidates to a list, pick the best

                candidates = [(scores[i+1, j] + insertCost[(i+1, j)], i+1, j, "INS")]
                for predIndex in self.prevIndices(node, nodeIDtoIndex):
                    candidates.append((scores[predIndex+1, j] + self.matchscore(sbase, pbase), predIndex+1, j, "MATCH"))
                    candidates.append((scores[predIndex+1, j+1] + deleteCost[(predIndex+1, j+1)], predIndex+1, j+1, "DEL"))

                scores[i+1, j+1], backGrphIdx[i+1, j+1], backStrIdx[i+1, j+1], movetype = max(candidates)
                if movetype == "INS":
                    insertCost[(i+1, j+1)] = self.__extendgapscore
                if movetype == "DEL":
                    deleteCost[(i+1, j+1)] = self.__extendgapscore

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

        nodeIDtoIndex, nodeIndexToID, m_scores, ins_scores, del_scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()
        scores = numpy.maximum(numpy.maximum(m_scores, ins_scores), del_scores)

        seqvec = numpy.array(list(self.sequence))

        # having the inner loop as a function improves performance
        # can use Cython, etc on this for significant further improvements
        def insertions(i, l2, match_scores, ins_scores, del_scores):
            matches = match_scores[i+1, 1:] + self.__opengapscore + self.__extendgapscore
            dels = del_scores[i+1, 1:] + self.__opengapscore + self.__extendgapscore
            best = numpy.where(matches > dels, matches, dels)
            for j in range(1, l2+1):
                ins = self.__extendgapscore + ins_scores[i, j-1]
                ins_scores[i+1, j] = max(ins, best[j-1])

        def deletions(p, l2, m_scores, ins_scores, del_scores):
            del_match = m_scores[p + 1, 1:] + self.__opengapscore + self.__extendgapscore
            del_ins = ins_scores[p + 1, 1:] + self.__opengapscore + self.__extendgapscore
            del_del = del_scores[p + 1, 1:] + self.__extendgapscore
            return numpy.maximum(numpy.maximum(del_match, del_ins), del_del)

        def matches(p, l2, m_scores, ins_scores, del_scores, matchpoints):
            match_match = m_scores[p + 1, :-1] + matchpoints
            match_ins = ins_scores[p + 1, :-1] + matchpoints
            match_del = del_scores[p + 1, :-1] + matchpoints
            return numpy.maximum(numpy.maximum(match_match, match_ins), match_del)

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            gbase = node.base
            predecessors = self.prevIndices(node, nodeIDtoIndex)

            # calculate all best deletions, matches in one go over all
            # predecessors.
            # First calculate for the first predecessor:
            deletescore = deletions(predecessors[0], l2, m_scores, ins_scores, del_scores)
            bestdelete = numpy.zeros((l2), dtype=numpy.int)+predecessors[0]+1

            matchpoints = self.matchscoreVec(gbase, seqvec)
            matchscore = matches(predecessors[0], l2, m_scores, ins_scores, del_scores, matchpoints)
            bestmatch = numpy.zeros((l2), dtype=numpy.int)+predecessors[0]+1

            # then, the remaining
            for predecessor in predecessors[1:]:
                newdeletescore = deletions(predecessor, l2, m_scores, ins_scores, del_scores)
                bestdelete     = numpy.where(newdeletescore > deletescore, predecessor+1, bestdelete)
                deletescore    = numpy.maximum(newdeletescore, deletescore)

                gbase = self.graph.nodeIdxToBase(predecessor)
                matchpoints = self.matchscoreVec(gbase, seqvec)
                newmatchscore = matches(predecessor, l2, m_scores, ins_scores, del_scores, matchpoints)
                bestmatch     = numpy.where(newmatchscore > matchscore, predecessor+1, bestmatch)
                matchscore    = numpy.maximum(newmatchscore, matchscore)

            del_scores[i+1, 1:] = deletescore
            m_scores[i+1, 1:] = matchscore

            # choose best options available of match, delete
            deleted       = deletescore > matchscore
            backGrphIdx[i+1, 1:] = numpy.where(deleted, bestdelete, bestmatch)
            backStrIdx [i+1, 1:] = numpy.where(deleted, numpy.arange(1, l2+1), numpy.arange(0, l2))

            # insertions: updated in place, don't depend on predecessors
            insertions(i, l2, m_scores, ins_scores, del_scores)
            best = numpy.maximum(deletescore, matchscore)
            inserted = ins_scores[i+1, 1:] > best

            backGrphIdx[i+1, 1:] = numpy.where(inserted, numpy.arange(l2), backGrphIdx[i+1, 1:])
            for j, ins in enumerate(inserted):
                if ins:
                    backStrIdx[i+1, j+1] = j

            scores[i+1, 1:] = numpy.maximum(numpy.maximum(ins_scores[i+1, 1:], deletescore), matchscore)

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
        minint = numpy.iinfo(numpy.int16).min

        match_scores = numpy.zeros((l1+1, l2+1), dtype=numpy.int)
        match_scores[1:l1+1, 0] = minint
        match_scores[0, 1:l2+1] = minint

        del_scores = numpy.zeros((l1+1, l2+1), dtype=numpy.int)
        del_scores[1:l1+1, 0] = minint
        del_scores[0, 1:l2+1] = 0

        ins_scores = numpy.zeros((l1+1, l2+1), dtype=numpy.int)
        ins_scores[1:l1+1, 0] = 0
        ins_scores[0, 1:l2+1] = minint

        # initialize insertion score
        # if global align, penalty for starting at head != 0
        if self.globalAlign:
            ni = self.graph.nodeiterator()
            ins_scores[0, 0] = self.__opengapscore - self.__extendgapscore
            for (index, node) in enumerate(ni()):
                prevIdxs = self.prevIndices(node, nodeIDtoIndex)
                best = ins_scores[prevIdxs[0]+1, 0]
                for prevIdx in prevIdxs:
                    best = max(best, del_scores[prevIdx+1, 0])
                ins_scores[index+1, 0] = best + self.__extendgapscore
            ins_scores[0, 0] = 0
            del_scores[1:l1+1, 0] = self.__opengapscore + numpy.arange(l1)*self.__extendgapscore

        # backtracking matrices
        backStrIdx = numpy.zeros((l1+1, l2+1), dtype=numpy.int)
        backGrphIdx = numpy.zeros((l1+1, l2+1), dtype=numpy.int)

        return nodeIDtoIndex, nodeIndexToID, match_scores, ins_scores, del_scores, backStrIdx, backGrphIdx

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
