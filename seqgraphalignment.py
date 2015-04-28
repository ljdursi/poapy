#!/usr/bin/env python
import poagraph
import numpy
import collections

class SeqGraphAlignment(object):
    __matchscore=4
    __mismatchscore=-2
    __opengapscore=-12
    __extendgapscore=-6

    def __init__(self, sequence, graph, fastMethod=True, globalAlign=False, *args, **kwargs):
        self.sequence   = sequence
        self.graph      = graph
        self.stringidxs = None
        self.nodeidxs   = None
        self.globalAlign= globalAlign
        if fastMethod:
            matches = self.alignStringToGraphFast(*args, **kwargs)
        else:
            matches = self.alignStringToGraphClear(*args, **kwargs)
        self.stringidxs, self.nodeidxs = matches

    def alignmentStrings(self):
        return ( "".join([self.sequence[i] if i is not None else "-" for i in self.stringidxs]), 
                 "".join([self.graph.nodedict[j].base if j is not None else "-" for j in self.nodeidxs]) )

    def matchscore(self,c1,c2):
        if c1 == c2:
            return self.__matchscore
        else:
            return self.__mismatchscore

    def matchscoreVec(self,c,v):
        res = numpy.where(v == c, self.__matchscore, self.__mismatchscore)
        return res
   
    def alignStringToGraphClear(self):
        """Align string to graph, following same approach as smith waterman
        example"""
        if not type(self.sequence) == str:
            raise TypeError("Invalid Type")
        
        l1 = self.graph.nNodes
        l2 = len(self.sequence)

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()

        # gap score "matrices"; will mostly be opening, just use dict to
        # indicate indices where a gap has oppened
        insertCost = collections.defaultdict(lambda : self.__opengapscore)
        deleteCost = collections.defaultdict(lambda : self.__opengapscore)

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i,node in enumerate(ni()):
            pbase = node.base

            for j,sbase in enumerate(self.sequence):
                # add all candidates to a list, pick the best
                    
                candidates = [(scores[i+1,j] + insertCost[(i+1,j)], i+1, j, "INS")]
                for predIndex in self.prevIndices(node, nodeIDtoIndex):
                    candidates.append( (scores[predIndex+1,j] + self.matchscore(sbase,pbase), predIndex+1, j, "MATCH") )
                    candidates.append( (scores[predIndex+1,j+1] + deleteCost[(predIndex+1,j+1)], predIndex+1, j+1, "DEL") )

                scores[i+1,j+1], backGrphIdx[i+1,j+1], backStrIdx[i+1,j+1], movetype = max(candidates)
                if movetype == "INS":
                    insertCost[(i+1,j+1)] = self.__extendgapscore
                if movetype == "DEL":
                    deleteCost[(i+1,j+1)] = self.__extendgapscore

                if not self.globalAlign and scores[i+1,j+1] < 0:
                    scores[i+1,j+1] = 0.
                    backGrphIdx[i+1,j+1] = -1
                    backStrIdx[i+1,j+1] = -1

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)


    def alignStringToGraphFast(self):
        """Align string to graph - using numpy to vectorize across the string
        at each iteration."""
        if not type(self.sequence) == str:
            raise TypeError("Invalid Type")
        
        l1 = self.graph.nNodes
        l2 = len(self.sequence)

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()

        seqvec = numpy.array(list(self.sequence))

        # gap score matrices, giving the score for opening/extending a gap _from_
        # this location
        insertCost = numpy.zeros((l1+1,l2+1),dtype=numpy.int)+self.__opengapscore
        deleteCost = numpy.zeros((l1+1,l2+1),dtype=numpy.int)+self.__opengapscore

        inserted   = numpy.zeros((l2+1),dtype=numpy.bool)

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i,node in enumerate(ni()):
            gbase = node.base
            predecessors = self.prevIndices(node, nodeIDtoIndex)

            matchpoints = self.matchscoreVec(gbase, seqvec)

            # calculate all best deletions, matches in one go over all
            # predecessors.
            # First calculate for the first predecessor:
            deletescore = scores[predecessors[0]+1,1:] + deleteCost[predecessors[0]+1,1:]
            bestdelete  = numpy.zeros((l2),dtype=numpy.int)+predecessors[0]+1

            matchscore = scores[predecessors[0]+1,:-1] + matchpoints
            bestmatch  = numpy.zeros((l2),dtype=numpy.int)+predecessors[0]+1

            for predecessor in predecessors[1:]:
                newdeletescore = scores[predecessor+1,1:] + deleteCost[predecessor+1,1:] 
                bestdelete     = numpy.where(newdeletescore > deletescore, predecessor+1, bestdelete)
                deletescore    = numpy.maximum(newdeletescore, deletescore)

                newmatchscore  = scores[predecessor+1,:-1] + matchpoints
                bestmatch      = numpy.where(newmatchscore > matchscore, predecessor+1, bestmatch)
                matchscore     = numpy.maximum(newmatchscore, matchscore)

            # choose best options available of match, delete
            deleted       = deletescore > matchscore
            bestmatchdel  = numpy.where  (deleted, bestdelete, bestmatch)

            scores     [i+1,1:] = numpy.maximum(deletescore, matchscore)
            backGrphIdx[i+1,1:] = bestmatchdel
            backStrIdx [i+1,1:] = numpy.where(deleted, numpy.arange(1,l2+1), numpy.arange(0,l2))
            deleteCost [i+1,1:] = numpy.where(deleted, self.__extendgapscore, deleteCost[i+1,1:])

            inserted[:] = False
            insscores = scores[i+1,:] + insertCost[i+1,:]
            for j in xrange(l2):
                if insscores[j] > scores[i+1,j+1]:
                    scores[i+1,j+1] = insscores[j]
                    backStrIdx[i+1,j+1] = j
                    inserted[j+1] = True
                    insscores[j+1] = scores[i+1,j+1] + insertCost[i+1,j+1]

            insertCost[i+1,inserted] = self.__extendgapscore
            deleteCost[i+1,inserted] = self.__opengapscore
            backGrphIdx[i+1,inserted] = i+1

            # if we're doing local alignment, don't let bad global alignment
            # drag us negative
            if not self.globalAlign:
                backGrphIdx[i+1,:] = numpy.where(scores[i+1,:] > 0, backGrphIdx[i+1,:], -1)
                backStrIdx [i+1,:] = numpy.where(scores[i+1,:] > 0, backStrIdx[i+1,:],  -1)
                scores[i+1,:]      = numpy.maximum(scores[i+1,:], 0)

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)

    def prevIndices(self, node, nodeIDtoIndex):
        """Return a list of the previous dynamic programming table indices 
           corresponding to predecessors of the current node."""
        prev = []
        for predID in node.inEdges.keys():
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

        nodeIDtoIndex= {}
        nodeIndexToID= {-1:None}
        # generate a dict of (nodeID) -> (index into nodelist (and thus matrix))
        ni = self.graph.nodeiterator()
        for (index,node) in enumerate(ni()):
            nodeIDtoIndex[node.ID] = index
            nodeIndexToID[index] = node.ID

        # Dynamic Programming data structures; scores matrix and backtracking
        # matrix
        scores = numpy.zeros((l1+1,l2+1))

        # initialize score
        if self.globalAlign:
            scores[0,:] = numpy.arange(0,l2+1)*self.__extendgapscore
            ni = self.graph.nodeiterator()
            for (index, node) in enumerate(ni()):
                prevIdxs = self.prevIndices(node, nodeIDtoIndex)
                best = 0 if len(prevIdxs) == 0 else scores[prevIdxs[0]+1,0]
                for prevIdx in prevIdxs:
                    best = max(best, scores[prevIdx+1,0])
                scores[index+1,0] = best + self.__extendgapscore

        # backtracking matrices
        backStrIdx = numpy.zeros((l1+1,l2+1),dtype=numpy.int)
        backGrphIdx = numpy.zeros((l1+1,l2+1),dtype=numpy.int)

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

        matches = []
        strindexes = []
        while (self.globalAlign or scores[besti,bestj] > 0) and not( besti == 0 and bestj == 0 ):
            nexti, nextj = backGrphIdx[besti,bestj], backStrIdx[besti,bestj]
            curstridx, curnodeidx = bestj-1, nodeIndexToID[besti-1]

            strindexes.insert(0, curstridx  if nextj != bestj else None)
            matches.insert   (0, curnodeidx if nexti != besti else None)

            besti, bestj = nexti, nextj

        return strindexes, matches

