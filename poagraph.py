#!/usr/bin/env python
try:
    from builtins import zip
    from builtins import str
    from builtins import object
except ImportError:
    pass
import numpy
import seqgraphalignment
import textwrap
import collections

class Node(object):
    def __init__(self, nodeID=-1, base='N'):
        self.ID = nodeID
        self.base = base
        self.inEdges = {}
        self.outEdges = {}
        self.alignedTo = []

    def __str__(self):
        return "(%d:%s)" % (self.ID, self.base)

    def _add_edge(self, edgeset, neighbourID, label):
        if neighbourID is None:
            return
        # already present? just update labels
        if neighbourID in edgeset:
            edgeset[neighbourID].addLabel(label)
        else:
            edge = Edge(neighbourID, self.ID, label)
            edgeset[neighbourID] = edge

    def addInEdge(self, neighbourID, label):
        self._add_edge(self.inEdges, neighbourID, label)

    def addOutEdge(self, neighbourID, label):
        self._add_edge(self.outEdges, neighbourID, label)

    def nextNode(self, label):
        """Returns the first (presumably only) outward neighbour
           having the given edge label"""
        nextID = None
        for e in self.outEdges:
            if label in self.outEdges[e].labels:
                nextID = e
        return nextID

    @property
    def inDegree(self):
        return len(self.inEdges)

    @property
    def outDegree(self):
        return len(self.outEdges)

    @property
    def labels(self):
        """Returns all the labels associated with an in-edge or an out edge."""
        labelset = set([])
        for e in list(self.inEdges.values()):
            labelset = labelset.union(e.labels)
        for e in list(self.outEdges.values()):
            labelset = labelset.union(e.labels)
        return list(labelset)


class Edge(object):
    def __init__(self, inNodeID=-1, outNodeID=-1, label=None):
        self.inNodeID  = inNodeID
        self.outNodeID = outNodeID
        if label is None:
            self.labels = []
        elif type(label) == list:
            self.labels = label
        else:
            self.labels = [label]

    def addLabel(self, newlabel):
        if not newlabel in self.labels:
            self.labels.append(newlabel)
        return

    def __str__(self):
        nodestr = "(%d) -> (%d) " % (self.inNodeID, self.outNodeID)
        if self.labels is None:
            return nodestr
        else:
            return nodestr + self.labels.__str__()


class POAGraph(object):
    def addUnmatchedSeq(self, seq, label=None, updateSequences=True):
        """Add a completely independant (sub)string to the graph, 
           and return node index to initial and final node"""
        if seq is None:
            return

        firstID, lastID = None, None
        neededSort = self.needsSort

        for base in seq:
            nodeID = self.addNode(base)
            if firstID is None:
                firstID = nodeID
            if lastID is not None:
                self.addEdge(lastID, nodeID, label)
            lastID = nodeID

        self.__needsort = neededSort # no new order problems introduced
        if updateSequences:
            self.__seqs.append(seq)
            self.__labels.append(label)
            self.__starts.append(firstID)
        return firstID, lastID

    def __init__(self, seq=None, label=None):
        self._nextnodeID = 0
        self._nnodes = 0
        self._nedges = 0
        self.nodedict = {}
        self.nodeidlist = []   # allows a (partial) order to be imposed on the nodes
        self.__needsort = False
        self.__labels = []
        self.__seqs = []
        self.__starts = []

        if seq is not None:
            self.addUnmatchedSeq(seq, label)

    def addNode(self, base):
        nid = self._nextnodeID
        newnode = Node(nid, base)
        self.nodedict[nid] = newnode
        self.nodeidlist.append(nid)
        self._nnodes += 1
        self._nextnodeID += 1
        self._needsSort = True
        return nid

    def addEdge(self, start, end, label):
        if start is None or end is None:
            return

        if not start in self.nodedict:
            raise KeyError('addEdge: Start node not in graph: '+str(start))
        if not end in self.nodedict:
            raise KeyError('addEdge: End node not in graph: '+str(end))
        
        oldNodeEdges = self.nodedict[start].outDegree + self.nodedict[end].inDegree

        self.nodedict[start].addOutEdge( end, label )
        self.nodedict[end].addInEdge( start, label )

        newNodeEdges = self.nodedict[start].outDegree + self.nodedict[end].inDegree

        if newNodeEdges != oldNodeEdges:
            self._nedges += 1

        self._needsSort = True
        return

    @property
    def needsSort(self):
        return self.__needsort

    @property
    def nNodes(self):
        return self._nnodes

    @property
    def nEdges(self):
        return self._nedges

    def toposort(self):
        """Sorts node list so that all incoming edges come from nodes earlier in the list."""
        sortedlist = []
        completed = set()

        def dfs(start, completed, sortedlist):
            stack, started = [start], set()
            while stack:
                nodeID = stack.pop()

                if nodeID in completed:
                    continue

                if nodeID in started:
                    completed.add(nodeID)
                    sortedlist.insert(0, nodeID)
                    started.remove(nodeID)
                    continue

                successors = [s for s in self.nodedict[nodeID].outEdges.keys() if not s in completed]
                started.add(nodeID)
                stack.append(nodeID)
                stack.extend(successors)

        while len(sortedlist) < self.nNodes:
            found = None
            for node in self.nodedict:
                if not node in completed:
                    found = node
                    break
            assert found is not None
            dfs(found, completed, sortedlist)

        self.nodeidlist = sortedlist
        self._needsSort = False
        #self.testsort()
        return

    def testsort(self):
        """ Test the nodeidlist to make sure it is topologically sorted:
            eg, all predecessors of a node preceed the node in the list"""
        if self.nodeidlist is None:
            return
        seen_nodes = set()
        for nodeidx in self.nodeidlist:
            node = self.nodedict[nodeidx]
            for in_neighbour in node.inEdges:
                assert in_neighbour in seen_nodes
            seen_nodes.add(nodeidx)
        return

    def nodeiterator(self):
        if self.needsSort:
            self.toposort()

        def nodegenerator():
            for nodeidx in self.nodeidlist:
                yield self.nodedict[nodeidx]

        return nodegenerator

    def __str__(self):
        selfstr = ""
        ni = self.nodeiterator()
        for node in ni():
            selfstr += node.__str__() + "\n"
            for outIdx in node.outEdges:
                selfstr += "        " + node.outEdges[outIdx].__str__() + "\n"
        return selfstr
    
    def incorporateSeqAlignment(self, alignment, seq, label=None):
        """Incorporate a SeqGraphAlignemnt into the graph."""
        newseq     = alignment.sequence
        stringidxs = alignment.stringidxs
        nodeidxs   = alignment.nodeidxs

        firstID = None
        headID = None
        tailID = None

        # head, tail of sequence may be unaligned; just add those into the
        # graph directly
        validstringidxs = [si for si in stringidxs if si is not None]
        startSeqIdx,endSeqIdx = validstringidxs[0], validstringidxs[-1]
        if startSeqIdx > 0:
            firstID, headID = self.addUnmatchedSeq(newseq[0:startSeqIdx], label, updateSequences=False)
        if endSeqIdx < len(newseq):
            tailID, __ = self.addUnmatchedSeq(newseq[endSeqIdx+1:],  label, updateSequences=False)

        # now we march along the aligned part. For each base, we find or create
        # a node in the graph:
        #   - if unmatched, the corresponding node is a new node
        #   - if matched:
        #       - if matched to a node with the same base, the node is that node
        #       - if matched to a node with a different base whch is in turn
        #         aligned to a node with the same base, that aligned node is
        #         the node
        #       - otherwise, we create a new node.
        # In all cases, we create edges (or add labels) threading through the
        # nodes.
        for sindex, matchID in zip(stringidxs, nodeidxs):
            if sindex is None:
                continue
            base = newseq[sindex]
            if matchID is None:
                nodeID = self.addNode(base)
            elif self.nodedict[matchID].base == base:
                nodeID = matchID
            else:
                otherAligns = self.nodedict[matchID].alignedTo
                foundNode = None
                for otherNodeID in otherAligns:
                    if self.nodedict[otherNodeID].base == base:
                        foundNode = otherNodeID
                if foundNode is None:
                    nodeID = self.addNode(base)
                    self.nodedict[nodeID].alignedTo = [matchID] + otherAligns
                    for otherNodeID in [matchID] + otherAligns:
                        self.nodedict[otherNodeID].alignedTo.append(nodeID)
                else:
                    nodeID = foundNode

            self.addEdge(headID, nodeID, label)
            headID = nodeID
            if firstID is None:
                firstID = headID

        # finished the unaligned portion: now add an edge from the current headID to the tailID.
        self.addEdge(headID, tailID, label)

        # resort
        self.toposort()

        self.__seqs.append(seq)
        self.__labels.append(label)
        self.__starts.append(firstID)
        return

    def consensus(self, excludeLabels=None):
        if excludeLabels is None:
            excludeLabels = []

        if self.needsSort:
            self.toposort()

        nodesInReverse = self.nodeidlist[::-1]
        maxnodeID = max(nodesInReverse)+1
        nextInPath = [-1]*maxnodeID
        scores = numpy.zeros((maxnodeID))

        for nodeID in nodesInReverse:
            bestWeightScoreEdge = (-1, -1, None)

            for neighbourID in self.nodedict[nodeID].outEdges:
                e = self.nodedict[nodeID].outEdges[neighbourID]
                weight = len([ l for l in e.labels if l not in excludeLabels ])
                weightScoreEdge = (weight, scores[neighbourID], neighbourID)

                if weightScoreEdge > bestWeightScoreEdge:
                    bestWeightScoreEdge = weightScoreEdge
                
            scores[nodeID] = sum(bestWeightScoreEdge[0:2])
            nextInPath[nodeID] = bestWeightScoreEdge[2]

        pos = numpy.argmax(scores)
        path   = []
        bases  = []
        labels = []
        while pos is not None and pos > -1:
            path.append(pos)
            bases.append(self.nodedict[pos].base)
            labels.append(self.nodedict[pos].labels)
            pos = nextInPath[pos]

        return path, bases, labels


    def allConsenses(self, maxfraction=0.5):
        allpaths = []
        allbases = []
        alllabels = []
        exclusions = []

        passno = 0
        lastlen = 1000
        maxpasses = 10

        while len(exclusions) < len(self.__labels) and lastlen >= 10 and passno < maxpasses:
            path, bases, labellists = self.consensus(exclusions)
            if len(path) > 0:
                allpaths.append(path)
                allbases.append(bases)
                alllabels.append(labellists)

                labelcounts = collections.defaultdict(int)
                for ll in labellists:
                    for l in ll:
                        labelcounts[l] += 1

                for label, seq in zip(self.__labels, self.__seqs):
                    if label in labelcounts:
                        if labelcounts[label] >= maxfraction*len(seq):
                            exclusions.append(label)

            lastlen = len(path)
            passno += 1

        return list(zip(allpaths, allbases, alllabels))

    def generateAlignmentStrings(self):
        """ Return a list of strings corresponding to the alignments in the graph """

        # Step 1: assign node IDs to columns in the output
        #  column_index[node.ID] is the position in the toposorted node list
        #    of the node itself, or the earliest node it is aligned to.
        column_index = {}
        current_column = 0

        ni = self.nodeiterator()
        for node in ni():
            other_columns = [column_index[other] for other in node.alignedTo
                                                  if other in column_index]
            if len(other_columns) > 0:
                found_idx = min(other_columns)
            else:
                found_idx = current_column
                current_column += 1

            column_index[node.ID] = found_idx

        ncolumns = current_column

        # Step 2: given the column indexes, populate the strings
        #   corresponding to the sequences inserted in the graph
        seqnames = []
        alignstrings = []
        for label, start in zip(self.__labels, self.__starts):
            seqnames.append(label)
            curnode_id = start
            charlist = ['-']*ncolumns
            while curnode_id is not None:
                node = self.nodedict[curnode_id]
                charlist[column_index[curnode_id]] = node.base
                curnode_id = node.nextNode(label)
            alignstrings.append("".join(charlist))

        # Step 3: Same as step 2, but with consensus sequences
        consenses = self.allConsenses()
        for i, consensus in enumerate(consenses):
            seqnames.append('Consensus'+str(i))
            charlist = ['-']*ncolumns
            for path, base in zip(consensus[0], consensus[1]):
                charlist[column_index[path]] = base
            alignstrings.append("".join(charlist))

        return list(zip(seqnames, alignstrings))


    def jsOutput(self):
        """returns a list of strings containing a a description of the graph for viz.js, http://visjs.org"""

        # get the consensus sequence, which we'll use as the "spine" of the
        # graph
        path, __, __ = self.consensus()
        pathdict = {}
        for i,nodeID in enumerate(path):
            pathdict[nodeID] = i*150

        lines = [ 'var nodes = [']

        ni = self.nodeiterator()
        count = 0
        for node in ni():
            line = '    {id:'+str(node.ID)+', label: "'+node.base+'"'
            if node.ID in pathdict and count % 5 == 0:
                line += ', allowedToMoveX: false, x: ' + str(pathdict[node.ID]) + ', y: 0 , allowedToMoveY: true },'
            else:
                line += '},'
            lines.append( line )

        lines[-1] = lines[-1][:-1]
        lines.append( '];' )

        lines.append( ' ' )

        lines.append( 'var edges = [' )
        nlabels = len(self.__labels)
        ni = self.nodeiterator()
        for node in ni():
            nodeID = str(node.ID)
            for edge in node.outEdges:
                target = str(edge)
                weight = str(len(node.outEdges[edge].labels)+1)
                lines.append('    {from: '+nodeID+', to: '+target+', value: '+weight+'},')
            for alignededge in node.alignedTo:
                # These edges indicate alignment to different bases, and are
                # undirected; thus make sure we only plot them once:
                if node.ID > alignededge:
                    continue
                target = str(alignededge)
                lines.append('    {from: '+nodeID+', to: '+target+', value: 1, style: "dash-line"},')
        lines[-1] = lines[-1][:-1]
        lines.append( '];' )
        return lines

    def htmlOutput(self, outfile):
        header="""
                  <!doctype html>
                  <html>
                  <head>
                    <title>POA Graph Alignment</title>

                    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/vis/3.11.0/vis.min.js"></script>
                  </head>

                  <body>

                  <div id="mynetwork"></div>

                  <script type="text/javascript">
                    // create a network
                  """
        outfile.write(textwrap.dedent(header[1:]))
        lines = self.jsOutput() 
        for line in lines:
            outfile.write(line+'\n')
        footer="""
                  var container = document.getElementById('mynetwork');
                  var data= {
                    nodes: nodes,
                    edges: edges,
                  };
                  var options = {
                    width: '100%',
                    height: '800px'
                  };
                  var network = new vis.Network(container, data, options);
                </script>

                </body>
                </html>
                """
        outfile.write(textwrap.dedent(footer))
