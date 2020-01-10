from readers import *
from alignment import Utils

class Edge():
    '''
    Defines an edge in an overlap graph, contains two Node objects, node_from,
    node_to, Overlap object, both scores (int) and a direction (str), which can be
    'to_right' or 'to_left'. The direction is considered during the path construction
    because only one type of direction can be used in a single path.
    '''
    def __init__(self,node_from,node_to,overlap,overlap_score,extension_score,direction):
        self.node_from=node_from
        self.node_to=node_to
        self.overlap=overlap
        self.overlap_score=overlap_score
        self.extension_score=extension_score
        self.direction=direction

    def construct_edges(overlap):
        '''
        From a given overlap, calculates both scores and constructs two edges
        for an overlap graph. For example, if the overlap is between nodes A - B,
        one contructed edge will be Edge(A,B,...,'to_right') and the other
        Edge(B,A,...,'to_left').
        '''
        #TODO: STO S TYPEOM? I REALTIVE STRANDOM?
        overlap_score=Utils.get_overlap_score(overlap)
        extension_score_left_to_right,extension_score_right_to_left=Utils.get_extension_scores(overlap,overlap_score)

        return [
                Edge(Node(overlap.query_name),Node(overlap.target_name),overlap,overlap_score,extension_score_left_to_right,'to_right'),
                Edge(Node(overlap.target_name),Node(overlap.query_name),overlap,overlap_score,extension_score_right_to_left,'to_left')
                ]

class Node():
    '''
    Defines a single node in an overlap graph.
    '''
    def __init__(self,name):
        self.name=name

    def __str__(self):
        return "name: %s"%(self.name)
    def __hash__(self):
        return hash(self.name)
    def __eq__(self,other):
        return self.name==other.name

class Graph():
    '''
    Defines an overlap graph. Contains edges in a dictionary self.edges that looks
    like {'to_right': {node1: [LIST_OF_THIS_NODE'S_TORIGHT_EDGES],node2: [LIST_OF_THIS_NODE'S_TORIGHT_EDGES]}...,
        'to_left': {node1: [LIST_OF_THIS_NODE'S_TOLEFT_EDGES],node2: [LIST_OF_THIS_NODE'S_TOLEFT_EDGES]}}.

    Contains the full str partial genomes of contigs in self.anchors and reads in self.extensions.
    '''
    def __init__(self,contig_path,read_path,contig_read_overlap_path,read_read_overlap_path):
        self.edges={'to_right':{},'to_left':{}}

        print('Loading data...')

        if contig_path.endswith('fasta'):
            self.anchors=FASTAReader(contig_path).reads
        else:
            self.anchors=FASTQReader(contig_path).reads
        if read_path.endswith('fasta'):
            self.extensions=FASTAReader(read_path).reads
        else:
            self.extensions=FASTQReader(read_path).reads

        cr_paf=PAFReader(contig_read_overlap_path)
        #cr_cleaned_overlaps=Utils.get_extension_overlaps(cr_paf.overlaps)
        cr_cleaned_overlaps=cr_paf.overlaps

        rr_paf=PAFReader(read_read_overlap_path)
        #rr_cleaned_overlaps=Utils.get_extension_overlaps(rr_paf.overlaps)
        rr_cleaned_overlaps=rr_paf.overlaps

        print('Data loaded.')

        all_overlaps=cr_cleaned_overlaps+rr_cleaned_overlaps
        for overlap in all_overlaps:
            new_edges=Edge.construct_edges(overlap)
            if new_edges[0].node_from not in self.edges['to_right']:
                self.edges['to_right'][new_edges[0].node_from ]=[]
            self.edges['to_right'][new_edges[0].node_from].append(new_edges[0])
            if new_edges[1].node_from not in self.edges['to_left']:
                self.edges['to_left'][new_edges[1].node_from]=[]
            self.edges['to_left'][new_edges[1].node_from].append(new_edges[1])

    def get_genome(self,node):
        '''
        Returns the string genome of a given node.
        '''
        if node.name in self.anchors:
            return self.anchors[node.name]
        if node.name in self.extensions:
            return self.extensions[node.name]
        print(node_name)
        print(self.anchors.keys())
        raise ValueError('Could not find the node name')

    def get_edges(self,node,direction):
        '''
        Returns the available edges from the given node in a given direction.
        If direction is None, returns edges from both directions (useful for the first
        node in a search when the direction is not yet decided).
        '''
        if direction is not None:
            try:
                return self.edges[direction][node]
            except:
                return []
        else:
            try:
                to_right_edges=self.edges['to_right'][node]
            except:
                to_right_edges=[]
            try:
                to_left_edges=self.edges['to_left'][node]
            except:
                to_left_edges=[]

            return to_right_edges+to_left_edges

    def reconstruct_path(self,end_state):
        '''
        From a given search State, returns the full reconstructed genome using backtracking.
        '''
        genome=''
        current_state=end_state
        count=0
        length=0
        all_length=0
        #print(current_state.node.name)
        #print(current_state.edge_from.overlap)
        if current_state.direction=='to_right':
            #print('PRAVA STRANA')
            node_genome=self.get_genome(current_state.node)
            genome=node_genome[current_state.edge_from.overlap.target_start:] #last node is from overlap start to the end of its get_genome
            switch_strand=False
            while current_state.previous_state.edge_from!=None:
                node_genome=self.get_genome(current_state.previous_state.node)
                all_length+=len(node_genome)
                node_genome=node_genome[current_state.previous_state.edge_from.overlap.target_end:current_state.edge_from.overlap.query_end]

                if current_state.edge_from.overlap.relative_strand=='-':
                    switch_strand=not switch_strand
                if switch_strand:
                    genome=Utils.reverse_complement(node_genome)+genome
                    #genome=node_genome+genome
                else:
                    genome=node_genome+genome

                current_state=current_state.previous_state
                count+=1
                length+=len(node_genome)
            node_genome=self.get_genome(current_state.previous_state.node)
            all_length+=len(node_genome)
            node_genome=node_genome[:current_state.edge_from.overlap.query_start]

            if current_state.edge_from.overlap.relative_strand=='-':
                switch_strand=not switch_strand
            if switch_strand:
                genome=Utils.reverse_complement(node_genome)+genome
                #genome=node_genome+genome
            else:
                genome=node_genome+genome

        else:
            #print('OBRATNO')
            node_genome=self.get_genome(current_state.node)
            genome=node_genome[:current_state.edge_from.overlap.query_start] #last node is from overlap start to the end of its get_genome
            switch_strand=False
            while current_state.previous_state.edge_from!=None:
                node_genome=self.get_genome(current_state.previous_state.node)
                all_length+=len(node_genome)
                node_genome=node_genome[current_state.edge_from.overlap.target_start:current_state.previous_state.edge_from.overlap.query_start]

                if current_state.edge_from.overlap.relative_strand=='-':
                    switch_strand=not switch_strand
                if switch_strand:
                    genome+=Utils.reverse_complement(node_genome)
                    #genome+=node_genome
                else:
                    genome+=node_genome

                current_state=current_state.previous_state
                count+=1
                length+=len(node_genome)
            node_genome=self.get_genome(current_state.previous_state.node)
            all_length+=len(node_genome)
            node_genome=node_genome[current_state.edge_from.overlap.target_start:]

            if current_state.edge_from.overlap.relative_strand=='-':
                switch_strand=not switch_strand
            if switch_strand:
                genome+=Utils.reverse_complement(node_genome)
                #genome+=node_genome
            else:
                genome+=node_genome
        '''
        print("broj nodeova")
        print(count)
        print('avg genome length')
        print(length/count)
        print('avg all length')
        print(all_length/count)
        print('genome length')
        print(len(genome))
        '''
        return genome

    def get_extension_genome(self,end_state):
        '''
        From a given search State, returns the full reconstructed genome using backtracking.
        '''
        genome=''
        current_state=end_state
        if current_state.node.name in self.anchors:
            return genome

        #print(current_state.node.name)
        #print(current_state.edge_from.overlap)
        if current_state.direction=='to_right':
            #print('PRAVA STRANA')
            node_genome=self.get_genome(current_state.node)
            genome=node_genome[current_state.edge_from.overlap.target_start:]
            while current_state.previous_state.node.name not in self.anchors:
                node_genome=self.get_genome(current_state.previous_state.node)
                node_genome=node_genome[current_state.previous_state.edge_from.overlap.target_end:current_state.edge_from.overlap.query_end]
                genome=node_genome+genome

                current_state=current_state.previous_state

        else:
            #print('OBRATNO')
            node_genome=self.get_genome(current_state.node)
            genome=node_genome[:current_state.edge_from.overlap.query_start]
            while current_state.previous_state.node.name not in self.anchors:
                node_genome=self.get_genome(current_state.previous_state.node)
                node_genome=node_genome[current_state.edge_from.overlap.target_start:current_state.previous_state.edge_from.overlap.query_start]
                genome+=node_genome

                current_state=current_state.previous_state
        return genome


class State():
    '''
    Defines a search State. Contains a node, state from which this state was
    found in self.previous_state, current path's score in self.score, direction of
    the edges, used_nodes in the search path (used to quickly check if all contigs
    were visited in this path) and edge_from, which is the edge that leads from
    previous_state to this state.
    '''
    def __init__(self,node,previous_state,score,direction,used_nodes,edge_from,anchors_found):
        self.node=node
        self.previous_state=previous_state
        self.score=score
        self.direction=direction
        self.used_nodes=used_nodes
        self.edge_from=edge_from
        self.anchors_found=anchors_found
    def __hash__(self):
        '''
        Nisam siguran treba li ovdje samo hash node.name-a za pamćenje visited stanja
        ili jos nesto, zasad ovak radi.
        '''
        #TODO: OVO JE DODANO
        return hash(self.node.name)

        #if self.previous_state==None:
        #    return hash((self.node.name,None,self.score,self.direction))
        #else:
        #    return hash((self.node.name,self.previous_state.node.name,self.score,self.direction))
    def __eq__(self,other):
        if other is None:
            return False
        if self.node==other.node and self.get_previous_node()==other.get_previous_node() and self.score==other.score and self.direction==other.direction:
            return True
        return False
    def __str__(self):
        return str(self.node)+' '+str(self.get_previous_node())
    def get_previous_node(self):
        '''
        Returns the previous state.
        '''
        if self.previous_state is None:
            return None
        return self.previous_state.node

    def compare(self,other):
        '''
        compares this State with another and returns the better state
        '''
        if not isinstance(other,State):
            return self
        if self.anchors_found>other.anchors_found:
            return self
        if self.score>other.score:
            return self
        return other


if __name__=='__main__':
    overlap_graph=Graph('overlaps_mjau.paf','read_overlaps.paf')
