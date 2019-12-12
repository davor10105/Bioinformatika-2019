from readers import PAFReader
from alignment import Utils

class Edge():
    def __init__(self,node_from,node_to,strand,overlap_score,extension_score,direction):
        self.node_from=node_from
        self.node_to=node_to
        self.strand=strand
        self.overlap_score=overlap_score
        self.extension_score=extension_score
        self.direction=direction

    def construct_edges(overlap):
        #TODO: STO S TYPEOM? I REALTIVE STRANDOM?
        overlap_score=Utils.get_overlap_score(overlap)
        extension_score_left_to_right,extension_score_right_to_left=Utils.get_extension_scores(overlap,overlap_score)

        return [
                Edge(Node(overlap.query_name),Node(overlap.target_name),overlap.relative_strand,overlap_score,extension_score_left_to_right,'to_right'),
                Edge(Node(overlap.target_name),Node(overlap.query_name),overlap.relative_strand,overlap_score,extension_score_right_to_left,'to_left')
                ]

class Node():
    def __init__(self,name):
        self.name=name

    def __str__(self):
        return "name: %s"%(self.name)
    def __hash__(self):
        return hash(self.name)
    def __eq__(self,other):
        return self.name==other.name

class Graph():
    def __init__(self,contig_read_overlap_path,read_read_overlap_path):
        self.anchors={}
        self.extensions={}
        self.edges={'to_right':{},'to_left':{}}

        cr_paf=PAFReader(contig_read_overlap_path)
        cr_cleaned_overlaps=Utils.get_extension_overlaps(cr_paf.overlaps)

        rr_paf=PAFReader(read_read_overlap_path)
        rr_cleaned_overlaps=Utils.get_extension_overlaps(rr_paf.overlaps)

        all_overlaps=cr_cleaned_overlaps+rr_cleaned_overlaps
        for overlap in all_overlaps:
            new_edges=Edge.construct_edges(overlap)
            if new_edges[0].node_from not in self.edges['to_right']:
                self.edges['to_right'][new_edges[0].node_from ]=[]
            self.edges['to_right'][new_edges[0].node_from].append(new_edges[0])
            if new_edges[1].node_from not in self.edges['to_left']:
                self.edges['to_left'][new_edges[1].node_from]=[]
            self.edges['to_left'][new_edges[1].node_from].append(new_edges[1])

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,node):
        self.nodes.remove(node)

    def get_edges(self,node,direction):
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


class State():
    def __init__(self,node,previous_state,score,direction):
        self.node=node
        self.previous_state=previous_state
        self.score=score
        self.direction=direction
    def __hash__(self):
        return hash((self.node,self.previous_state,self.score,self.direction))
    def __eq__(self,other):
        if self.node==other.node and self.previous_state.node==other.previous_state.node and self.score==other.score and self.direction==other.direction:
            return True
        return False


if __name__=='__main__':
    overlap_graph=Graph('overlaps_mjau.paf','read_overlaps.paf')
