from readers import PAFReader
from alignment import Utils

class Edge():
    def __init__(self,node_from,node_to,strand,overlap_score,extension_score):
        self.node_from=node_from
        self.node_to=node_to
        self.strand=strand
        self.overlap_score=overlap_score
        self.extension_score=extension_score

    def construct_edges(overlap):
        #TODO: STO S TYPEOM? I REALTIVE STRANDOM?
        overlap_score=Utils.get_overlap_score(overlap)
        extension_score_left_to_right,extension_score_right_to_left=Utils.get_extension_scores(overlap,overlap_score)

        return [
                Edge(Node(overlap.query_name),Node(overlap.target_name),overlap.relative_strand,overlap_score,extension_score_left_to_right),
                Edge(Node(overlap.target_name),Node(overlap.query_name),overlap.relative_strand,overlap_score,extension_score_right_to_left)
                ]

class Node():
    def __init__(self,name):
        self.name=name

    def __str__(self):
        return "name: %s"%(self.name)

class Graph():
    def __init__(self,contig_read_overlap_path,read_read_overlap_path):
        self.anchors={}
        self.extensions={}
        self.edges={}

        cr_paf=PAFReader(contig_read_overlap_path)
        cr_cleaned_overlaps=Utils.get_extension_overlaps(cr_paf.overlaps)

        rr_paf=PAFReader(read_read_overlap_path)
        rr_cleaned_overlaps=Utils.get_extension_overlaps(rr_paf.overlaps)

        edges=[]
        all_overlaps=cr_cleaned_overlaps+rr_cleaned_overlaps
        for overlap in all_overlaps:
            edges+=Edge.construct_edges(overlap)

        print(rr_cleaned_overlaps[0].target_length,rr_cleaned_overlaps[0].target_start,rr_cleaned_overlaps[0].target_end)
        print(rr_cleaned_overlaps[0].query_length,rr_cleaned_overlaps[0].query_start,rr_cleaned_overlaps[0].query_end)

        print(edges[0].node_from,edges[0].node_to,edges[0].overlap_score,edges[0].extension_score)
        print(edges[1].node_from,edges[1].node_to,edges[1].overlap_score,edges[1].extension_score)
        exit()

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,node):
        self.nodes.remove(node)



if __name__=='__main__':
    overlap_graph=Graph('overlaps_mjau.paf','read_overlaps.paf')
