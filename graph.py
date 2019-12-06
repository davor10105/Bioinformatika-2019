from readers import PAFReader
from alignment import Overlap


class Edge():
    def __init__(self,node_from,node_to,type,strand,overlap_score,extension_score):
        self.node_from=node_from
        self.node_to=node_to
        self.type=type
        self.strand=strand
        self.overlap_score=overlap_score
        self.extension_score=extension_score

    def construct_edge(overlap,type):
        overlap_score=Overlap.get_overlap_score(overlap)

        if overlap.target_start==0:
            start_node=Node(overlap.query_name)
            end_node=Node(overlap.target_name)

            #TODO: ZASTO JE UOPCE RAZLIKA IZMEDU EXTENSION SCOREOVA?
            extension_score=Overlap.get_extension_score(overlap,overlap_score)
        else:
            start_node=Node(overlap.target_name)
            end_node=Node(overlap.query_name)

            extension_score=Overlap.get_extension_score(overlap,overlap_score)

        return Edge(start_node,end_node,type,overlap.relative_strand,overlap_score,extension_score)


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
        cr_cleaned_overlaps=Overlap.get_extension_overlaps(cr_paf.overlaps)

        rr_paf=PAFReader(read_read_overlap_path)
        rr_cleaned_overlaps=Overlap.get_extension_overlaps(rr_paf.overlaps)


        print(rr_cleaned_overlaps[0].target_length,rr_cleaned_overlaps[0].target_start,rr_cleaned_overlaps[0].target_end)
        print(rr_cleaned_overlaps[0].query_length,rr_cleaned_overlaps[0].query_start,rr_cleaned_overlaps[0].query_end)
        exit()



        '''with open(query_path,'r') as file:
            contig_lines=file.readlines()
            contigs={}
            for i in range(0,len(contig_lines),2):
                contigs[contig_lines[i][1:].strip()]=contig_lines[i+1]

            for overlap in paf.overlaps:
                sequence=contigs[overlap.query_name][overlap.query_start:overlap.query_end]
                node=Node(sequence)

                self.add_node(node)

        print(len(self.nodes))'''

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,node):
        self.nodes.remove(node)



if __name__=='__main__':
    overlap_graph=Graph('overlaps.paf','read_overlaps.paf')
