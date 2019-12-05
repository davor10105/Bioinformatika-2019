from readers import PAFReader


class Edge():
    def __init__(self,node_from,node_to,type):
        self.node_from=node_from
        self.node_to=node_to
        self.type=type

class Node():
    def __init__(self,sequence,edges=[]):
        self.sequence=sequence
        self.edges=edges #edges is a list of Edge

    def __str__(self):
        return "sequence: %s, connected nodes: %s"%(self.sequence,len(self.connected_nodes))

class Graph():
    def __init__(self,overlap_path,query_path):
        self.nodes=[]

        paf=PAFReader(overlap_path)
        with open(query_path,'r') as file:
            contig_lines=file.readlines()
            contigs={}
            for i in range(0,len(contig_lines),2):
                contigs[contig_lines[i][1:].strip()]=contig_lines[i+1]

            for overlap in paf.overlaps:
                sequence=contigs[overlap.query_name][overlap.query_start:overlap.query_end]
                node=Node(sequence)

                self.add_node(node)

        print(len(self.nodes))

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,node):
        self.nodes.remove(node)


if __name__=='__main__':
    overlap_graph=Graph('overlaps_mjau.paf','data/ecoli/contigs.fasta')
