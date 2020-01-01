from graph import *
from readers import *

class CostSearch():
    def __init__(self,graph):
        self.graph=graph

    def search(self,start_node,end_node):
        start_state=State(start_node,None,0,None,{start_node},None)

        visited=set()
        open=[start_state]
        current_state=start_state

        anchor_names=set(self.graph.anchors.keys())
        while len(anchor_names-set([node.name for node in current_state.used_nodes]))!=0:
            #print(current_state.node)
            #print(end_node)
            #print(len(open))
            if len(open)==0:
                raise ValueError('No path to end node')

            current_state=open.pop(0)
            if current_state in visited:
                continue
            visited.add(current_state)

            for edge in self.graph.get_edges(current_state.node,current_state.direction):
                if edge.node_to not in current_state.used_nodes:
                    new_used_nodes=set(current_state.used_nodes)
                    new_used_nodes.add(edge.node_to)
                    new_state=State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction,new_used_nodes,edge)
                    if new_state not in visited:
                        open.append(new_state)

            open=sorted(open,key=lambda state: state.score,reverse=True)

        #print(current_state.used_nodes)
        return current_state

if __name__=='__main__':
    print('Creating the overlap graph...')
    overlap_graph=Graph('./data/ecoli/contigs.fasta','./data/ecoli/reads.fasta','overlaps_mjau.paf','read_overlaps.paf')
    print('Overlap graph done.')
    print('Starting the search...')
    search=CostSearch(overlap_graph)

    for name in overlap_graph.anchors:
        try:
            state=search.search(Node(name),Node('ctg3'))
            print('Search done.')
            print('Reconstructing the genome...')
            genome=overlap_graph.reconstruct_path(state)
            print('Reconstruction done.')
            print('Saving genome...')
            FASTAReader.save('kita',genome,name+'kurcina.fasta')
            print('Save done.')
        except:
            print('Nije %s'%(name))
