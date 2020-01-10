from graph import *
from readers import *

class CostSearch():
    '''
    Implementation of the cost search algorithm.
    '''
    def __init__(self,graph):
        self.graph=graph

    def search(self,start_node):
        '''
        Starts the search from Node start_node, adding each new state to open list
        that is sorted by score (either overlap or extension) from higher to lower values
        (we are trying to maximize the score) until we reach the State that visited
        all contigs (anchor nodes)
        '''
        start_state=State(start_node,None,0,None,{start_node},None,1)
        '''
        start_state is (starting_node, previous state and edge are None, starting score is 0
        and direction is None so that the path can go either way in the beginning).
        '''

        visited=set()
        open=[start_state]
        current_state=start_state

        upper_length_threshold=0.15*sum([len(self.graph.anchors[anchor]) for anchor in self.graph.anchors])

        anchor_names=set(self.graph.anchors.keys())
        #while len(anchor_names-set([node.name for node in current_state.used_nodes]))!=0:
        while len(anchor_names)-current_state.anchors_found!=0:
            '''
            Check if current state visited all contigs
            '''
            #print(current_state.node)
            #print(end_node)
            #print(len(open))
            if len(open)==0:
                '''
                if open list is empty, the path cannot be found
                '''
                raise ValueError('No path to end node')

            current_state=open.pop(0)

            #intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
            print('Number of found contigs:')
            #print(len(intersecting_anchors))
            print(current_state.anchors_found)


            if current_state.direction!=None:
                path_length=len(self.graph.reconstruct_path(current_state))
                intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
                anchor_length=sum([len(self.graph.anchors[anchor]) for anchor in intersecting_anchors])
                if path_length-anchor_length>upper_length_threshold:
                    print('Number of found contigs:')
                    print(len(intersecting_anchors))
                    while current_state.previous_state.node.name not in anchor_names:
                        current_state=current_state.previous_state

                    open=sorted(open,key=lambda state:len(anchor_names.intersection(set([node.name for node in state.used_nodes]))),reverse=True)
                    print('Node s najvecim brojem contiga:')
                    print(len(anchor_names.intersection(set([node.name for node in open[0].used_nodes]))))

                    temp_state=current_state
                    while temp_state.previous_state!=None:
                        if temp_state.node.name in anchor_names:
                            print(temp_state.node.name,end=' ')
                        temp_state=temp_state.previous_state
                    print(temp_state.node.name,end=' ')
                    print()

                    return current_state
                    #raise ValueError('Path is too long')


            #if current_state in visited:
            #    continue
            #visited.add(current_state)

            for edge in self.graph.get_edges(current_state.node,current_state.direction):
                '''
                get all edges from current_state in current_state's 'preferred' direction
                '''
                if edge.node_to not in current_state.used_nodes:
                    '''
                    if the node_to wasn't visited before, add the new state to open list
                    (prevents going through the same node twice)
                    '''
                    new_used_nodes=set(current_state.used_nodes)
                    new_used_nodes.add(edge.node_to)
                    new_anchor_num=current_state.anchors_found
                    if edge.node_to.name in anchor_names:
                        new_anchor_num+=1
                    new_state=State(edge.node_to,current_state,current_state.score+edge.extension_score,edge.direction,new_used_nodes,edge,new_anchor_num)
                    '''
                    Ovdje se moze mijenjati iz +edge.overlap_score u +edge.extension score
                    kako bi se koristio jedan ili drugi scoring, ali mislim da je skoro
                    pa svejedno, za ecoli se dobije isto.
                    '''
                    #if new_state not in visited:
                    open.append(new_state)

            '''
            Mislim da bi se ovdje nekako moglo registrirati da se put ne moze naci
            prije nego se memorija popuni (jer open lista ima ogroman broj stateova).
            '''
            open=sorted(open,key=state_key,reverse=True)
            #print(len(open))

        #print(current_state.used_nodes)
        return current_state

def state_key(state):
    return state.anchors_found,state.score

if __name__=='__main__':
    print('Creating the overlap graph...')
    '''
    needs contig.fasta reads.fasta overlaps_contig_reads and overlaps_reads_reads files
    to construct the overlap graph
    '''
    overlap_graph=Graph('./data/cjejuni/contigs.fasta','./data/cjejuni/reads.fastq','./data/cjejuni/contig_read.paf','./data/cjejuni/read_read.paf')
    print('Overlap graph done.')
    print('Starting the search...')
    search=CostSearch(overlap_graph)

    '''
    searchu se predaje Node koji je onda zapravo starting_node, u ovoj for petlji
    je kako se to radi iterativno za svaki contig, onda samo umjesto 'ctg3' stavite
    name kako bi se krenulo od tog Nodea.
    '''
    current_best_state=None
    for name in overlap_graph.anchors:
        try:
            state=search.search(Node(name))
            current_best_state=state.compare(current_best_state)
        except:
            print('Nije %s'%(name))

    print('Search done.')
    print('Reconstructing the genome...')
    genome=overlap_graph.reconstruct_path(current_best_state)
    print('Reconstruction done.')
    print('Saving genome...')

    '''
    saves the found genome in file kita2.fasta and names the genome 'kita'
    '''
    FASTAReader.save('kita',genome,'cjejuni.fasta')
    print('Save done.')
