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

        upper_length_threshold=0.75*sum([len(self.graph.anchors[anchor]) for anchor in self.graph.anchors])/len(self.graph.anchors)

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

            intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
            print('Number of found contigs:')
            print(intersecting_anchors)
            print(len(current_state.used_nodes))
            #print(current_state.anchors_found)

            if len(current_state.used_nodes)>3000:
                while current_state.previous_state.node.name not in anchor_names:
                    current_state=current_state.previous_state
                return current_state
            '''
            if current_state.direction!=None:
                path_length=len(self.graph.reconstruct_path(current_state))
                intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
                anchor_length=sum([len(self.graph.anchors[anchor]) for anchor in intersecting_anchors])
                if len(self.graph.get_extension_genome(current_state))>upper_length_threshold:
                    print('Number of found contigs:')
                    print(len(intersecting_anchors))
                    #while current_state.previous_state.node.name not in anchor_names:
                    #    current_state=current_state.previous_state

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
            '''


            #if current_state in visited:
            #    continue
            #visited.add(current_state)
            edges=sorted(self.graph.get_edges(current_state.node,current_state.direction),key=lambda edge: edge.overlap_score,reverse=True)
            anchor_edges=[]
            read_edges=[]
            for edge in edges:
                if edge.node_to.name in anchor_names:
                    anchor_edges.append(edge)
                else:
                    read_edges.append(edge)
            anchor_edges=sorted(anchor_edges,key=lambda edge: edge.overlap_score,reverse=True)
            read_edges=sorted(read_edges,key=lambda edge: edge.overlap_score,reverse=True)
            if current_state.node.name not in anchor_names and len(anchor_edges)>0:
                store_edges=anchor_edges
            else:
                store_edges=read_edges

            for edge in edges[:50]:
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
                    new_state=State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction,new_used_nodes,edge,new_anchor_num)
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
    return state.score,state.anchors_found
    #return state.anchors_found,state.score
    return state.score

class DFSSearch():
    def __init__(self,graph):
        self.graph=graph

    def get_maxes(open):
        max_used_nodes=len(max(open,key=lambda state:len(state.used_nodes)).used_nodes)
        max_score=max(open,key=lambda state:state.score).score

        return max_used_nodes,max_score

    def state_cmp(max_used_nodes,max_score,used_node_weight,score_weight):
        def f(state):
            if max_used_nodes==0:
                cur_used_score=0
            else:
                cur_used_score=len(state.used_nodes)/max_used_nodes
            if max_score==0:
                cur_score_score=0
            else:
                cur_score_score=state.score/max_score
            return state.anchors_found,-used_node_weight*cur_used_score+score_weight*cur_score_score

        return f

    def search(self,start_node,found_anchors,max_node_length=1000,max_open_len=200,max_no_change=1000,used_node_weight=1.,score_weight=1.2):
        start_state=State(start_node,None,0,None,{start_node},None,1)
        open=[start_state]
        current_state=start_state

        upper_length_threshold=0.25*sum([len(self.graph.anchors[anchor]) for anchor in self.graph.anchors])/len(self.graph.anchors)

        anchor_names=set(self.graph.anchors.keys())
        current_best_state=current_state
        no_change=0
        while len(anchor_names)-current_state.anchors_found!=0:

            if len(open)==0:
                raise ValueError('No path to end node')

            current_state=open.pop(0)

            if current_best_state.anchors_found>=current_state.anchors_found:
                no_change+=1
            else:
                current_best_state=current_state
                no_change=0

            #if current_state.direction!=None:
            #    if len(self.graph.get_extension_genome(current_state))>upper_length_threshold:
            #        continue

            if len(current_state.used_nodes)>max_node_length:
                continue

            intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
            print('Number of found contigs:')
            print(intersecting_anchors)
            print(len(current_state.used_nodes))
            print(no_change)

            edges=sorted(self.graph.get_edges(current_state.node,current_state.direction),key=lambda edge: edge.overlap_score,reverse=True)
            anchor_edges=[]
            read_edges=[]
            for edge in edges:
                if edge.node_to.name in anchor_names:
                    anchor_edges.append(edge)
                else:
                    read_edges.append(edge)
            anchor_edges=sorted(anchor_edges,key=lambda edge: edge.overlap_score,reverse=True)[:max(1,int(len(anchor_edges)*0.5))]
            read_edges=sorted(read_edges,key=lambda edge: edge.overlap_score,reverse=True)[:max(1,int(len(read_edges)*0.5))]
            if current_state.node.name not in anchor_names and len(anchor_edges)>0:
                store_edges=anchor_edges
            else:
                store_edges=read_edges

            index=0
            #for edge in edges[:50]:
            for edge in edges:
                if edge.node_to not in current_state.used_nodes and edge.node_to.name not in found_anchors:
                    new_used_nodes=set(current_state.used_nodes)
                    new_used_nodes.add(edge.node_to)
                    new_anchor_num=current_state.anchors_found
                    if edge.node_to.name in anchor_names:
                        new_anchor_num+=1
                    new_state=State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction,new_used_nodes,edge,new_anchor_num)
                    open.insert(index,new_state)
                    index+=1

            max_used_nodes,max_score=DFSSearch.get_maxes(open)
            open=sorted(open,key=DFSSearch.state_cmp(max_used_nodes,max_score,used_node_weight,score_weight),reverse=True)[:max_open_len]
            if no_change>max_no_change:
                best_state=open[0]
                while best_state.node.name not in anchor_names:
                    best_state=best_state.previous_state
                return best_state

        return current_state


if __name__=='__main__':
    MAX_NODE_LENGTH=1000
    MAX_OPEN_LEN=1000
    MAX_NO_CHANGE=500
    USED_NODE_WEIGHT=1
    SCORE_WEIGHT=1.2

    ORGANISM_NAME='cjejuni'
    '''
    ecoli

    MAX_NODE_LENGTH=1000
    MAX_OPEN_LEN=1000
    MAX_NO_CHANGE=500
    USED_NODE_WEIGHT=1
    SCORE_WEIGHT=1.2

    0.25 cutoff
    svi edgeovi
    '''

    '''
    cjejuni

    MAX_NODE_LENGTH=1000
    MAX_OPEN_LEN=1000
    MAX_NO_CHANGE=500
    USED_NODE_WEIGHT=1
    SCORE_WEIGHT=1.2

    0.25 cutoff
    svi edgeovi
    '''

    '''
    za bgram

    MAX_NODE_LENGTH=1000
    MAX_OPEN_LEN=200
    MAX_NO_CHANGE=1000
    USED_NODE_WEIGHT=2.
    SCORE_WEIGHT=0.5
    '''

    print('Creating the overlap graph...')
    '''
    needs contig.fasta reads.fasta overlaps_contig_reads and overlaps_reads_reads files
    to construct the overlap graph
    '''
    overlap_graph=Graph('./data/cjejuni/contigs.fasta','./data/cjejuni/reads.fastq','./data/cjejuni/contig_read.paf','./data/cjejuni/read_read.paf')
    print('Overlap graph done.')
    print('Starting the search...')
    search=DFSSearch(overlap_graph)

    '''
    If not all contigs were found.
    '''
    all_anchors=set(overlap_graph.anchors.keys())
    found_anchors=set()
    whole_path=[]
    while len(overlap_graph.anchors)-len(found_anchors)>0:
        best_states=[]
        for name in set(overlap_graph.anchors.keys()).difference(found_anchors):
            try:
                state=search.search(Node(name),found_anchors=found_anchors,max_node_length=MAX_NODE_LENGTH,max_open_len=MAX_OPEN_LEN,max_no_change=MAX_NO_CHANGE,used_node_weight=USED_NODE_WEIGHT,score_weight=SCORE_WEIGHT)
                best_states.append(state)
                genome=overlap_graph.reconstruct_path(state)
                #current_best_state=state.compare(current_best_state)
                FASTAReader.save('kita',genome,name+ORGANISM_NAME+'.fasta')
                print('Got one')
            except:
                print('Nije %s'%(name))
        max_used_nodes,max_score=DFSSearch.get_maxes(best_states)
        best_state=sorted(best_states,key=DFSSearch.state_cmp(max_used_nodes,max_score,USED_NODE_WEIGHT,SCORE_WEIGHT),reverse=True)[0]
        new_found_anchors=all_anchors.intersection(set([node.name for node in best_state.used_nodes]))
        found_anchors=found_anchors.union(new_found_anchors)

        whole_path.append(best_state)

    for path in whole_path:
        print(all_anchors.intersection(set([node.name for node in path.used_nodes])))


    print('Search done.')
    print('Reconstructing the genome...')
    genomes=[]
    starting_switch_strand=False
    for path in whole_path:
        genome,starting_switch_strand=overlap_graph.reconstruct_path(path,starting_switch_strand=starting_switch_strand,return_strand=True)
        genomes.append(genome)
    print('Reconstruction done.')
    print('Saving genome...')

    '''
    saves the found genome in file kita2.fasta and names the genome 'kita'
    '''
    FASTAReader.save(['ctg'+str(i) for i in range(len(genomes))],genomes,ORGANISM_NAME+'.fasta')
    print('Save done.')
