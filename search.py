from graph import *
from readers import *
import argparse

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

    def state_cmp(max_used_nodes,max_score):
        def f(state):
            return state.anchors_found,-1.*len(state.used_nodes)/max_used_nodes+1.2*state.score/max_score

        return f

    def search(self,start_node):
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

            if len(current_state.used_nodes)>1000:
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
                if edge.node_to not in current_state.used_nodes:
                    new_used_nodes=set(current_state.used_nodes)
                    new_used_nodes.add(edge.node_to)
                    new_anchor_num=current_state.anchors_found
                    if edge.node_to.name in anchor_names:
                        new_anchor_num+=1
                    new_state=State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction,new_used_nodes,edge,new_anchor_num)
                    open.insert(index,new_state)
                    index+=1

            max_used_nodes,max_score=DFSSearch.get_maxes(open)
            open=sorted(open,key=DFSSearch.state_cmp(max_used_nodes,max_score),reverse=True)[:5000]
            if no_change>500:
                return open[0]

        return current_state

'''
Dodano za kasnije
    potencijalno premjestiti negdje
    ili staviti u main
'''
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('contig_path',type=str,help="Contig file path")
    parser.add_argument('read_path',type=str,help="Read file path")
    parser.add_argument('contig_read_overlap_path',type=str,help="Contig read overlap file path")
    parser.add_argument('read_read_overlap_path',type=str,help="Contig read overlap file path")
    args=parser.parse_args()
    #print(args.contig_path)
    #print(args.read_path)
    return args

if __name__=='__main__':
    '''
    ZA ECOLI RADI ZA 500 NO_CHANGE
    1 ZA LEN
    1.2 ZA SCORE
    5000 cutoff open liste
    ubacivalo se samo 50 nodeova iz edges (radi i za sve)
    '''

    '''
    ZA CJEJUNI OKI RADI ZA 1000 NO_CHANGE ili 2000
    -2 ZA LEN
    0. ZA SCORE ili 0.5
    100 cutoff open liste
    ubacivali su se svi nodeovi iz edges
    '''

    print('Creating the overlap graph...')
    '''
    needs contig.fasta reads.fasta overlaps_contig_reads and overlaps_reads_reads files
    to construct the overlap graph
    '''
    #args=parse_arguments()
    #overlap_graph=Graph(args.contig_path,args.read_path,args.contig_read_overlap_path,args.read_read_overlap_path)
    overlap_graph=Graph('./data/ecoli/contigs.fasta','./data/ecoli/reads.fasta','./data/ecoli/contig_read.paf','./data/ecoli/read_read.paf')
    print('Overlap graph done.')
    print('Starting the search...')
    search=DFSSearch(overlap_graph)

    '''
    searchu se predaje Node koji je onda zapravo starting_node, u ovoj for petlji
    je kako se to radi iterativno za svaki contig, onda samo umjesto 'ctg3' stavite
    name kako bi se krenulo od tog Nodea.
    '''
    best_states=[]
    for name in overlap_graph.anchors:
        try:
            state=search.search(Node(name))
            best_states.append(state)
            genome=overlap_graph.reconstruct_path(state)
            #current_best_state=state.compare(current_best_state)
            FASTAReader.save('kita',genome,name+'ecoli.fasta')
            print('Got one')
        except:
            print('Nije %s'%(name))

    print("Duljine nadenih puteva:")
    for state in best_states:
        print(state.anchors_found)

    max_used_nodes,max_score=DFSSearch.get_maxes(best_states)
    best_state=sorted(best_states,key=DFSSearch.state_cmp(max_used_nodes,max_score),reverse=True)[0]

    print('Search done.')
    print('Reconstructing the genome...')
    genome=overlap_graph.reconstruct_path(best_state)
    print('Reconstruction done.')
    print('Saving genome...')

    '''
    saves the found genome in file kita2.fasta and names the genome 'kita'
    '''
    FASTAReader.save('kita',genome,'ecoli.fasta')
    print('Save done.')
