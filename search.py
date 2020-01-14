from graph import *
from readers import *

def state_key(state):
    return state.score,state.anchors_found
    return state.score

class CostSearch():
    def __init__(self,graph):
        self.graph=graph

    def get_maxes(open):
        '''
        Calculates number of used nodes of a state with largest number of used nodes
        and max score value for all states.
        '''
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

        anchor_names=set(self.graph.anchors.keys())
        current_best_state=current_state
        no_change=0

        while len(anchor_names)-current_state.anchors_found!=0:

            if len(open)==0:
                raise ValueError('No path to end node')

            current_state=open.pop(0)

            '''
            Check number of current found anchors and set change
            '''
            if current_best_state.anchors_found>=current_state.anchors_found:
                no_change+=1
            else:
                current_best_state=current_state
                no_change=0

            if len(current_state.used_nodes)>max_node_length:
                continue

            intersecting_anchors=anchor_names.intersection(set([node.name for node in current_state.used_nodes]))
            print('Number of found contigs:')
            print(intersecting_anchors)
            print(len(current_state.used_nodes))
            print(no_change)

            '''
            Sort edges of a current node based on overlap score
            '''
            edges=sorted(self.graph.get_edges(current_state.node,current_state.direction),key=lambda edge: edge.overlap_score,reverse=True)

            index=0
            for edge in edges:
                if edge.node_to not in current_state.used_nodes and edge.node_to.name not in found_anchors:
                    '''
                    Create new state for each new node connected to node in current state
                    '''
                    new_used_nodes=set(current_state.used_nodes)
                    new_used_nodes.add(edge.node_to)
                    new_anchor_num=current_state.anchors_found
                    if edge.node_to.name in anchor_names:
                        new_anchor_num+=1
                    new_state=State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction,new_used_nodes,edge,new_anchor_num)
                    open.insert(index,new_state)
                    index+=1

            max_used_nodes,max_score=CostSearch.get_maxes(open)
            '''
            Sort all states based on score and number of used nodes
            '''
            open=sorted(open,key=CostSearch.state_cmp(max_used_nodes,max_score,used_node_weight,score_weight),reverse=True)[:max_open_len]

            if no_change>max_no_change:
                '''
                Return last state with anchor node from sorted list of states
                '''
                best_state=open[0]
                while best_state.node.name not in anchor_names:
                    best_state=best_state.previous_state
                return best_state

        return current_state
