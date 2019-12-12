from graph import *

class CostSearch():
    def __init__(self,graph):
        self.graph=graph

    def search(self,start_node,end_node):
        start_state=State(start_node,None,0,None)

        visited=set()
        open=[start_state]
        current_state=start_state
        while current_state.node!=end_node:
            print(len(open))
            if len(open)==0:
                raise ValueError('No path to end node')

            current_state=open.pop(0)
            if current_state in visited:
                continue
            visited.add(current_state)

            for edge in self.graph.get_edges(current_state.node,current_state.direction):
                open.append(State(edge.node_to,current_state,current_state.score+edge.overlap_score,edge.direction))

            open=sorted(open,key=lambda state: state.score,reverse=True)

        return current_state

if __name__=='__main__':
    overlap_graph=Graph('overlaps_mjau.paf','read_overlaps.paf')
    search=CostSearch(overlap_graph)
    print(search.search(Node('ctg1'),Node('ctg2')))
