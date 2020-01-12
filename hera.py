from search import *
import argparse

class HERA():
    def run(contigs_path,reads_path,contig_read_overlap_path,read_read_overlap_path,
            output_path,OVERLAP_CONFIDENCE=0.25,MAX_NODE_LENGTH=1000,MAX_OPEN_LEN=1000,
            MAX_NO_CHANGE=500,USED_NODE_WEIGHT=1,SCORE_WEIGHT=1.2):

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
        overlap_graph=Graph(contigs_path,reads_path,contig_read_overlap_path,read_read_overlap_path,overlap_confidence=OVERLAP_CONFIDENCE)
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
                    #FASTAReader.save('kita',genome,name+ORGANISM_NAME+'.fasta')
                    #print('Got one')
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
        FASTAReader.save(['ctg'+str(i) for i in range(len(genomes))],genomes,output_path)
        print('Save done.')

def parse_arguments():
    # def run(contigs_path,reads_path,contig_read_overlap_path,read_read_overlap_path,
    #        output_path,OVERLAP_CONFIDENCE=0.25,MAX_NODE_LENGTH=1000,MAX_OPEN_LEN=1000,
    #        MAX_NO_CHANGE=500,USED_NODE_WEIGHT=1,SCORE_WEIGHT=1.2):

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs_path',type=str,help="Contig file path")
    parser.add_argument('reads_path',type=str,help="Read file path")
    parser.add_argument('contig_read_overlap_path',type=str,help="Contig read overlap file path")
    parser.add_argument('read_read_overlap_path',type=str,help="Contig read overlap file path")
    parser.add_argument('output_path',type=str,help="Output file path")
    parser.add_argument('--OVERLAP_CONFIDENCE',type=float,default=0.25,help="Specify overlape confidence")
    parser.add_argument('--MAX_NODE_LENGTH',type=int,default=1000,help="Specify max node length")
    parser.add_argument('--MAX_OPEN_LEN',type=int,default=1000,help="")
    parser.add_argument('--MAX_NO_CHANGE',type=int,default=500,help="")
    parser.add_argument('--USED_NODE_WEIGHT',type=float,default=1,help="")
    parser.add_argument('--SCORE_WEIGHT',type=float,default=1.2,help="")
    args=parser.parse_args()
    return args

if __name__=='__main__':
    '''
    args=parse_arguments()
    HERA.run(args.contigs_path, args.reads_path, args.contig_read_overlap_path, args.read_read_overlap_path,
            args.output_path, args.OVERLAP_CONFIDENCE, args.MAX_NODE_LENGTH,args. MAX_OPEN_LEN,
            args.MAX_NO_CHANGE, args.USED_NODE_WEIGHT, args.SCORE_WEIGHT)
    '''
    HERA.run('./data/cjejuni/contigs.fasta','./data/cjejuni/reads.fastq','./data/cjejuni/contig_read.paf','./data/cjejuni/read_read.paf','cjejuni.fasta')
    #HERA.run('EC_contigs.fasta','EC_reads.fasta','EC_cr_overlaps.paf','EC_rr_overlaps.paf','ec_novo.fasta')
