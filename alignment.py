from readers import Overlap

class Utils():
    def get_extension_overlaps(overlaps):
        '''
        Returns only the overlaps that aren't contained, for example:
        this is OK:
            node1:    AAATTT
            node2:      ATTTAA
        this is not OK:
            node1:    AAATTT
            node2:      ATT

        Additionally, checks the confidence of the overlaps and returns only
        those above the threshold. (OVO NISAM SIGURAN)
        '''
        extension_overlaps=[]
        unconfident=0
        for overlap in overlaps:
            left_right_overlap=Utils.check_if_contained(overlap)
            if left_right_overlap is not None:
                if Utils.check_confidence(overlap):
                    extension_overlaps.append(left_right_overlap)
                else:
                    unconfident+=1

        print("Unconfident")
        print(unconfident)
        print("Confident")
        print(len(overlaps)-unconfident)

        return extension_overlaps

    def check_if_contained(overlap):
        '''
        Checks if the overlap is contained and returns None if it is.
        Otherwise, rearranges the the overlapping nodes so that they always go left
        to right in the Overlap object (query is the left node and target is the
        right node). For example:

        query:    contig1:    AAATTTTTTTTTAA
        target:    read1:  AAAAAATTT

        will become:
        query:    read1:   AAAAAATTT
        target:    contig1:   AAATTTTTTTTTAA

        Da malo bolje pojasnim, kada se u minimap ubace contig.fasta i reads.fasta,
        on ce overlape uvijek dati u tom redoslijedu, znaci query ce uvijek biti node
        iz contiga, a u targetu uvijek node iz reada. E sad, to bas nije prakticno
        provjeravati unutar koda drugdje pa se ovdje kao query uvijek stavi lijeviji
        node, a u target desniji node bez obzira je li contig ili read.
        '''
        #TODO: MOZDA TREBA DRUGACIJE
        if overlap.query_length>overlap.target_length:
            longer_sequence=(overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name)
            shorter_sequence=(overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name)
        else:
            shorter_sequence=(overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name)
            longer_sequence=(overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name)

        if shorter_sequence[1]-0 > longer_sequence[1]-0:
            #left overhang
            overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name=shorter_sequence[0],shorter_sequence[1],shorter_sequence[2],shorter_sequence[3]
            overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name=longer_sequence[0],longer_sequence[1],longer_sequence[2],longer_sequence[3]
            return overlap
        elif shorter_sequence[0]-shorter_sequence[2] > longer_sequence[0]-longer_sequence[2]:
            #right overhang
            overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name=longer_sequence[0],longer_sequence[1],longer_sequence[2],longer_sequence[3]
            overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name=shorter_sequence[0],shorter_sequence[1],shorter_sequence[2],shorter_sequence[3]
            return overlap
        return None

    def check_confidence(overlap,cutoff=0.):
        '''
        Checks if the overlap is within a certain threshold.
        '''
        #TODO: NEMAM POJMA RACUNA LI SE TO OVAKO
        #if overlap.block_length<40000:
        #   return False
        if overlap.num_matches/overlap.block_length<cutoff:
            return False
        return True

    def get_overlap_score(overlap):
        '''
        Calculates the overlap score.
        '''
        OL1=overlap.query_end-overlap.query_start
        OL2=overlap.target_end-overlap.target_start
        SI=overlap.block_length

        return (OL1+OL2)*SI/2

    def get_extension_scores(overlap,overlap_score):
        '''
        Calculates the extension score.
        '''
        OH1=overlap.query_length-overlap.query_end
        OH2=overlap.target_start-0
        EL1=overlap.query_start-0
        EL2=overlap.target_length-overlap.target_end

        ES1=overlap_score+EL1/2-(OH1+OH2)/2
        ES2=overlap_score+EL2/2-(OH1+OH2)/2

        return ES2,ES1
