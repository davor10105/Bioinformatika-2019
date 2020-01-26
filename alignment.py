
class Utils():
    """
    Contains methods for calculating various properties of sequences
    and overlaps.
    """

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
        '''
        if overlap.query_length>overlap.target_length:
            longer_sequence=(overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name)
            shorter_sequence=(overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name)
        else:
            shorter_sequence=(overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name)
            longer_sequence=(overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name)

        if overlap.relative_strand=='+':
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
        else:
            #shorter_sequence=(shorter_sequence[0],shorter_sequence[0]-shorter_sequence[2],shorter_sequence[0]-shorter_sequence[1],shorter_sequence[3])
            if shorter_sequence[0]-shorter_sequence[2]-0 > longer_sequence[1]-0:
                #left overhang
                overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name=longer_sequence[0],longer_sequence[1],longer_sequence[2],longer_sequence[3]
                overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name=shorter_sequence[0],shorter_sequence[1],shorter_sequence[2],shorter_sequence[3]

                return overlap
            elif shorter_sequence[0]-(shorter_sequence[0]-shorter_sequence[1]) > longer_sequence[0]-longer_sequence[2]:
                #right overhang
                overlap.query_length,overlap.query_start,overlap.query_end,overlap.query_name=shorter_sequence[0],shorter_sequence[1],shorter_sequence[2],shorter_sequence[3]
                overlap.target_length,overlap.target_start,overlap.target_end,overlap.target_name=longer_sequence[0],longer_sequence[1],longer_sequence[2],longer_sequence[3]

                return overlap
            return None

    def check_confidence(overlap,cutoff=0.25):
        '''
        Checks if the overlap is within a certain threshold.
        '''
        if overlap.num_matches/overlap.block_length<cutoff:
            return False
        return True

    def get_overlap_score(overlap):
        '''
        Calculates the overlap score.
        :return: overlap score
        '''
        OL1=overlap.query_end-overlap.query_start
        OL2=overlap.target_end-overlap.target_start
        SI=overlap.block_length

        return (OL1+OL2)*SI/2

    def get_extension_scores(overlap,overlap_score):
        '''
        Calculates the extension score given overlap and overlap score.
        :param overlap_score:
        '''
        OH1=overlap.query_length-overlap.query_end
        OH2=overlap.target_start-0
        EL1=overlap.query_start-0
        EL2=overlap.target_length-overlap.target_end

        ES1=overlap_score+EL1/2-(OH1+OH2)/2
        ES2=overlap_score+EL2/2-(OH1+OH2)/2

        return ES2,ES1

    def reverse_complement(sequence):
        '''
        Returns a reverse complement sequence
        '''
        sequence = Utils.complement(sequence)
        return sequence[::-1]

    def complement(sequence):
        '''
        Returns a complement sequence for a given orginal sequence.
        '''
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        bases = list(sequence)
        bases = [complement[base] for base in bases]
        return ''.join(bases)


if __name__=='__main__':
    print(Utils.reverse_complement('AATTTGGCCCTT'))
