from readers import Overlap

class Utils():
    def get_extension_overlaps(overlaps):
        extension_overlaps=[]
        for overlap in overlaps:
            left_right_overlap=Utils.check_if_contained(overlap)
            if left_right_overlap is not None:
                extension_overlaps.append(left_right_overlap)

        return extension_overlaps

    def check_if_contained(overlap):
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

    def get_overlap_score(overlap):
        OL1=overlap.query_end-overlap.query_start
        OL2=overlap.target_end-overlap.target_start
        SI=overlap.block_length

        return (OL1+OL2)*SI/2

    def get_extension_scores(overlap,overlap_score):
        OH1=overlap.query_length-overlap.query_end
        OH2=overlap.target_start-0
        EL1=overlap.query_start-0
        EL2=overlap.target_length-overlap.target_end

        ES1=overlap_score+EL1/2-(OH1+OH2)/2
        ES2=overlap_score+EL2/2-(OH1+OH2)/2

        return ES2,ES1
