class Overlap():
    def get_extension_overlaps(overlaps):
        extension_overlaps=[]
        for overlap in overlaps:
            if not Overlap.check_if_contained(overlap.query_length,overlap.query_start,overlap.query_end) and not Overlap.check_if_contained(overlap.target_length,overlap.target_start,overlap.target_end):
                extension_overlaps.append(overlap)

        return extension_overlaps

    def check_if_contained(sequence_length,overlap_start,overlap_end):
        if overlap_end-overlap_start<sequence_length:
            return False
        return True

    def get_overlap_score(overlap):
        OL1=overlap.query_end-overlap.query_start
        OL2=overlap.target_end-overlap.target_start
        SI=overlap.block_length

        return (OL1+OL2)*SI/2

    def get_extension_score(overlap,overlap_score,extends_left):
        pass
