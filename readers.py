
class Overlap():
    def __init__(self,args):
        self.args=args[:12]
        arg_names=[('query_name',str),('query_length',int),('query_start',int),('query_end',int),
                    ('relative_strand',str),('target_name',str),('target_length',int),('target_start',int),
                    ('target_end',int),('num_matches',int),('block_length',int),('mapping_quality',int)
                    ]

        for block,arg in zip(arg_names,args):
            setattr(self,block[0],block[1](arg))


class PAFReader():
    def __init__(self,filepath):
        self.overlaps=[]
        with open(filepath,'r') as file:
            for line in file.readlines():
                args=[arg.strip() for arg in line.split('\t')]
                self.overlaps.append(Overlap(args))

if __name__=='__main__':
    overlaps=PAFReader('overlaps.paf')
    print(overlaps.overlaps[69].args)
