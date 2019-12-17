
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

class FASTAReader():
    def __init__(self,filepath):
        self.reads={}
        with open(filepath,'r') as file:
            lines=file.readlines()
            for i in range(len(lines)//2):
                read_name=lines[i*2][1:].strip()
                read=lines[i*2+1].strip()
                self.reads[read_name]=read
    def save(node_name,genome,filepath):
        with open(filepath,'w') as file:
            file.write('>'+node_name+'\n')
            file.write(genome)

if __name__=='__main__':
    overlaps=PAFReader('overlaps.paf')
    print(overlaps.overlaps[69].args)
