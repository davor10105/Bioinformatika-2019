
class Overlap():
    '''
    Defines an overlap, basically PAF format overlap in a class.
    '''
    def __init__(self,args):
        self.args=args[:12]
        arg_names=[('query_name',str),('query_length',int),('query_start',int),('query_end',int),
                    ('relative_strand',str),('target_name',str),('target_length',int),('target_start',int),
                    ('target_end',int),('num_matches',int),('block_length',int),('mapping_quality',int)
                    ]

        for block,arg in zip(arg_names,args):
            setattr(self,block[0],block[1](arg))

    def __str__(self):
        return '(%s,%s,%s) (%s,%s,%s)'%(self.query_length,self.query_start,self.query_end,self.target_length,self.target_start,self.target_end)


class PAFReader():
    '''
    Reads PAF formatted file.
    '''
    def __init__(self,filepath):
        self.overlaps=[]
        with open(filepath,'r') as file:
            for line in file.readlines():
                args=[arg.strip() for arg in line.split('\t')]
                self.overlaps.append(Overlap(args))

#TODO: TREBA I FASTAQ READER
class FASTAReader():
    '''
    Reads a FASTA file and constructs a dictionary:
    {'node_name': NODE_GENOME (string)}
    '''
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
