from alignment import Utils


class Overlap():
    '''
    Defines an overlap, basically PAF format overlap in a class.
    '''

    def __init__(self, args):
        self.args = args[:12]
        arg_names = [('query_name', str), ('query_length', int), ('query_start', int), ('query_end', int),
                     ('relative_strand', str), ('target_name', str), ('target_length', int), ('target_start', int),
                     ('target_end', int), ('num_matches', int), ('block_length', int), ('mapping_quality', int)
                     ]

        for block, arg in zip(arg_names, args):
            setattr(self, block[0], block[1](arg))

    def __str__(self):
        return '(%s,%s,%s) (%s,%s,%s)' % (
        self.query_length, self.query_start, self.query_end, self.target_length, self.target_start, self.target_end)


class PAFReader():
    '''
    Reads PAF formatted file.
    '''

    def __init__(self, filepath, overlap_confidence):
        self.overlaps = []
        with open(filepath, 'r') as file:
            count = 0
            for line in file.readlines():
                count += 1
                args = [arg.strip() for arg in line.split('\t')]
                new_overlap = Utils.check_if_contained(Overlap(args))
                if new_overlap != None:
                    if Utils.check_confidence(new_overlap, overlap_confidence):
                        self.overlaps.append(new_overlap)

                if count % 10000 == 0:
                    print("Broj obradenih overlapa:")
                    print(count)

    def job(self, lines, index_from_inclusive, index_to_exclusive):
        '''
        Loads overlaps for a given lines of a file in interval [a, b>.
        :param lines: read lines
        :param index_from_inclusive:
        :param index_to_exclusive:
        :return:
        '''
        for line in lines[index_from_inclusive:min(index_to_exclusive, len(lines))]:
            args = [arg.strip() for arg in line.split('\t')]
            new_overlap = Utils.check_if_contained(Overlap(args))
            if new_overlap != None:
                if Utils.check_confidence(new_overlap):
                    self.overlaps.append(new_overlap)
        self.finished += 1

    def multithreaded_read(self, n_threads, lines):
        '''
        Multithreaded version of loading overlaps into the list.
        :param n_threads:
        :param lines:
        :return:
        '''
        import threading
        job_size = len(lines) // n_threads
        ts = [threading.Thread(target=self.job, args=(lines, i * job_size, (i + 1) * job_size)) for i in range(n_threads)]
        [(t.start(), t.join()) for t in ts]


class FASTAReader():
    '''
    Reads a FASTA file and constructs a dictionary:
    {'node_name': NODE_GENOME (string)}
    '''

    def __init__(self, filepath):
        self.reads = {}
        with open(filepath, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines) // 2):
                read_name = lines[i * 2][1:].split()[0].strip()
                read = lines[i * 2 + 1].strip()
                self.reads[read_name] = read

    def save(node_names, genomes, filepath):
        with open(filepath, 'w') as file:
            for node_name, genome in zip(node_names, genomes):
                file.write('>' + node_name + '\n')
                file.write(genome + '\n')


class FASTQReader():
    '''
    Reads a FASTAQ file and constructs a dictionary:
    {'node_name': NODE_GENOME (string)}
    '''

    def __init__(self, filepath):
        self.reads = {}
        with open(filepath, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines) // 4):
                read_name = lines[i * 4][1:].strip()
                read = lines[i * 4 + 1].strip()
                self.reads[read_name] = read


if __name__ == '__main__':
    overlaps = PAFReader('overlaps.paf')
    print(overlaps.overlaps[69].args)
