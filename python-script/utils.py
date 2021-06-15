from past.builtins import xrange
from Bio    import SeqIO

#import RNA

complement_table = str.maketrans('ATGCU', 'TACGA')

def stream_fasta_seq_list(fasta_filename):
    with open(fasta_filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield str(record.seq)

def get_fasta_seq_list(fasta_filename):
    return list(stream_fasta_seq_list(fasta_filename))

def stream_txt_seq_list(text_filename):
    with open(text_filename) as infile:
        for line in infile:
            yield line.strip()

def get_txt_seq_list(text_filename):
    return list(stream_txt_seq_list(text_filename))

def uniquify_background_list(background_list):
    uniq_background_set = set()
    while background_list:
        uniq_background_set.add(background_list.pop())
    background_list = []
    while uniq_background_set:
        background_list.append(uniq_background_set.pop())
    return background_list

def stream_kmers(seq, k):
    if k >= len(seq):
        return [seq]
    return (seq[i:i+k] for i in xrange(len(seq)-k+1))

def get_comp(seq):
    return seq.translate(complement_table)

def get_revcomp(seq):
    return get_comp(seq)[::-1]

def stream_min_kmers(seq, k):
    for kmer in stream_kmers(seq, k):
        rmer = get_revcomp(kmer)
        yield min(rmer, kmer)

if __name__ == '__main__':
    pass