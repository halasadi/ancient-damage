from pyfaidx import Fasta
import csv
import pysam
import argparse as arg


MIN_BQ_SCORE = 20
MIN_MQ_SCORE = 30

def find_substitutions_wtags(aligned_pairs):
    """
    For every mutation, the mutation position
    along the read (mpos) and the position of the mutation in the genome (posg)
    """
    mpos = ()
    posg = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                mpos = mpos + (aligned_pairs[i][0],)
                posg = posg + (aligned_pairs[i][1],)
                
    return((mpos, posg))

def find_substitutions(read, chr, reference):
    ref_pos = read.get_reference_positions(full_length=True)
    align_qualities = read.query_qualities
    seq = read.query_sequence

    posg = ()
    mpos = ()
    for i in range(len(ref_pos)):
        if (ref_pos[i] is None or seq[i] is None):
            continue
        if (align_qualities[i] is None or align_qualities[i] < MIN_BQ_SCORE):
            continue
        if (i < read.qstart or i > read.qend):
            continue
        ref = reference[(chr-1)][ref_pos[i]]
        bp = seq[i]
        if (ref != bp):
            posg = posg + (ref_pos[i],)
            mpos = mpos + (i,)
    return((mpos, posg))


def get_bases_strandbreaks(read, chr, reference):
    ref_pos = read.get_reference_positions()
    leftbreak = reference[(chr-1)][ref_pos[0]-1]
    rightbreak = reference[(chr-1)][ref_pos[-1]+1]
    return((leftbreak, rightbreak))

if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    parser.add_argument("-f", "--fasta", required=True, help = "reference file")
    parser.add_argument("-o", "--out", required=True, help="out file")
    parser.add_argument("--add-chr", help = "add chr prefix?, you can find out by running samtools idxstats <your bamfile> | head -1 ", default = False, action = 'store_true', dest = "add_chr")
    parser.add_argument("--use-tags", help = "use the MD and NM tags in your bam file?", default = True, action = 'store_true', dest = "use_tags")
    
    args = parser.parse_args()

    ## ../data/T004_all_chr.bam
    samfile = pysam.AlignmentFile(args.bam, "rb")
    ## "/project/jnovembre/data/external_public/reference_genomes/hs37d5.fa"
    fastafile = Fasta(args.fasta, as_raw = True)
    patternsDict = {}

    chrs = [i for i in range(1,23)]
    
    for chr in chrs:
        
        for read in samfile.fetch(('chr' + str(chr)) if args.add_chr else str(chr)):
            if (read.mapping_quality < MIN_MQ_SCORE or read.is_duplicate):
                continue

            if (args.use_tags and read.get_tag('NM') == 0):
                continue

            seq = read.query_sequence
            
            if (args.use_tags):
                aligned_pairs = read.get_aligned_pairs(with_seq=True)
                (mutPos, posg) = find_substitutions_wtags(aligned_pairs)
            else:
                (mutPos, posg) = find_substitutions(read, chr, fastafile)
            
            mapq = read.query_qualities

            if (read.is_reverse):
                strando = '-'
            else:
                strando = '+'

            try:
                (leftbreak, rightbreak) = get_bases_strandbreaks(read, chr, fastafile)
            except IndexError:
                leftbreak = 'N'
                rightbreak = 'N'
                
            for i in range(len(posg)):
                pos = mutPos[i]
                mut = seq[pos]
        
                if (mapq[pos] < MIN_BQ_SCORE):
                    continue

                # we don't count mutation in soft clipped areas
                if (pos < read.qstart or pos > read.qend or mut == 'N'):
                    continue

                start = posg[i]-1
                end = posg[i]+1
                ref = fastafile[(chr-1)][start:(end+1)]
                patt = ref[0:2] + '->' + mut + ref[2:4]

                mutStart = pos - read.qstart
                # read.qend is not 0-based
                mutEnd = (read.qend-1) - pos

                val = (patt, mutStart, mutEnd, leftbreak, rightbreak, strando)
                if val in patternsDict:
                    patternsDict[val] += 1
                else:
                    patternsDict[val] = 1

    # write to file
    with open(args.out + '.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key[0], key[1], key[2], key[3], key[4], key[5], value])



