from pyfaidx import Fasta
import csv
import pysam
import argparse as arg


MIN_BQ_SCORE = 20
MIN_MQ_SCORE = 30

def find_substitutions(aligned_pairs):
    refs = ()
    mpos = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                mpos = mpos + (aligned_pairs[i][0],)
                refs = refs + (aligned_pairs[i][2].upper(),)
    return((refs, mpos))


if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    parser.add_argument("-f", "--fasta", required=True, help = "reference file")
    parser.add_argument("-o", "--out", required=True, help="out file")
    parser.add_argument("--add-chr", help = "add chr prefix?, you can find out by running samtools idxstats <your bamfile> | head -1 ", default = True, action = 'store_true')
    
    
    args = parser.parse_args()

    print(args.addChr)
    
    samfile = pysam.AlignmentFile(args.fname, "rb")
    ## "/project/jnovembre/data/external_public/reference_genomes/hs37d5.fa"
    fastfile = Fasta(args.fasta, as_raw = True)

    patternsDict = {}

    chrs = [str(i) for i in range(1,23)]
    
    for chr in chrs:
        
        for read in samfile.fetch(('chr' + chr) if parser.add_chr else chr):
            if (read.get_tag('NM') == 0 or read.mapping_quality < MIN_MQ_SCORE):
                continue

            seq = read.query_sequence
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            (refs, mutPos) = find_substitutions(aligned_pairs)
            mapq = read.query_qualities
    
            for i in range(len(mutPos)):

                pos = mutPos[i]
                mut = seq[pos]
        
                if (mapq[pos] < MIN_MP_SCORE):
                    continue

                # we don't count mutation in soft clipped areas
                if (pos < read.qstart or pos > read.qend or mut == 'N'):
                    continue

                start = pos-2
                end = pos+2
                # 0-based
                ref = fastafile[(chr-1)][(start-1):end]
                patt = ref[0:2] + mut + ref[3:5]

                mutStart = pos - read.qstart
                # read.qend is not 0-based
                mutEnd = (read.qend-1) - pos
                val = (patt, mutStart, mutEnd)

                if val in patternsDict:
                    patternsDict[val] += 1
                else:
                    patternsDict[val] = 1

    # write to file
    with open(args.out, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key[0], key[1], key[2], value])


