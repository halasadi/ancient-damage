from pyfaidx import Fasta
import csv
import pysam
import argparse as arg


MIN_BQ_SCORE = 20
MIN_MQ_SCORE = 30

def find_substitutions(read, chr, reference):
    ref_pos = read.get_reference_positions(full_length=True)
    align_qualities = read.query_qualities
    seq = read.query_sequence

    mut = ()
    posg = ()

    for i in range(len(ref_pos)):
        if (ref_pos[i] is None or seq[i] is None):
            continue
        if (aligned_qualities[i] is None or align_qualities[i] < MIN_BQ_SCORE):
            continue
        if (i < read.qstart or i > read.qend):
            continue
        ref = reference[(chr-1)][(ref_pos[i]-1)]
        seq.bp = seq[i]
        if (ref != seq.bp):
            mut  = mut + (seq.bp,)
            posg = posg + (ref,)
    return((mut, posg))

if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    parser.add_argument("-f", "--fasta", required=True, help = "reference file")
    parser.add_argument("-o", "--out", required=True, help="out file")
    parser.add_argument("--add-chr", help = "add chr prefix?, you can find out by running samtools idxstats <your bamfile> | head -1 ", default = False, action = 'store_true', dest = "add_chr")
    args = parser.parse_args()

    ## ../data/T004_all_chr.bam
    samfile = pysam.AlignmentFile(args.bam, "rb")
    ## "/project/jnovembre/data/external_public/reference_genomes/hs37d5.fa"
    fastafile = Fasta(args.fasta, as_raw = True)

    patternsDict = {}

    chrs = [i for i in range(1,23)]
    
    for chr in chrs:
        
        for read in samfile.fetch(('chr' + str(chr)) if args.add_chr else str(chr)):
            if (read.get_tag('NM') == 0 or read.mapping_quality < MIN_MQ_SCORE):
                continue

            for i in range(len(mutPos)):

            
                start = posg[i]-2
                end = posg[i]+2
                ref = fastafile[(chr-1)][start:(end+1)]
                patt = ref[0:3] + '->' + mut + ref[3:5]

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


