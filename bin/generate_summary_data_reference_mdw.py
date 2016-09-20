import pysam
import csv
import argparse as arg

# returns index of aligned_pairs that contains the substitution (0 based)
def find_substitutions(aligned_pairs):
    inds = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                inds = inds + (i,)
    return(inds)


MIN_BQ_SCORE = 20

if __name__ == '__main__':

    parser = arg.ArgumentParser()
    parser.add_argument("-bf", "--bamfile", required=True, help="bam file")
    parser.add_argument("-ff", "--fastafile", required=True, help="fasta file")
    parser.add_argument("-n", "--nflanking", required=True, help="number of flanking base pairs")
    parser.add_argument("-o", "--out", required=True, help="out file")
    args = parser.parse_args()
    
    samfile = pysam.AlignmentFile(args.bamfile, "rb") 
    fastafile = pysam.FastaFile(args.fastafile)
    
    nFlanking = int(args.nflanking)

    ### for now use a dictionary
    ### "AAG->TCC" is a key
    ### and the amount of such occurences is the value
    patternsDict = {}


    chrs = [str(i) for i in range(1, 23, 1)]
    for chr in chrs:
        for read in samfile.fetch(chr):
            NM = read.get_tag('NM')
            if (NM == 0 or (read.mapping_quality < 30)):
                continue
            seq = read.seq
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            inds = find_substitutions(aligned_pairs)
            for ind in inds:
                mut = seq[aligned_pairs[ind][0]]
                quality_scores = read.query_qualities
                # seems to be some sort of bug with pysam
                # need to check if they are of the same length
                if (len(quality_scores) != len(aligned_pairs)):
                    continue

                # throw away mutation on the first and last two base pairs
                if (mut == 'N' or quality_scores[ind] < MIN_BQ_SCORE or ind < 2 or ind > (len(seq)-3) or ind < read.qstart or ind > read.qend):
                    continue
                pos = aligned_pairs[ind][1]
                
                ref = fastafile.fetch('chr' + chr, pos-nFlanking, pos+nFlanking+1)
                ref = ref.upper()
                patt = ref[:(nFlanking+1)] + '->' + mut + ref[-nFlanking:]
                if patt in patternsDict:
                    patternsDict[patt] += 1
                else:
                    patternsDict[patt] = 1

    # write to file
    with open(args.out, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key, value])
