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
    fasta = pysam.FastaFile(args.fastafile)
    
    nFlanking = int(args.nflanking)

    ### for now use a dictionary
    ### "AAG->TCC" is a key
    ### and the amount of such occurences is the value
    patternsDict = {}

    chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
        'chr19', 'chr20', 'chr21', 'chr22']

    for chr in chrs:
        for read in samfile.fetch(chr):
            NM = read.get_tag('NM')
            if (NM == 0):
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
                if (mut == 'N' or quality_scores[ind] < MIN_BQ_SCORE):
                    continue
                pos = aligned_pairs[ind][1]
                leftFlank = fasta.fetch(chr, pos-nFlanking, pos)
                ref = fasta.fetch(chr, pos, pos+1)
                rightFlank = fasta.fetch(chr, pos+1, pos+nFlanking+1)
                patt = leftFlank.upper() + ref.upper() + '->' + mut + rightFlank.upper()
                if patt in patternsDict:
                    patternsDict[patt] += 1
                else:
                    patternsDict[patt] = 1

    # write to file
    with open(args.out + '.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key, value])
