import pysam
import csv
import argparse as arg

def find_substitutions(aligned_pairs):
    inds = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                inds = inds + (i,)
    return(inds)


### "AAG->TCC" 
def get_pattern(ind, aligned_pairs, seq):
    ref = aligned_pairs[ind][2].upper()
    mut = seq[ind]
    patt = seq[ind-2] + seq[ind-1] + ref + '->' + mut + seq[ind+1] + seq[ind+2]
    return(patt)


MIN_BQ_SCORE = 20

if __name__ == '__main__':

    parser = arg.ArgumentParser()
    parser.add_argument("-f", "--fname", required=True, help="bam file")
    parser.add_argument("-o", "--out", required=True, help="out file")
    args = parser.parse_args()
    
    samfile = pysam.AlignmentFile(args.fname, "rb") 

    ### for now use a dictionary
    ### "AAG->TCC" is a key
    ### and the amount of such occurences is the value
    patternsDict = {}

    for read in samfile.fetch():
        NM = read.get_tag('NM')
        if (NM > 0):
            seq = read.seq
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            inds = find_substitutions(aligned_pairs)
            for ind in inds:
                # make sure that mutation has two flanking base pairs
                # and that it's a substitution
                mut = seq[aligned_pairs[ind][0]]
                if (ind < 2 or ind > (len(seq)-3) or mut == 'N'):
                    continue
                quality_scores = read.query_qualities
                if (quality_scores[ind-2] < MIN_BQ_SCORE or quality_scores[ind-1] < MIN_BQ_SCORE or quality_scores[ind] < MIN_BQ_SCORE or quality_scores[ind+1] < MIN_BQ_SCORE or quality_scores[ind+2] < MIN_BQ_SCORE):
                    continue
                patt = get_pattern(ind, aligned_pairs, seq)
                if patt in patternsDict:
                    patternsDict[patt] += 1
                else:
                    patternsDict[patt] = 1

    # write to file
    with open(args.out + '.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key, value])
