import pysam
import csv
import argparse as arg

MIN_BQ_SCORE = 20


# re-write code so instead of find_substitutions function, replace it with a function that finds positions of all mutations, reference, and mutation

# returns the position of the mutation along the read (0 based)
def find_substitutions(aligned_pairs):
    inds = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                inds = inds + (i,)
    return(inds)

### "AAG->TCC" 
def concat_pattern(ind, aligned_pairs, seq, TYPE):
    ref = aligned_pairs[ind][2].upper()
    mutPos = aligned_pairs[ind][0]
    mut = seq[mutPos]
    if (TYPE == '0L'):
        patt = 'XX' + ref + '->' + mut + seq[mutPos+1] + seq[mutPos+2]
    elif (TYPE == '1L'):
        patt = 'X' + seq[mutPos-1] + ref + '->' + mut + seq[mutPos+1] + seq[mutPos+2]
    elif (TYPE == '0R'):
        patt = seq[mutPos-2] + seq[mutPos-1] + ref + '->' + mut + 'XX'
    elif (TYPE == '1R'):
        patt = seq[mutPos-2] + seq[mutPos-1] + ref + '->' + mut + seq[mutPos+1] + 'X'
    else:
        patt = seq[mutPos-2] + seq[mutPos-1] + ref + '->' + mut + seq[mutPos+1] + seq[mutPos+2]
    return(patt)

def is_quality(qs):
    if (len([x for x in qs if x < MIN_BQ_SCORE]) > 0):
        return(False)
    return(True)

if __name__ == '__main__':
    parser = arg.ArgumentParser()
    parser.add_argument("-f", "--fname", required=True, help="bam file")
    parser.add_argument("-o", "--out", required=True, help="out file")
    args = parser.parse_args()
    
    samfile = pysam.AlignmentFile(args.fname, "rb") 

    patternsDict = {}

    for read in samfile.fetch():

        if (read.get_tag('NM') == 0 or read.mapping_quality < 30):
            continue
        
        seq = read.query_sequence
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        mutInds = find_substitutions(aligned_pairs)

        #ind must be less than read.qend
        
        for ind in mutInds:
            
            mutPos = aligned_pairs[ind][0]            
            mut = seq[mutPos]
            TYPE = 'N'
            
            # we don't count mutation in soft clipped areas
            if (mutPos < read.qstart or mutPos > read.qend or mut == 'N'):
                continue

            # below indices includes the mutation
            # to check if all the base quality scores are greater than 30

            # no base pair flanking to the left
            if (mutPos < (read.qstart+1)):
                TYPE = '0L'
                start = read.qstart
                end = read.qstart + 1
                
            # one base pair flanking to the left
            elif (mutPos < (read.qstart+2)):
                TYPE = '1L'
                start = read.qstart
                end = read.qstart + 2

            # qend is not 0 based (it is the length of the read)
            # no base pair flanking to the right
            elif (mutPos > (read.qend-2)):
                TYPE = '0R'
                start = read.qend - 1
                end = read.qend 

            # one base pair flanking to the right
            elif (mutPos > (read.qend-3)):
                TYPE = '1R'
                start = read.qend-2
                end = read.qend
            else:
                start = ind-2
                end = ind+3

            qualityScores = read.query_qualities 
            if (not is_quality(qualityScores[start:end])):
                continue

            mutStart = mutPos - read.qstart
            mutEnd   = (read.qend-1) - mutPos

            patt = concat_pattern(ind, aligned_pairs, seq, TYPE)
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
