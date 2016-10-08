import pysam
import csv
import argparse as arg

MIN_BQ_SCORE = 20

# returns the position of the mutation along the read (0 based)
def find_substitutions(aligned_pairs):
    inds = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                inds = inds + (i,)
    return(inds)

### "AAG->TCC" 
def concat_pattern(ind, start, end, aligned_pairs, seq, TYPE):
    ref = aligned_pairs[ind][2].upper()
    mut = seq[ind]
    if (TYPE == '0L'):
        patt = 'XX' + ref + '->' + mut + seq[ind+1] + seq[ind+2]
    elif (TYPE == '1L'):
        patt = 'X' + seq[ind-1] + ref + '->' + mut + seq[ind+1] + seq[ind+2]
    elif (TYPE == '0R'):
        patt = seq[ind-2] + seq[ind-1] + ref + '->' + mut + 'XX'
    elif (TYPE == '1R'):
        patt = seq[ind-2] + seq[ind-1] + ref + '->' + mut + seq[ind+1] + 'X'
    else:
        patt = seq[ind-2] + seq[ind-1] + ref + '->' + mut + seq[ind+1] + seq[ind+2]
    return(patt)

def is_quality(qs):
    if (length([x for x in qs if x < MIN_BQ_SCORE]) > 0):
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
        
        seq = read.seq
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        mutInds = find_substitutions(aligned_pairs)
        
        for ind in mutInds:
            
            mut = seq[aligned_pairs[ind][0]]
            TYPE = 'N'
            
            # we don't count mutation in soft clipped areas
            if (ind < read.qstart or ind > read.qend or mut == 'N'):
                continue

            # no base pair flanking to the left
            if (ind < (read.qstart+1)):
                TYPE = '0L'
                start = read.qstart
                end = read.qstart + 1
                
            # one base pair flanking to the left
            elif (ind < (read.qstart+2)):
                TYPE = '1L'
                start = read.qstart
                end = read.qstart + 2

            # no base pair flanking to the right
            elif (ind > (read.qend-2)):
                TYPE = '0R'
                start = read.qend 
                end = read.qend + 1

            # one base pair flanking to the right
            elif (ind > (read.qend-3)):
                TYPE = '1R'
                start = read.qend-1
                end = read.end + 1
            else:
                start = ind-2
                end = ind+3

            qualityScores = read.query_qualities 
            if (!is_quality(qualityScores[start:end])):
                continue

            patt = concat_pattern(ind, start, end, aligned_pairs, read.seq, TYPE)
            
            # check this line
            val = (patt, start, end)              
            if patt in patternsDict:
                patternsDict[patt] += 1
            else:
                patternsDict[patt] = 1

    # write to file
    with open(args.out, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in patternsDict.items():
            writer.writerow([key[0], key[1], key[2], value])
