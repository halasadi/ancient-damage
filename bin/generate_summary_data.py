import pysam
import csv
import argparse as arg

MIN_BQ_SCORE = 20
MIN_MP_SCORE = 30

def find_substitutions(aligned_pairs):
    refs = ()
    mpos = ()
    for i in range(0, len(aligned_pairs)):
        if (aligned_pairs[i][2] is not None):
            if (aligned_pairs[i][2].islower()):
                mpos = mpos + (aligned_pairs[i][0],)
                refs = refs + (aligned_pairs[i][2].upper(),)
    return((refs, mpos))

### "AAG->TCC" 
def concat_pattern(pos, ref, seq, TYPE):
    mut = seq[pos]
    if (TYPE == '0L'):
        patt = 'XX' + ref + '->' + mut + seq[pos+1] + seq[pos+2]
    elif (TYPE == '1L'):
        patt = 'X' + seq[pos-1] + ref + '->' + mut + seq[pos+1] + seq[pos+2]
    elif (TYPE == '0R'):
        patt = seq[pos-2] + seq[pos-1] + ref + '->' + mut + 'XX'
    elif (TYPE == '1R'):
        patt = seq[pos-2] + seq[pos-1] + ref + '->' + mut + seq[pos+1] + 'X'
    else:
        patt = seq[pos-2] + seq[pos-1] + ref + '->' + mut + seq[pos+1] + seq[pos+2]
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

        if (read.get_tag('NM') == 0 or read.mapping_quality < MIN_MQ_SCORE):
            continue

        seq = read.query_sequence
        aligned_pairs = read.get_aligned_pairs(with_seq=True)

        (refs, mutPos) = find_substitutions(aligned_pairs)

        for i in range(len(mutPos)):

            pos = mutPos[i]
            ref = refs[i]
            mut = seq[pos]
            TYPE = 'N'
            
            # we don't count mutation in soft clipped areas
            if (pos < read.qstart or pos > read.qend or mut == 'N'):
                continue

            # no base pair flanking to the left
            if (pos < (read.qstart+1)):
                TYPE = '0L'
                start = read.qstart
                end = read.qstart + 1
                
            # one base pair flanking to the left
            elif (pos < (read.qstart+2)):
                TYPE = '1L'
                start = read.qstart
                end = read.qstart + 2

            # qend is not 0 based (it is the length of the read)
            # no base pair flanking to the right
            elif (pos > (read.qend-2)):
                TYPE = '0R'
                start = read.qend - 1
                end = read.qend 

            # one base pair flanking to the right
            elif (pos > (read.qend-3)):
                TYPE = '1R'
                start = read.qend-2
                end = read.qend
            else:
                start = pos-2
                end = pos+3

            qualityScores = read.query_qualities 
            if (not is_quality(qualityScores[start:end])):
                continue

            mutStart = pos - read.qstart
            # read.qend is not 0-based
            mutEnd   = (read.qend-1) - pos

            patt = concat_pattern(pos, ref, seq, TYPE)
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
