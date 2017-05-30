How to make the csv file to input in R

(1) index the file with samtools

`samtools index example_ancient_WG.bam`

(2) Run the python scipt with default options

`python generate_summary_bams.py -b example_ancient_WG.bam -f /path/to/reference/hs37d5.fa -o example_ancient_WG.csv --add-chr


To run view all the options, please type

`python --help generate_summary_bams.py`
