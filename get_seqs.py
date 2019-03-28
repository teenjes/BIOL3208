
"""
The goal for this program is to provide an iterable function that can be run
across a reference genome to isolate and return sequences that lie between two
primer sequences provided, given certain variables.
The function takes the following files:

    The reference genome path string 
        Either as a fasta or fastq file ending in that file name

    The primers file that contains both the forward and reverse primer sequences

    Parameters for determining whether the region selected is ITS or TEF
        Implemented as a maximum and minimum length of the extracted sequence
            This allows for extraction of other regions

"""

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import subprocess
import os
from shutil import copyfile
import argparse


parser = argparse.ArgumentParser(description="""
The goal for this program is to provide an iterable function that can be run
across a reference genome to isolate and return sequences that lie between two
primer sequences provided, given certain variables.
The function takes the following files:

    The reference genome path string 
        Either as a fasta or fastq file ending in that file name

    The primers file that contains both the forward and reverse primer sequences

    Parameters for determining whether the region selected is ITS or TEF
        Implemented as a maximum and minimum length of the extracted sequence
            This allows for extraction of other regions

""")
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true")
group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("reference_genome", help="The reference genome to scan for the particular sequence")
parser.add_argument("primer", help="The primer sequences on either side of the desired sequence")
parser.add_argument("minimum", type=int, help="Minimum length of the sequence")
parser.add_argument("maximum", type=int, help="Maximum length of the sequence")
args = parser.parse_args()

print(args.reference_genome)

if args.reference_genome[-5:] == "fastq":
    reference = SeqIO.convert((args.reference_genome), "fastq", args.reference_genome[:-5]+'fasta', "fasta")
elif args.reference_genome[-5:] == "fasta":
    reference = args.reference_genome
else:
    print("ERROR: File not supported. Please use a fasta or fastq file with a file ending of either fasta or fastq respectively")
    exit()

if args.primer[-5:] == "fasta":
    pass
else:
    print("ERROR: File not supported. Please use a fasta file with a .fasta ending")
    exit()


new_reference_file = "references/"+reference[11:-16]+"/"+reference[11:]

cmd0 = 'copyfile(reference, new_reference_file)'
subprocess.getoutput(cmd0)

database = "references/"+reference[11:-16]+"/"+reference[11:-6]+'_db'
outfmt6 = "references/"+reference[11:-16]+"/"+reference[11:-6]+'.outfmt6'

cmd = 'makeblastdb -in %s -dbtype nucl -out %s' % (reference, database)
subprocess.getoutput(cmd)

cmd2 = 'blastn -query %s -db %s -evalue=100000 -task "blastn-short" -outfmt 6 > %s' % (args.primer, database, outfmt6)
subprocess.getoutput(cmd2)

df = pd.read_csv(outfmt6, sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

forward_df = df.loc[(df['qseqid']=='forward')]
forward_df = forward_df.loc[df['evalue']==min(forward_df['evalue'])]
forward_df = forward_df.reset_index()

reverse_df = df.loc[(df['qseqid']=='reverse')]
reverse_df = reverse_df.loc[df['evalue']==min(reverse_df['evalue'])]
reverse_df = reverse_df.reset_index()

forward_bed = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd'])
for i in range(0, len(forward_df)):
    if forward_df['sstart'][i] <= forward_df['send'][i]:
        forward_bed = forward_bed.append({'chrom': forward_df['sseqid'][i],'chromStart': (forward_df['sstart'][i])-1, 'chromEnd': forward_df['send'][i]}, ignore_index=True)
    else:
        forward_bed = forward_bed.append({'chrom': forward_df['sseqid'][i],'chromStart': (forward_df['send'][i])-1, 'chromEnd': forward_df['sstart'][i]}, ignore_index=True)
forward_bed.sort_values(['chromStart'])

reverse_bed = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd'])
for i in reverse_df.index:
    if reverse_df['sstart'][i] < reverse_df['send'][i]:
        reverse_bed = reverse_bed.append({'chrom': reverse_df['sseqid'][i],'chromStart': (reverse_df['sstart'][i])-1, 'chromEnd': reverse_df['send'][i]}, ignore_index=True)
    else:
        reverse_bed = reverse_bed.append({'chrom': reverse_df['sseqid'][i],'chromStart': (reverse_df['send'][i])-1, 'chromEnd': reverse_df['sstart'][i]}, ignore_index=True)
reverse_bed.sort_values(['chromStart'])

intervals_list = []
for row in forward_bed.itertuples(index=True, name='Pandas'):
    counter = 0
    for rows in reverse_bed.itertuples(index=True, name='Pandas'):
        # Here lies the difference between sequences other than primers used - expected sequence length
        if getattr(row, 'chrom') == getattr(rows, 'chrom') and args.minimum < np.absolute(getattr(row, 'chromEnd') - getattr(rows, 'chromStart')) < args.maximum and getattr(row, 'chromEnd') > getattr(rows, 'chromStart'):
            counter += 1
            intervals_list.append([forward_bed['chrom'][getattr(row, 'Index')], reverse_bed['chromStart'][getattr(rows, 'Index')], forward_bed['chromEnd'][getattr(row, 'Index')]])
        elif getattr(row, 'chrom') == getattr(rows, 'chrom') and args.minimum < np.absolute(getattr(row, 'chromEnd') - getattr(rows, 'chromStart')) < args.maximum and getattr(row, 'chromEnd') < getattr(rows, 'chromStart'):
            counter += 1
            intervals_list.append([forward_bed['chrom'][getattr(row, 'Index')], forward_bed['chromEnd'][getattr(row, 'Index')], reverse_bed['chromStart'][getattr(rows, 'Index')]])
    if counter > 1:
        print("WARNING: two reverse primers identified for one forward primer")
        print(row)

intervals_frame = pd.DataFrame(data=intervals_list, columns=['chrom', 'chromStart', 'chromEnd'])

bedfile = "references/"+reference[11:-16]+"/"+reference[11:-6]+'.ITS.bedfile'

intervals_frame.to_csv(bedfile, sep='\t', header=False, index=False)

bedoutput = "references/"+reference[11:-16]+"/"+reference[11:-6]+'.ITS.bedoutput.fasta'
cmd3 = 'bedtools getfasta -fo %s -fi %s -bed %s' % (bedoutput, reference, bedfile)
subprocess.getoutput(cmd3)



visual = SeqIO.parse(bedoutput, "fasta")

    
    
if args.quiet:
    print("\n%i sequences found\n" % (len(intervals_frame)))
    counter = 1
    for record in visual:
        if "N" not in record.seq:
            print("\nSequence %i\n" % counter)
            counter += 1
        else:
            print("\nSequence %i\n" % counter)
            print("ERROR: sequence contains 'N' ")
            counter += 1
elif args.verbose:
    print("\n%i sequences found, regions shown below\n" % (len(intervals_frame)))
    print(intervals_frame)
    print("\nSequences shown below")
    counter = 1
    for record in visual:
        if "N" not in record.seq:
            print("\nSequence %i\n" % counter)
            print(record.seq)
            counter += 1
            print("\n")
        else:
            print("\nSequence %i\n" % counter)
            print("ERROR: sequence contains 'N' ")
            counter += 1
            print("\n")
else:
    print("\n%i sequences found, sequences shown below" % (len(intervals_frame)))
    counter = 1
    for record in visual:
        if "N" not in record.seq:
            print("\nSequence %i\n" % counter)
            print(record.seq)
            counter += 1
            print("\n")
        else:
            print("\nSequence %i\n" % counter)
            print("ERROR: sequence contains 'N' ")
            counter += 1
            print("\n")
