#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/python

import pysam
import argparse
parser = argparse.ArgumentParser(description='median size of bam fragments?')
parser.add_argument('--bam', type=str, dest="bam", help='input bam')

args = parser.parse_args()
ibam = args.bam

length = []
bamF = pysam.AlignmentFile(ibam, 'rb')
for read in bamF.fetch(until_eof=True):
    if read.is_paired:
        if not read.is_unmapped:
            tt = read.template_length
            length.append(tt)
bamF.close()
length = [i for i in length if i > 0]
length.sort()
oput1 = sum(length) / len(length) ### mean
oput2 = length[round(len(length)/2)+1]

print("mean: ", round(oput1))
print("median: ", round(oput2))