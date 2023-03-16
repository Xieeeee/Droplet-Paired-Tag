#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='filter bam based on barcodes field')
parser.add_argument('--cluster', type=str, dest="cluster", help='four columns meta: sample.bam, library, barcode, cluster')
parser.add_argument('--indir', type=str, dest="indir", help='directory path of input bam')
parser.add_argument('--field', type=str, dest="fild", default="CB", help='bam field storing barcodes. Will add library prefix post spliting.')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output directory')

args = parser.parse_args()

import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc

def run():
    """Split bam and add initial"""
    start_time = pc()
    """ init input files """
    clusterf = args.cluster
    dirPrefix = str(args.indir)
    outPrefix = args.outprfx
    fild = args.fild
    print("filtering bam files by clusters...")
    generate_bams(clusterf, dirPrefix, outPrefix, fild)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bams(clusterf, dirPrefix, outPrefix, field):
    clust_dat = pd.read_csv(clusterf, sep="\t") ### four columns meta: sample.bam, library, barcode, cluster
    clust_dat["uniq_barcode"] = clust_dat[['library', 'barcode']].apply(lambda x: ':'.join(x), axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values, index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
    sample_id = pd.unique(clust_dat["sample"])
    templateBamF = pysam.AlignmentFile("".join((dirPrefix, "/", sample_id[0])))
    for cls_id in pd.unique(clust_dat["cluster"]):
        outsamfname = "".join((outPrefix, ".cluster.", str(cls_id), ".sam"))
        opysam = pysam.AlignmentFile(outsamfname, "wh", template=templateBamF)
        outf_dict[cls_id] = opysam
    if len(sample_id) == 1:
        sample = sample_id[0]
        library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
        library = library_id[0]
        generate_bam_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, field)
    else:
        for sample in sample_id:
            library_id = pd.unique(clust_dat.loc[clust_dat['sample'] == sample]['library'])
            library = library_id[0]
            generate_bam_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, field)

def generate_bam_worker(clust_dat_dict, outf_dict, sample, library, dirPrefix, outPrefix, field):
    print("extract reads from", sample)
    bamfname = "".join((dirPrefix, "/", sample))
    bamF = pysam.AlignmentFile(bamfname, "rb")
    for b in bamF.fetch(until_eof=True):
        try:
            old_barcode = b.get_tag(field, with_value_type=False)
        except KeyError:
            continue
        ## old_barcode = b.query_name.split(':')[0]
        new_barcode = ":".join((library, old_barcode))
        if new_barcode in clust_dat_dict:
            # rest_qname = ":".join(b.query_name.split(':')[1:])
            # new_qname = ":".join((str(new_barcode),rest_qname))
            # b.query_name = str(new_qname)
            b.set_tag(field, new_barcode, replace=True)
            cls = clust_dat_dict[new_barcode]
            osam = outf_dict[cls]
            osam.write(b)
        else:
            continue
    osam.close()
    bamF.close()

if __name__ == "__main__":
    """filter bam based on QNAMES"""
    run()


