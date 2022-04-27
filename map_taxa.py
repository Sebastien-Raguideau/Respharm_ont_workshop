#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from os.path import dirname,basename
from collections import defaultdict,Counter
from ete3 import NCBITaxa
import numpy as np
import pysam
import glob
import os
import argparse


def matrix_write(matrix,file_name,col_names,row_names) :
    with open(file_name,"w") as handle:
        handle.write("/\t%s\n"%"\t".join(col_names))
        handle.writelines('%s\t%s\n'%(row_names[index],"\t".join(list(map(str,line)))) for index,line in enumerate(matrix))


def filter_minimap2_bam(bam_file, pid_min, breadth_min):
    contig_to_map = defaultdict(list)
    contig_to_map_pid = defaultdict()
    contig_to_map_bid = defaultdict()
    samfile = pysam.AlignmentFile(bam_file, "rb")
    database = samfile.references
    # check mapping id is at 99%
    # check that length of alignment it at least 10% of contig
    fail = defaultdict(list)
    for refmap in samfile:

        fungus = database[refmap.reference_id]
        fungus_len = refmap.reference_length

        contig = refmap.qname
        contig_len = refmap.qlen
        if contig_len==0:
            continue

        alen = refmap.alen
        bcv = alen/fungus_len
        matchs = sum(el[1] for el in refmap.cigartuples if el[0]==0)
        pid = matchs/float(alen)
        breath_cov = alen/float(contig_len)
        if (pid>=pid_min)&(breath_cov>=breadth_min):
            contig_to_map[contig].append([fungus,fungus_len,contig_len,alen,matchs,pid,bcv,refmap.qstart,refmap.qend,refmap.reference_start,refmap.reference_end])
        else:
            fail[contig].append([fungus,fungus_len,contig_len,alen,matchs,pid,bcv,refmap.qstart,refmap.qend,refmap.reference_start,refmap.reference_end])
    return contig_to_map,fail

def get_consistent_taxid(tax):
    ncbi =  NCBITaxa()
    all_taxid = {el[0] for el in ncbi.get_name_translator(tax).values()}
    i = 1
    candidate = [el[0] for el in ncbi.get_name_translator([tax[-i]]).values()]
    while not candidate:
        i+=1
        candidate = [el[0] for el in ncbi.get_name_translator([tax[-i]]).values()]
    taxa_id = max(candidate,key=lambda x:len(set(ncbi.get_lineage(x))&all_taxid))
    ranks_dict = defaultdict(str)
    ranks_dict.update({rank:list(ncbi.get_taxid_translator([taxid]).values())[0] for taxid,rank in ncbi.get_rank(ncbi.get_lineage(taxa_id)).items()})
    ranks = [ranks_dict[rank] for rank in ['superkingdom',"kingdom","phylum", "class", "order","family", "genus", "species"]]
    return ranks

def main(SILVA,input_folder,out,pid_min,breadth_min):
    SAMPLES = {basename(file).split("_mapped_sorted")[0]:file for file in glob.glob("%s/*.bam"%input_folder)}

    header_to_taxa = {header.rstrip().split(" ")[0]:header.rstrip().split(" ")[1].split(";") for header,seq in sfp(open(SILVA))}

    # sample, reads_tot, reads 95%
    stats = []
    # load contig infos
    sorted_samples = sorted(SAMPLES)
    counts = defaultdict(lambda:np.zeros(len(SAMPLES)))
    for sample,bam_file in SAMPLES.items():
        sample_index = sorted_samples.index(sample) 
        # get breadth of coverage info 
        pid_min = 0.95
        breadth_min = 0.50
        read_to_map,fail = filter_minimap2_bam(bam_file,pid_min,breadth_min)
        stats.append([sample,len(set(read_to_map.keys())|set(read_to_map.keys())),len(set(read_to_map.keys()))])

        # just do some best hit
        read_to_best = {key:max(val,key=lambda x:x[5]*x[6]) for key,val in read_to_map.items()}

        # just do some best hit
        for contig,hit in read_to_best.items():
            counts[hit[0]][sample_index]+=1
    sorted_otu = sorted(counts.keys())
    mat = np.zeros((len(sorted_otu),len(SAMPLES)))
    for index,otu in enumerate(sorted_otu):
       mat[index,:]+=counts[otu]
    
    # output cnt matrix
    output = "%s/otu_cnts.tsv"%OUT
    matrix_write(mat.astype(int),output,sorted_samples,sorted_otu)

    # output taxa table
    header = ['superkingdom',"kingdom","phylum", "class", "order","family", "genus", "species"]
    output = "%s/otu_taxa.tsv"%OUT
    with open(output,"w") as handle:
        handle.write("otu\t%s\n"%("\t".join(header)))
        for otu in sorted_otu:
            taxa = get_consistent_taxid(header_to_taxa[otu])
            handle.write("%s\t%s\n"%(otu,"\t".join(taxa)))




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="folder containing a unique .bam file per sample")
    parser.add_argument("silva_db", help="path to silva database")
    parser.add_argument("output", help="output folder")
    parser.add_argument("-i",default="0.95", help="minimum percent identity")
    parser.add_argument("-b",default="0.50", help="minimum breadth of coverage")
    args = parser.parse_args()
    main(args.silva_db,args.folder,args.output,float(arg.i),float(arg.b))
