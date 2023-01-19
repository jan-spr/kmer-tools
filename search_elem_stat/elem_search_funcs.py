#Functions to investigate loci in different Domains, Organisms for varying genome elements/areas

#create files which store coverage, number, coverage density, num density for each element type on each chromosome
#need functions to :
#       a) convert different values
#       b) merge values over chromosomes
#       c) read and write values from file

#interesting categories:
# pc, npc, pseudo, integenic
# subgroups:
# intron, exon

#exotic interesting elements:
#a) RNAs: lncRNA, miRNA, snoRNA, snRNA, precursorRNA
#b) LTR, L1
#c) npc: tRNA,rRNA

#organisms:

import sys
from tabnanny import verbose
sys.path.insert(1, r'E:\Genomik\Oligo')                     #Windows Oligo
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')           #Linux Oligo
sys.path.append('../..')                                    #add parent folder

import Oligo
import kmer_tools
from tqdm import tqdm
import os
import numpy as np
import pandas as pd

seqs_folder = "../download_genbank/"
k_data_folder   = "../../data/k_data/"
elem_data_folder = "../../data/elem_data/"


#def store_coverage_chromosomal(organisms,functions,function_names,domain,):
    #for i,func in tqdm(enumerate(functions)):
        #f=pd.read_csv(k_data_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,function_names[i],1,3), sep="\t",header =0)
        #loci = func(chromosome=chromo)

def get_elem_props(chromosome,func):
    loci = func(chromosome=chromosome)
    length = chromosome.length
    coverage = Oligo.Loci.coverage_sum(loci)

    props = {}
    props["organism"]       = chromosome.orga.name
    props["chromosome"]     = chromosome.name
    props["chromo_str"]     = str(chromosome)
    props["id"]             = chromosome.id
    props["length"]         = length
    props["num"]            =  len(loci)
    props["coverage"]       = coverage
    props["num_density"]    = len(loci) / length
    props["coverage_density"] = coverage / length
    return props



def write_elem_file(orga_names,domain,function,func_name):
    with open(elem_data_folder+domain+'/temp.txt', 'w') as f:
        #1st create header
        entries = ["organism","chromosome","chromo_str","id","length","num","coverage","num_density","coverage_density"]

        #2nd write header
        for label in entries:
            f.write(label+"\t")
        f.write("\n")

        #3rd: write one line for each chromosome (create all data -> store in dict -> write into file)
        for name in orga_names:
            genome = Oligo.File.read_genome(name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose = 0)
            if(len(genome)==0):
                print("Error: orga not found")
                return None
            for j,chromo in enumerate(genome):
                #f=pd.read_csv(k_data_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,function_names[i],1,3), sep="\t",header =0)
                props = get_elem_props (chromo, function)
                for label in entries:
                    f.write(str(props[label])+"\t")
                f.write("\n")
    os.replace(elem_data_folder+domain+'/temp.txt',elem_data_folder+domain+'/elem_data_%s_%s.el_dat' % (domain,func_name))

def read_prop_chrom(chromosome,prop,domain,func_name):
    f=pd.read_csv(elem_data_folder+domain+'/elem_data_%s_%s.el_dat' % (domain,func_name), sep="\t",header =0)
    for index, row in f.iterrows():
        if (str(chromosome) ==row["chromosome"]):
            return row[prop]

def read_prop_orga(orga_name,prop,domain,func_name, weighted = True):
    f=pd.read_csv(elem_data_folder+domain+'/elem_data_%s_%s.el_dat' % (domain,func_name), sep="\t",header =0)
    length_arr = []
    prop_arr = []
    for index, row in f.iterrows():
        if (orga_name ==row["organism"]):
            length_arr.append(row["length"])
            prop_arr.append(row[prop])
    if weighted: 
        return np.average(prop_arr,weights=length_arr)
    else : 
        return np.average(prop_arr)

def func_loop_elem(orga_names,domain,functions,function_names):
    for i,func in tqdm(enumerate(functions)):
        write_elem_file(orga_names,domain,func,function_names[i])