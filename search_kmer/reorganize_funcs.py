#useful functions for reorganizing data in kmer analysis
#1) functions for reading and writing to k_data files 

import sys

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo
from tqdm import tqdm
import pandas as pd
import os

sys.path.append('../..')   

from kmer_tools.search_kmer.count_func_lib import *


seqs_folder     = "../download_genbank/"
spec_folder_normal     = "../../data/k_spectra/"
spec_folder_strand     = "../../data/k_spectra_strand/"
spec_folder = spec_folder_strand
k_data_folder_normal    = "../../data/k_data/"
k_data_folder_strand    = "../../data/k_data_strand/"
k_data_folder           = k_data_folder_normal
monomers = ["A","T","C","G"]


def k_words_expand(word_list):
    monomers = ["A","T","C","G"]
    new_word_list = []
    for word in word_list:
        for letter in monomers:
            new_word_list.append(word+letter)
    return new_word_list

def kmer_list(k_value):
    word_list = []
    monomers = ["A","T","C","G"]
    word_list = monomers
    for i in range(k_value-1):
        word_list = k_words_expand(word_list)
    return word_list

def get_spec_coverage(spec):
    k=spec.k
    kmers = [kmer for kmer in spec.get_kmers()]
    num_tot =spec.sum()
    return num_tot + k - 1

def read_kdat_chromSpec(kdata_folder,domain,func_name,k_value,count_func=word_num,chromosome = None):
    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,func_name,1,5), sep="\t",header =0)
    words = kmer_list(k_value)
    counts=err = 0
    for index, row in f.iterrows():
        if (str(chromosome) ==row["chromosome"]):
            counts,err= count_func(row,words)
            break
    return Oligo.Kmer.KmerSpectrum(kmers=words,values=counts)

def read_kdat_orgaRow(kdata_folder,domain,func_name,orga_name = None,verbose=False):
    if verbose:
        print("opening: "+kdata_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,6))
    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,6), sep="\t",header =0)
    for index, row in f.iterrows():
        if (orga_name ==row["organism"]):
            return row
    return None

"""
def read_kdat_dimerRow(kdata_folder,domain,func_name,orga_name = None,verbose=False,chromo=False):
    label = "organism"
    if chromo: label = "chromosome"
    if verbose:
        print("opening: "+k_data_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,2))
    f=pd.read_csv(k_data_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,2), sep="\t",header =0)
    for index, row in f.iterrows():
        if (orga_name ==row["organism"]):
            return row
    return None
"""


def read_kdat_dimerRow(kdata_folder,domain,func_name,orga_name = None,verbose=False,chromo=False):
    label = "organism"
    if chromo: label = "chromosome"
    if verbose:
        print("opening: "+k_data_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,2))
    f=pd.read_csv(k_data_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,2), sep="\t",header =0)
    if chromo: f=pd.read_csv(k_data_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,func_name,1,2), sep="\t",header =0)
    for index, row in f.iterrows():
        if (orga_name ==row[label]):
            return row
    return None

def read_kdat_orgaSpec(kdata_folder,domain,func_name,k_value,count_func=word_num,orga_name = None,verbose=False):
    if verbose:
        print("opening: "+kdata_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,6))
    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,6), sep="\t",header =0)
    words = kmer_list(k_value)
    counts=err = 0
    for index, row in f.iterrows():
        if (orga_name ==row["organism"]):
            counts,err= count_func(row,words)
            index=-1
            break
    if verbose:
        print(orga_name," index: ",index)
        print(words,counts)
    return Oligo.Kmer.KmerSpectrum(kmers=words,values=counts)

def read_kdat_orgaFunc(kdata_folder,domain,func_name,k_value,count_func=word_num,orga_name = None):
    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,1,5), sep="\t",header =0)
    words = kmer_list(k_value)
    for index, row in f.iterrows():
        if (orga_name ==row["organism"]):
            returnRow =row
            index=-1
            break
    #print(orga_name," index: ",index)
    counts,err = count_func(row,words)
    return counts,err

def read_kdat_spec_Opt(kdata_folder,domain,chromo,func_name,k_value,count_func=word_num): #optimised version (now outdated, replaced by organism version)
    #1st: locate row of wanted entry
    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,func_name,1,5), sep="\t",header =0,usecols=["chromosome"])
    for index, row in f.iterrows():
        if (str(chromo) ==row["chromosome"]):
            break
    #2nd load whole row of just the wanted chromosome
    if (index == 0): skip_inds = 0
    elif (index > 0) : skip_inds = np.arange(0+1,index+1)

    f=pd.read_csv(kdata_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,func_name,1,5), sep="\t",header =0,skiprows=skip_inds,nrows=1)
    words = kmer_list(k_value)
    counts=err = 0
    row = f.iloc[0]
    counts,err= count_func(row,words)
    return Oligo.Kmer.KmerSpectrum(kmers=words,values=counts)

def name_to_orga(orga_names,organisms):
    orga_list=[]
    for name in orga_names:
        organism = None
        for orga in organisms:
            if(orga.name == name):
                organism = orga
        orga_list.append(organism)
    return orga_list

def name_to_domain(orga_names,organisms):
    orga_list=[]
    for name in orga_names:
        organism = None
        for orga in organisms:
            if(orga.name == name):
                organism = orga
        orga_list.append(organism)
    return orga_list
            
def store_kmer_chromosomal(organisms,function_names,domain,k_values=[1,2,3]): 
    # 1st: read in k_mer data from files
    # 2nd: write into new combined file (line by line for each chromosome)
    for i,func_name in tqdm(enumerate(function_names)):
        print(func_name)
        entries = ["chromosome","length"]
        with open(k_data_folder+domain+'/temp.txt', 'w') as f:
            #1st: generate list of categories
            entries = ["chromosome","length","G+C"]
            k_words_all = []
            for k_value in k_values:
                k_words_all += kmer_list(k_value)
            entries += k_words_all

            #2nd: write header
            for label in entries:
                f.write(label+"\t")
            f.write("\n")
            for orga in tqdm(organisms):
                genome = Oligo.File.read_genome(orga.name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose = 0)
                if(len(genome)==0):
                    print("Error: orga not found")
                    return None
                #3rd: write one line for each chromosome (read all data -> write all data)
                for chromo in genome:
                    values = {}
                    values["chromosome"] = str(chromo)
                    values["length"] = chromo.length
                    for k_word in k_words_all: #initiializing dictionary with 0-entries
                        values[str(k_word)] = 0
                    for k in k_values:
                        #print("opening: "+(spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k)).replace('?','').replace(':','_'))
                        spectrum = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k)).replace(':','_').replace('?',"_"),verbose = 0)
                        k_mers_sorted   = sorted(spectrum.get_kmers())
                        bins_sorted     = spectrum.get_bins()
                        for i in range(len(k_mers_sorted)):
                            values[str(k_mers_sorted[i])] = bins_sorted[i].count
                    values["G+C"] = str(int(values["G"])+int(values["C"]))
                    #print(values)
                    for label in entries:
                        f.write(str(values[label])+"\t")
                    f.write("\n")
        #rename
        os.replace(k_data_folder+domain+'/temp.txt',k_data_folder+domain+'/kdata_%s_%s_k=%d-%d.k_dat' % (domain,func_name,k_values[0],k_values[-1]))

def store_kmer_orga(organisms,function_names,domain,k_values=[1,2,3],verbose=0): 
    # 1st: read in k_mer data from files
    # 2nd: write into new combined file (line by line for each chromosome)
    for i,func_name in tqdm(enumerate(function_names)):
        print(func_name)
        entries = ["chromosome","length"]
        with open(k_data_folder+domain+'/temp.txt', 'w') as f:
            #1st: generate list of categories
            entries = ["organism","length","G+C"]
            k_words_all = []
            for k_value in k_values:
                k_words_all += kmer_list(k_value)
            entries += k_words_all

            #2nd: write header
            for label in entries:
                f.write(label+"\t")
            f.write("\n")
            for orga in tqdm(organisms):
                genome = Oligo.File.read_genome(orga.name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose = 0)
                #3rd: write one line for each chromosome (read all data -> write all data)
                values = {}
                values["organism"] = orga.name
                for k_word in k_words_all: #initiializing dictionary with 0-entries
                    values[str(k_word)] = 0
                length_arr=[]
                spectra_chrom=[]
                for k in k_values:
                    length = 0
                    for chromo in genome:
                        if verbose: print("opening: "+(spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k)).replace('?','').replace(':','_'))
                        spectra_chrom.append(Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k)).replace('?','_').replace(':','_'),verbose = 0))
                        length += chromo.length
                    spectrum = Oligo.Kmer.combine(spectra_chrom,verbose=0)
                    k_mers_sorted   = sorted(spectrum.get_kmers())
                    bins_sorted     = spectrum.get_bins()
                    for i in range(len(k_mers_sorted)):
                        values[str(k_mers_sorted[i])] = bins_sorted[i].count
                values["length"] = length
                values["G+C"] = str(int(values["G"])+int(values["C"]))
                #print(values)
                for label in entries:
                    f.write(str(values[label])+"\t")
                f.write("\n")
        os.replace(k_data_folder+domain+'/temp.txt',k_data_folder+domain+'/kdata_%s_orga_%s_k=%d-%d.k_dat' % (domain,func_name,k_values[0],k_values[-1]))