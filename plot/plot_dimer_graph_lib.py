from cmath import nan
import sys
from turtle import color

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import seaborn as sns
from tqdm import tqdm
from kmer_tools.plot.dimer_func_lib import *
from kmer_tools.plot.plot_hm_lib import plot_heatmap,add_wrap_split

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data/"
dimer_bar_folder = "../../plots/dimer/dimer_bar/"
mean_dev_folder = "../../plots/dimer/dimer_bar/mean_dev/"
dimer_hm_folder = "../../plots/dimer/dimer_hm/"
dimer_graph_folder = "../../plots/dimer/dimer_graph/"

colors = [
    'red',
    'blue',
    'green',
    'orange',
    'yellow',
    'purple',
    'red',
    'blue'
]

linestyles = [
  "dashed",
  "dotted",
  "dashdot",
  "-.",
  "-",
  "--",
  ":",
  "solid",
]

def nan_filter(matrix):
    nan_arr = [np.isnan(row).any() for row in matrix]
    print(nan_arr)
    nan_ind = [i for i,bool in zip(np.arange(len(nan_arr)), nan_arr) if not bool]
    matrix_filtered = [spec for spec, bool in zip(matrix, nan_arr) if not bool]
    return nan_ind, matrix_filtered

def plot_dimer_cum_bia_OrgaGroups(organism_names, domain_list, group_names, orga_to_group_list, function_names, function_titles, colors=colors, file_label=""):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    for i,func_name in tqdm(enumerate(function_names)):
        mean_values = []
        err_values = []
        for group_name in group_names:
            print("group_name = ",group_name,":")
            #orga_names = []
            cum_vals = []
            group_orgas = [name for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            group_idxs = [i for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            for j,orga in zip(group_idxs,group_orgas):
                domain = domain_list[j]
                print(domain)
                func_spec_pair =[]
                dist = dimer_dist(orga,func_name,domain=domain)
                bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
                cum_val = np.sum(np.abs(bias_vals))
                cum_err =  np.sqrt(np.sum(bias_errs**2))
                cum_vals.append(cum_val)
            
        
            print("cum vals:",func_name," - ",cum_vals)
            cum_vals = [val for val in cum_vals if val == val]
            m = np.nanmean(cum_vals)
            e = np.nanstd(cum_vals)
            mean_values.append(m)
            err_values.append(e)
        print(func_name,mean_values)
        print(len(group_names),len(mean_values),len(err_values))
        print(group_names,mean_values,err_values)
        ax.errorbar(x=np.arange(len(group_names)), y=np.array(mean_values), yerr=np.array(err_values), label=function_titles[i], 
        fmt='o--',  linestyle=linestyles[i],alpha=0.6,capsize=5)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('cumulative Deviation (b.E.)')
    ax.set_xticks(np.arange(len(group_names)))
    ax.set_xticklabels(group_names)
    #ax.set_ylim((-1,1))
    plt.title(file_label + " cumulative dimer biases")
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "cumDimerBias_graph" + '_v1.pdf'
    else:
        filename = "cumDimerBias_graph" + '_v1.pdf'
    plt.savefig(dimer_graph_folder+filename)
    plt.savefig(dimer_graph_folder+filename[:-3]+"png")


def plot_dimer_TpA_CpG_scatter(organism_names, domains, domain_list, group_names, function_name, function_title, colors=colors, file_label=""):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    TA_ind = dimer_dist.dimer_words.index("TA")
    CG_ind = dimer_dist.dimer_words.index("CG")

    for i,domain in tqdm(enumerate(domains)):
        TA_vals = []
        TA_errs = []
        CG_vals = []
        CG_errs = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            dist = dimer_dist(orga,function_name,domain=domain)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            TA_vals.append(bias_vals[TA_ind])
            CG_vals.append(bias_vals[CG_ind])
            TA_errs.append(bias_errs[TA_ind])
            CG_errs.append(bias_errs[CG_ind])
        
        print(domain,zip(TA_vals, TA_errs, CG_vals, CG_errs))
        ax.scatter(x=TA_vals, y=CG_vals, label=domain, c=colors[i])
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('CpG bias')
    ax.set_xlabel('TpA bias')
    ax.axhline()
    ax.axvline()
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    plt.title(file_label + " CpG over TpA bias in "+function_title)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "TpA_CpG_corr_"+ function_name + '_v1'
    else:
        filename = "TpA_CpG_corr_"+ function_name+ '_v1'
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")

def plot_dimer_TpA_CpG_scatter_func(organism_names, domain, domain_list, group_names, function_names, function_titles, colors=colors, file_label=""):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    TA_ind = dimer_dist.dimer_words.index("TA")
    CG_ind = dimer_dist.dimer_words.index("CG")

    for i,function_name in tqdm(enumerate(function_names)):
        TA_vals = []
        TA_errs = []
        CG_vals = []
        CG_errs = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            dist = dimer_dist(orga,function_name,domain=domain)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            TA_vals.append(bias_vals[TA_ind])
            CG_vals.append(bias_vals[CG_ind])
            TA_errs.append(bias_errs[TA_ind])
            CG_errs.append(bias_errs[CG_ind])
        
        print(domain,zip(TA_vals, TA_errs, CG_vals, CG_errs))
        ax.errorbar(x=TA_vals,xerr=TA_errs, y=CG_vals,yerr=CG_errs, label=function_titles[i], c=colors[i],fmt='.',elinewidth=0.5,ecolor="gray",capsize=1)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('CpG bias')
    ax.set_xlabel('TpA bias')
    ax.axhline()
    ax.axvline()
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    ax.set_xlim((-0.8,0.2))
    plt.title(file_label + " CpG over TpA bias in "+domain)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "TpA_CpG_corr_"+ domain + '_v1'
    else:
        filename = "TpA_CpG_corr_"+ function_name+ '_v1'
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")

def plot_dimer_word_scatter_func(word1,word2,organism_names, domain, domain_list, group_names, function_names, function_titles, colors=colors, file_label=""):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    w1_ind = dimer_dist.dimer_words.index(word1)
    w2_ind = dimer_dist.dimer_words.index(word2)

    for i,function_name in tqdm(enumerate(function_names)):
        w1_vals = []
        w1_errs = []
        w2_vals = []
        w2_errs = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            dist = dimer_dist(orga,function_name,domain=domain,chromo=chromo)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            w1_vals.append(bias_vals[w1_ind])
            w2_vals.append(bias_vals[w2_ind])
            w1_errs.append(bias_errs[w1_ind])
            w2_errs.append(bias_errs[w2_ind])
        
        print(domain,zip(w1_vals, w1_errs, w2_vals, w2_errs))
        ax.scatter(x=w1_vals,y=w2_vals, label=function_titles[i], c=colors[i])
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel(word2[0]+'p'+word2[1]+' bias')
    ax.set_xlabel(word1[0]+'p'+word1[1]+' bias')
    ax.axhline()
    ax.axvline()
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    #ax.set_xlim((-0.8,0.2))
    plt.title(file_label +" "+ word2[0]+'p'+word2[1]+' over '+word1[0]+'p'+word1[1]+' bias in '+domain)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + word1+"_"+word2+"_corr_"+ domain + '_v1'
    else:
        filename = + word1+"_"+word2+"_corr_"+ function_name+ '_v1'
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")

def plot_dimer_word_scatter_func_chromo(word1,word2,organism_names, domain, domain_list, group_names, function_names, function_titles, colors=colors, file_label="",chromo=True):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    w1_ind = dimer_dist.dimer_words.index(word1)
    w2_ind = dimer_dist.dimer_words.index(word2)

    for i,function_name in tqdm(enumerate(function_names)):
        w1_vals = []
        w1_errs = []
        w2_vals = []
        w2_errs = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            genome= Oligo.File.read_genome(orga,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose = 0)
            for chromo in genome:
                dist = dimer_dist(str(chromo),function_name,domain=domain,chromo=True)
                bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
                w1_vals.append(bias_vals[w1_ind])
                w2_vals.append(bias_vals[w2_ind])
                w1_errs.append(bias_errs[w1_ind])
                w2_errs.append(bias_errs[w2_ind])
        
        print(domain,zip(w1_vals, w1_errs, w2_vals, w2_errs))
        ax.scatter(x=w1_vals,y=w2_vals, label=function_titles[i], c=colors[i],s=2.5**2)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel(word2[0]+'p'+word2[1]+' bias')
    ax.set_xlabel(word1[0]+'p'+word1[1]+' bias')
    ax.axhline()
    ax.axvline()
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    #ax.set_xlim((-0.8,0.2))
    plt.title(file_label +" "+ word2[0]+'p'+word2[1]+' over '+word1[0]+'p'+word1[1]+' bias in '+domain)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + word1+"_"+word2+"_corr_"+ domain + '_v1'
    else:
        filename = + word1+"_"+word2+"_corr_"+ function_name+ '_v1'
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")

def plot_dimer_CpA_CpG_scatter(organism_names, domains, domain_list, group_names, function_name, function_title, colors=colors, file_label=""):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    CA_ind = dimer_dist.dimer_words.index("CA")
    CG_ind = dimer_dist.dimer_words.index("CG")

    for i,domain in tqdm(enumerate(domains)):
        CA_vals = []
        CA_errs = []
        CG_vals = []
        CG_errs = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            dist = dimer_dist(orga,function_name,domain=domain)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            CA_vals.append(bias_vals[CA_ind])
            CG_vals.append(bias_vals[CG_ind])
            CA_errs.append(bias_errs[CA_ind])
            CG_errs.append(bias_errs[CG_ind])
        
        print(domain,zip(CA_vals, CA_errs, CG_vals, CG_errs))
        ax.scatter(x=CA_vals, y=CG_vals, label=domain, c=colors[i])
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('CpG bias')
    ax.set_xlabel('CpA bias')
    ax.axhline()
    ax.axvline()
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    plt.title(file_label + " CpG over CpA bias in "+function_title)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "CpA_CpG_corr_"+ function_name + '_v1'
    else:
        filename = "CpA_CpG_corr_"+ function_name+ '_v1'
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")

def plot_dimer_CpG_vs_GC_scatter(organism_names, domains, domain_list, group_names, function_name, function_title, colors=colors, file_label="",TpA=False):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    TA_ind = dimer_dist.dimer_words.index("TA")
    CG_ind = dimer_dist.dimer_words.index("CG")

    for i,domain in tqdm(enumerate(domains)):
        TA_vals = []
        TA_errs = []
        CG_vals = []
        CG_errs = []
        GC_cont = []
        domain_orgas = [name for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        domain_idxs = [i for i,name in enumerate(organism_names) if (domain_list[i] == domain)]
        for j,orga in zip(domain_idxs,domain_orgas):
            print(domain)
            dist = dimer_dist(orga,function_name,domain=domain)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            TA_vals.append(bias_vals[TA_ind])
            CG_vals.append(bias_vals[CG_ind])
            TA_errs.append(bias_errs[TA_ind])
            CG_errs.append(bias_errs[CG_ind])
            row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,function_name,orga_name = orga, verbose=False)
            GC_cont.append( (row["G"] + row["C"]) / (row["G"] + row["C"] + row["A"] + row["T"]))
        
        print(CG_vals,"\n",GC_cont)
        if not TpA: ax.scatter(x=GC_cont, y=CG_vals, label=domain, c=colors[i])
        else:       ax.scatter(x=GC_cont, y=TA_vals, label=domain, c=colors[i])
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('CpG bias')
    if TpA: ax.set_ylabel('TpA bias')
    ax.set_xlabel('GC content')
    ax.axhline()
    ax.axvline(x=0.5)
    ax.grid(linestyle='-')
    #ax.set_ylim((-1,1))
    plt.title(file_label + " CpG bias over GC content in "+function_title)
    if TpA:     plt.title(file_label + " TpA bias over GC content in "+function_title)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    filename = "CpG_GCcont_corr_"+ function_name+ '_v1'
    if TpA: filename = "TpA" + filename [3:]
    if file_label is not None:
        filename = file_label + " " + filename
        
     
    plt.savefig(dimer_graph_folder+filename+".pdf")
    plt.savefig(dimer_graph_folder+filename+".png")