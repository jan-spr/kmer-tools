from cmath import nan
import sys

from sqlalchemy import false

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

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data/"
heatmap_folder  = "../../plots/kmer_correlation/domain_heatmaps/"
heatmap_combined_folder  = "../../plots/kmer_correlation/domain_heatmaps/all_orgas/"
heatmap_repeats_folder  = "../../plots/kmer_correlation/domain_heatmaps/repeats/"
graph_folder    = "../../plots/kmer_correlation/k_correlation_orga/"
GC_bar_folder = "../../plots/elem_data/GC_cont/"
dimer_bar_folder = "../../plots/dimer/dimer_bar/"
mean_dev_folder = "../../plots/dimer/dimer_bar/barplot mean_bias/"
orga_bar_folder = "../../plots/dimer/dimer_bar/orga_bar/"

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

def plot_dimer_bars(bar_val_arr,bar_err_arr,labels,title,figsize=(10,5),label=None,width=1):
    labelsize="x-large"
    fig = plt.figure(figsize=figsize)
    ax=plt.subplot(1,1,1)
    n=len(labels)
    for i,label in enumerate(labels):
        vals = bar_val_arr[i]
        errs = bar_err_arr[i]
        if np.sum(np.abs(errs))==0: errs = None
        ax.bar(np.arange(16)-width/2+i*width/n,vals,yerr=errs,linewidth=0.5,alpha=1,label =label ,width=width/n,error_kw=dict(ecolor='gray', lw=1, capsize=2, capthick=1))
    ax.set_xticks(np.arange(16))
    ax.set_xticklabels(dimer_dist.dimer_words,fontsize=labelsize)
    ax.set_ylabel("dimer bias",fontsize=labelsize)
    ax.set_title(title,fontsize=labelsize)
    plt.tight_layout()
    plt.legend()
    return fig,ax




def plot_dimer_compare_model(orga_name,domain,func_name,folder=None,filename=None):
    plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1)
    distribution = dimer_dist(orga_name=orga_name,domain=domain,func_name=func_name)
    width = 1
    models=["mono_expectation", "count"]
    n=len(models)

    for i,model in enumerate(models):
        print(model)
        vals,err = distribution.get_dist(model)
        vals,err = distribution.normalize_count_arr(vals,err)
        print("vals = ",vals)
        print("err = ",err)
        if np.sum(np.abs(err))==0: err = None
        ax.bar(np.arange(16)-width/2+i*width/n,vals,yerr=err,linewidth=0.5,alpha=1,label =model ,width=width/n,error_kw=dict(ecolor='gray', lw=1, capsize=2, capthick=1))

    if folder is None:
        folder = dimer_bar_folder
    if filename is None:
        filename = "dimer_bar_compare "+orga_name+" "+func_name+".pdf"
    ax.set_xticks(np.arange(16))
    ax.set_xticklabels(distribution.dimer_words)
    ax.set_ylabel("normalized dimer count")
    ax.set_title(orga_name+" - "+"dimer counts and expectations in "+func_name)
    plt.tight_layout()
    plt.legend()
    plt.savefig(folder+filename)
    plt.close()

def plot_dimer_compare_elems(orga_name,domain,func_names,func_titles=None,folder=None,filename=None):
    if func_titles is None: func_titles = func_names
    plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1)
    distribution = dimer_dist(orga_name=orga_name,domain=domain,func_name=func_names[0])
    width = 1
    n=len(func_names)

    for i,func_name in enumerate(func_names):
        distribution.update_func(func_name)
        vals,err = distribution.get_dist("count")
        vals,err = distribution.normalize_count(vals,err)
        print("vals = ",vals)
        print("err = ",err)
        if np.sum(np.abs(err))==0: err = None
        ax.bar(np.arange(16)-width/2+i*width/n,vals,yerr=err,linewidth=0.5,alpha=1,label =func_titles[i] ,width=width/n,error_kw=dict(ecolor='gray', lw=1, capsize=2, capthick=1))

    if folder is None:
        folder = "../../plots/dimer/dimer_bar/"
    if filename is None:
        filename = "dimer_bar_compare "+kmer_tools.download_genbank.orga_name_normal(orga_name)+".png"
    ax.set_xticks(np.arange(16))
    ax.set_xticklabels(distribution.dimer_words)
    ax.set_ylabel("normalized dimer count")
    ax.set_title(kmer_tools.download_genbank.orga_name_normal(orga_name)+" - "+" dimer counts in different elements")
    plt.tight_layout()
    plt.legend()
    plt.savefig(folder+filename)
    plt.savefig(folder+filename+".pdf")
    plt.close()

def plot_dimer_mean_devs(orga_names,func_names,func_labels=None,domain=None,label=None):
    if domain is None: dom = True
    else: dom= False
    if func_labels is None: func_labels = func_names
    bar_vals_arr = []
    err_vals_arr = []
    for func_name in func_names:
        #print(orga_name)
        
        val_arrs = []
        err_arrs = []
        print(orga_names)
        for orga_name in orga_names:
            print(orga_name,domain)
            if dom: domain = kmer_tools.download_genbank.name_to_domain(orga_name)
            print(orga_name,domain)
            distribution = dimer_dist(orga_name,func_name,domain=domain,total=False)
            vals,errs = distribution.get_func_diff("count","mono_expectation")
            #vals,errs = distribution.normalize_count(distribution.get_dist("count"))
            #print("vals = ",vals)
            val_arrs.append(vals)
            err_arrs.append(errs)
        #print("val_arrs = ",val_arrs)
        mean_vals = np.mean(val_arrs,axis=0)
        std_vals  = np.std(val_arrs,axis=0)
        err_vals  = np.sqrt(np.sum(np.array(err_arrs)**2,axis=0))/len(orga_names)
        err_vals  = np.sqrt(std_vals**2 + err_vals**2)
        #print([(x,y) for x,y in zip(mean_vals,err_vals)])
        bar_vals_arr.append(mean_vals)
        err_vals_arr.append(err_vals)
    folder = mean_dev_folder
    title = "dimer mean deviation " + domain
    if label is not None: title = "dimer mean deviation " + label
    fig,ax = plot_dimer_bars(bar_vals_arr, err_vals_arr, labels=func_labels, title=title, figsize=(10,5), label=None, width=1)
    plt.savefig(folder+title+".png")
    plt.savefig(folder+title+".pdf")

def plot_dimer_sum_diffs(orga_names,func_name_pairs,domain=None,label=None):
    bar_vals_arr = []
    err_vals_arr = []
    for func_pair in func_name_pairs:
        #print(orga_name)
        
        val_arrs = []
        err_arrs = []
        for orga_name in tqdm(orga_names):
            distribution1 = dimer_dist(orga_name,func_pair[0],domain=domain,total=False)
            distribution2 = dimer_dist(orga_name,func_pair[1],domain=domain,total=False)
            vals1,errs1 = distribution1.get_func_diff("count","mono_expectation")
            vals2,errs2 = distribution2.get_func_diff("count","mono_expectation")
            vals_f = vals1 - vals2
            errs_f = np.sqrt(errs1**2 + errs2**2)
            #vals,errs = distribution.normalize_count(distribution.get_dist("count"))
            #print("vals = ",vals)
            val_arrs.append(vals_f)
            err_arrs.append(errs_f)
        #print("val_arrs = ",val_arrs)
        mean_vals = np.mean(val_arrs,axis=0)
        std_vals  = np.std(val_arrs,axis=0)
        err_vals  = np.sqrt(np.sum(np.array(err_arrs)**2,axis=0))/len(orga_names)
        err_vals  = np.sqrt(std_vals**2 + err_vals**2)
        #print([(x,y) for x,y in zip(mean_vals,err_vals)])
        bar_vals_arr.append(mean_vals)
        err_vals_arr.append(err_vals)
    folder = mean_dev_folder
    title = "dimer_region_diff_" + domain+".png"
    fig,ax = plot_dimer_bars(bar_vals_arr, err_vals_arr, labels=["-".join(func_pair) for func_pair in func_name_pairs], title=title, figsize=(10,5), label="coding regons comparison", width=1)
    plt.savefig(folder+title)


def plot_dimer_orga_categ_bias(orga_name,func_names,func_name_titles,domain=None,label=None):
    if domain is None: 
        domain = kmer_tools.download_genbank.name_to_domain(orga_name)
    bar_vals_arr = []
    err_vals_arr = []
    val_arrs = []
    err_arrs = []
    for func_name in func_names:
        #print(orga_name)
        
        distribution = dimer_dist(orga_name,func_name,domain=domain,total=False)
        vals,errs = distribution.get_func_diff("count","mono_expectation")
        #vals,errs = distribution.normalize_count(distribution.get_dist("count"))
        #print("vals = ",vals)
        val_arrs.append(vals)
        err_arrs.append(errs)
        #print("val_arrs = ",val_arrs)
    folder = dimer_bar_folder
    filename = domain+ " dimer_orga_categ_bias " + orga_name +".png"
    title = orga_name + " dimer bias comparison"
    fig,ax = plot_dimer_bars(val_arrs, err_arrs, labels=func_name_titles, title=title, figsize=(10,5), label=None, width=1)
    plt.savefig(folder+filename)

def plot_dimer_categ_bias_box(organism_names,func_names,func_name_titles,domain=None,label=None):
    if domain is None: 
        domain = kmer_tools.download_genbank.name_to_domain(organism_names[0])
    bar_vals_arrs = []
    err_vals_arrs = []
    for func_name in func_names:
        #print(orga_name)
        val_arrs = []
        err_arrs = []
        for orga_name in tqdm(organism_names):
            distribution = dimer_dist(orga_name,func_name,domain=domain,total=False)
            vals,errs = distribution.get_func_diff("count","mono_expectation")
            #vals,errs = distribution.normalize_count(distribution.get_dist("count"))
            #print("vals = ",vals)
            val_arrs.append(vals)
            err_arrs.append(errs)
        #print("val_arrs = ",val_arrs)
        #print([(x,y) for x,y in zip(mean_vals,err_vals)])
        bar_vals_arrs.append(val_arrs)
        err_vals_arrs.append(err_arrs)
    folder = dimer_bar_folder
    filename = domain+ " dimer_categ_bias_box.png"
    title = domain + " dimer bias comparison"
    kmer_tools.plot.plot_double_boxplot(bar_vals_arrs[0],bar_vals_arrs[1],labels=func_name_titles,xlabel_array=dimer_dist.dimer_words,y_label="dimer bias (norm.)",title=title,filepath=folder+filename)


def plot_dimer_word_bias_orgas(orga_names,func_name,func_name_title,Word,label=None,
    domain_list=None,sorted=True,filter_empty=False,filer_Nan=False,split_class=False,domain=None,mark_wrap=False):
    fontsize='x-large'
    dimer_words=dimer_dist.dimer_words
    index = dimer_words.index(Word)
    vals_arr = []
    err_arr = []
    for orga_name in tqdm(orga_names):
        distribution1 = dimer_dist(orga_name,func_name,domain=domain,total=False)
        dimer_vals,dimer_errs = distribution1.get_func_diff("count","mono_expectation")

        #vals,errs = distribution.normalize_count(distribution.get_dist("count"))
        #print("vals = ",vals)
        vals_arr.append(dimer_vals[index])
        err_arr.append(dimer_errs[index])
    #print("val_arrs = ",val_arrs)


    folder = orga_bar_folder
    title = "dimer Orga Bias" + Word+" "+func_name+".png"
    fig = plt.figure(figsize=(8,4.5))
    ax=plt.subplot(1,1,1)
    n=len(orga_names)
    width=0.85
    
    ax.bar(np.arange(n),vals_arr,yerr=err_arr,linewidth=0.5,alpha=1,label =label ,width=width,error_kw=dict(ecolor='gray', lw=1, capsize=2, capthick=1))
    ax.set_xticks(np.arange(n))
    #ax.set_xticklabels([name.split(" ")[0]+" "+name.split(" ")[1][0]+"." for name in orga_names],rotation=90)
    ax.set_ylabel("dimer bias ", fontsize=fontsize)
    ax.set_title("dimer bias of "+Word+" in "+func_name_title+" over organisms", fontsize=fontsize)
    plt.tight_layout()
    ylim = ax.get_ylim()
    #plt.legend()

    #kmer_tools.plot.add_domain_split(orga_names,orga_names,ax,split_x=True,split_y=False,labels=False)
    kmer_tools.plot.add_sublabel_split(orga_names,orga_names,fig,ax,split_x=True,split_y=False,show_labels=True,class_split=False,mark_wrap=False)
    ax.tick_params(axis='both', which='major', labelsize='large')
    ax.tick_params(axis='both', which='minor', labelsize='large')

    ax.set_ylim(ylim)
    ax.set_xlim((0,len(orga_names)))

    plt.savefig(folder+title)
    #return fig,ax
    


function_names=[
    'genes',
    'exons',
    'introns',
    'cds',
    'mRNAs',
    'intergenics',
    'protein_coding_genes',
    'pseudo_genes',
    #'pseudo_genes_exons',
    #'pseudo_genes_introns',
    'non_protein_coding_genes',
    'ncRNAs',
    'lncRNAs',
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'tRNAs',
    'rRNAs',
    'precursor_RNAs',
    #'gaps',
]
repeat_funcs_labels = kmer_tools.search_kmer.repeat_funcs_labels
for func_name in function_names+repeat_funcs_labels:
    #plot_dimer_compare_model(orga_name="Rattus norvegicus",domain = "animalia",func_name=func_name)
    pass