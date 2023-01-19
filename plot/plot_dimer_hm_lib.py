from cmath import nan
import sys

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import seaborn as sns
from tqdm import tqdm
from kmer_tools.plot.dimer_func_lib import *
from kmer_tools.plot.plot_hm_lib import plot_heatmap,add_wrap_split,fit_cbar_height

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data/"
dimer_bar_folder = "../../plots/dimer/dimer_bar/"
mean_dev_folder = "../../plots/dimer/dimer_bar/mean_dev/"
dimer_hm_folder = "../../plots/dimer/dimer_hm/"

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

def nan_filter(matrix):
    nan_arr = [np.isnan(row).any() for row in matrix]
    print(nan_arr)
    nan_ind = [i for i,bool in zip(np.arange(len(nan_arr)), nan_arr) if not bool]
    matrix_filtered = [spec for spec, bool in zip(matrix, nan_arr) if not bool]
    return nan_ind, matrix_filtered

def plot_dimer_word_count(orga_names,func_name,func_label,domain=None,domain_list=None,label=None,filter_nan=False,dimer_bias=False,fix_axis=True):
    matrix=[]
    for orga_name in tqdm(orga_names):
        row=[]
        dist = dimer_dist(orga_name,func_name,domain=domain)
        count_vals,_ = dist.normalize_count(dist.get_dist("count"))
        if dimer_bias:
            count_vals,_ = dist.get_func_diff("count","mono_expectation")
        matrix.append(count_vals)
    nan_ind = np.arange(len(orga_names))
    matrix = np.array(matrix)
    matrix [matrix == 0] = np.nan 
    if filter_nan:
        nan_ind, matrix = nan_filter(matrix)
    print(matrix)
    labels_orga = [kmer_tools.download_genbank.orga_name_normal(name) for name in orga_names]
    labels_orga = [labels_orga[i] for i in nan_ind]
    #labels_orga =[name.split(" ")[0]+" "+name.split(" ")[1][0]+"." for name in labels_orga]
    labels_words = dimer_dist.dimer_words
    folder = dimer_hm_folder
    filename = func_label + " " + "dimer_count_hm_"
    if label is not None: filename += " "+ label + " "
    filename += domain + ".png"

    v_min = v_max = None
    if fix_axis:
        if dimer_bias:
            v_min = -0.6
            v_max = 0.6
        else:
            v_min =0
            v_max = 2.5

    filter_ind=np.arange(len(matrix))
    if True:
        filter_ind,matrix = nan_filter(matrix)
        print(filter_ind,len(filter_ind))
        if len(filter_ind)==0: return 0
    labels_orga = [labels_orga[i] for i in filter_ind]
    fig_height= 9/62*len(labels_orga)
    if len(labels_orga) < 40:
        fig_height= 9/62*len(labels_orga)*1.5
    if fig_height<3:fig_height=3

    fig,ax = plot_heatmap(matrix,figsize=(6,fig_height),labels_y=labels_orga,labels_x=labels_words,output_file=folder+filename,square=False,split_y=True,split_class=True,cmap='rainbow',v_min=v_min,v_max=v_max)
    #add_wrap_split(fig,ax,orga_names,labels_words,func_names1,func_names2,split_x=True,split_y=True,show_labels=False,class_split=False,mark_wrap=False,add_title=True)
    
    plt.savefig(folder+filename)

def plot_dimer_cum_bias(orga_names,func_names,func_labels,domain=None,domain_list=None,label=None,filter_nan=False,total_exp= False):
    matrix=[]
    for orga_name in tqdm(orga_names):
        row=[]
        for func_name in func_names:
            dist = dimer_dist(orga_name,func_name,domain=domain)
            bias_vals,bias_errs = dist.get_func_diff("count","mono_expectation")
            if total_exp : bias_vals,bias_errs = dist.get_func_diff("count","total_expectation")
            cum_val = np.sum(np.abs(bias_vals))
            cum_err =  np.sqrt(np.sum(bias_errs**2))
            row.append(cum_val)
        matrix.append(row)
    nan_ind = np.arange(len(orga_names))
    matrix = np.array(matrix)
    matrix [matrix == 0] = np.nan 
    if filter_nan:
        nan_ind, matrix = nan_filter(matrix)
    print(matrix)
    labels_orga = [kmer_tools.download_genbank.orga_name_normal_dict[name] for name in orga_names]
    labels_orga = [labels_orga[i] for i in nan_ind]
    labels_func = func_labels
    folder = dimer_hm_folder
    filename = "dimer_cum_bias_hm_"
    if label is not None: filename += " "+ label + " "
    if total_exp: filename += " total_exp "
    filename += domain + ".png"
    fig,ax = plot_heatmap(matrix,figsize=(5.5,10),labels_y=labels_orga,labels_x=labels_func,output_file=folder+filename,square=False,split_y=True,split_class=True,cmap='rainbow',v_max=8)
    plt.savefig(folder+filename)


def plot_dimer_cum_bias_diff(orga_names,func_pair_names,func_pair_labels,domain=None,domain_list=None,label=None,filter_nan=True):
    matrix=[]
    for orga_name in tqdm(orga_names):
        row=[]
        for name_pair in func_pair_names:
            val = dimer_bias_dist_cum(orga_name,name_pair[0],name_pair[1],domain=domain)[0]
            row.append(val)
        matrix.append(row)
    if filter_nan:
        nan_ind, matrix = nan_filter(matrix)
    print(matrix)
    labels_orga = [kmer_tools.download_genbank.orga_name_normal_dict[name] for name in orga_names]
    labels_orga = [labels_orga[i] for i in nan_ind]
    labels_func = [" - ".join(names) for names in func_pair_labels]
    folder = dimer_hm_folder
    filename = "dimer_dev_bias_hm_"
    if label is not None: filename += " "+ label + " "
    filename += domain + ".png"
    fig,ax = plot_heatmap(matrix,figsize=(5,8),labels_y=labels_orga,labels_x=labels_func,output_file=folder+filename,square=False,split_y=True,split_class=True,cmap='rainbow')
    plt.savefig(folder+filename)

def plot_dimer_covariance_heatmap(orga_names,func_names,func_labels=None,domain=None,domain_list=None,organisms_label=None,label=None,filter_nan=False,verbose=0):
    if func_labels is  None: func_labels = func_names
    for i,func_name in enumerate(func_names):
        if verbose: print(func_labels[i])
        dimer_words= dimer_dist.dimer_words
        dimer_bias_vals=[]
        for j,orga_name in enumerate(orga_names):
            if domain_list is not None: domain = domain_list[j]
            dist=dimer_dist(orga_name,func_name,domain=domain)
            count_vals,_ = dist.get_func_diff("count","mono_expectation")
            dimer_bias_vals.append(count_vals)
        dimer_bias_vals = np.transpose(dimer_bias_vals)
        matrix = np.corrcoef(dimer_bias_vals)
        if verbose: print(domain)
        if verbose: print(DataFrame(matrix,index=dimer_words,columns=dimer_words))
        if organisms_label is None: organisms_label = domain
        fig,ax = plot_heatmap(matrix,figsize=(6,6),labels_y=dimer_words,labels_x=dimer_words,square=True,split_y=False,split_class=False,cmap='RdBu',v_min=-1,v_max=1)
        ax[0].set_title("correlation of dimer biases in "+func_name+ " of " +organisms_label,fontsize="large")
        
        folder = dimer_hm_folder
        filename = "dimer_covariance_hm_"
        if label is not None: filename += " "+ label + " "
        filename += organisms_label + " " + func_name + ".png"
        plt.tight_layout()
        fit_cbar_height(ax[0],ax[1])
        ax[1].set_frame_on(True)
        plt.savefig(folder+filename)
        
def plot_dimer_covariance_heatmap_combined(orga_names,func_names,func_labels=None,domain=None,domain_list=None,organisms_label=None,label=None,filter_nan=False,verbose=0):
    if func_labels is  None: func_labels = func_names
    for i,func_name in enumerate(func_names):
        if verbose: print(func_labels[i])
        dimer_words= dimer_dist.dimer_words
        dimer_bias_vals=[]
        for j,orga_name in enumerate(orga_names):
            if domain_list is not None: domain = domain_list[j]
            dist=dimer_dist(orga_name,func_name,domain=domain)
            count_vals,_ = dist.get_func_diff("count","mono_expectation")
            dimer_bias_vals.append(count_vals)
    dimer_bias_vals = np.transpose(dimer_bias_vals)
    matrix = np.corrcoef(dimer_bias_vals)
    if verbose: print(domain)
    if verbose: print(DataFrame(matrix,index=dimer_words,columns=dimer_words))
    if organisms_label is None: organisms_label = domain
    fig,ax = plot_heatmap(matrix,figsize=(6,6),labels_y=dimer_words,labels_x=dimer_words,square=True,split_y=False,split_class=False,cmap='RdBu',v_min=-1,v_max=1)
    ax[0].set_title("correlation of dimer biases over "+ label + " of " +organisms_label)
    
    folder = dimer_hm_folder
    filename = "dimer_covariance_hm_"
    if label is not None: filename += " "+ label + " "
    filename += organisms_label + " " + func_name + ".png"
    plt.tight_layout()
    fit_cbar_height(ax[0],ax[1])
    ax[1].set_frame_on(True)
    plt.savefig(folder+filename)

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