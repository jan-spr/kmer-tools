from cmath import nan
import sys

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

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data/"
heatmap_folder  = "../../plots/kmer_correlation/domain_heatmaps/"
heatmap_combined_folder  = "../../plots/kmer_correlation/domain_heatmaps/all_orgas/"
heatmap_repeats_folder  = "../../plots/kmer_correlation/domain_heatmaps/repeats/"
graph_folder    = "../../plots/kmer_correlation/k_correlation_orga/"
coverage_folder = "../../plots/elem_data/coverage/"
coverage_hm_folder = "../../plots/elem_data/coverage_hm/"
elem_plot_folder = "../../plots/elem_data/"
box_compare_folder = "../../plots/elem_data/box_compare/"

def plot_boxplot(value_array,label_array,x_label,y_label,title,filepath):
    fig, ax = plt.subplots(figsize=(12,4))

    bp = ax.boxplot(value_array,positions=np.arange(len(label_array)),labels = label_array)
    ax.set_ylabel(y_label)
    ax.axhline(y=0,color="black",linestyle="dotted")  
    ax.set_title(title,fontsize=20)
    #plt.tight_layout()
    plt.savefig(filepath)
    plt.close()

def plot_double_boxplot(value_array1,value_array2,labels,xlabel_array,x_label=None,y_label=None,title=None,filepath="",colors=['blue','red']):
    fig, ax = plt.subplots(figsize=(12,4))
    bp=[]
    offset=0/(4*3)
    medianprops =  dict(linestyle='--', linewidth=0, color='purple')


    meta_xlabel_array = [xlabel_array,None]
    for k,value_array in enumerate([value_array1,value_array2]):
        
        value_array  = np.transpose(value_array)
        print(value_array,meta_xlabel_array[k])
        print(len(value_array),[len(i) for i in value_array])
        print(colors[k])
    value_array1  = np.transpose(value_array1)
    value_array2  = np.transpose(value_array2)
    bp.append(ax.boxplot(value_array1.tolist(),positions=np.arange(len(value_array1))+offset,labels=xlabel_array,medianprops=medianprops,boxprops=dict(color=colors[0]),showfliers=False))
    bp.append(ax.boxplot(value_array2.tolist(),positions=np.arange(len(value_array2))-offset,labels=xlabel_array,medianprops=medianprops,boxprops=dict(color=colors[1]),showfliers=False,manage_ticks=False))
        
        
    ax.legend([bp[0]["boxes"][0], bp[1]["boxes"][0]], labels, loc='upper right')
    ax.set_ylabel(y_label)
    ax.axhline(y=0,color="black",linestyle="dotted") 
    ax.set_title(title,fontsize=20)
    #plt.xlim(0,1.4)
    #plt.legend()
    #plt.tight_layout()
    plt.savefig(filepath)
    plt.close()

def plot_n_boxplot(value_arrays,labels,xlabel_array,x_label,y_label,title,filepath,colors=['blue','red']):
    fig, ax = plt.subplots(figsize=(10,5))
    bp=[]
    offset=0/(4*3)
    medianprops =  dict(linestyle='--', linewidth=0, color='purple')


    meta_xlabel_array = [xlabel_array,None]
    for k,value_array in enumerate(value_arrays):
        
        value_array  = np.transpose(value_array)
        print(value_array,meta_xlabel_array[k])
        print(len(value_array),[len(i) for i in value_array])
        print(colors[k])
        value_array  = np.transpose(value_array)
    
        
        bp.append(ax.boxplot(value_array.tolist(),positions=np.arange(len(value_array))-offset,labels=xlabel_array,medianprops=medianprops,boxprops=dict(color=colors[k]),showfliers=False,manage_ticks=k==0))
        
    ax.legend([bp[0]["boxes"][0], bp[1]["boxes"][0]], labels, loc='upper right')
    ax.legend([bp[i]["boxes"][0] for i in range(len(value_arrays))], labels, loc='upper right')
    ax.set_ylabel(y_label)
    ax.axhline(y=0,color="black",linestyle="dotted") 
    ax.set_title(title,fontsize=20)
    #plt.xlim(0,1.4)
    #plt.legend()
    #plt.tight_layout()
    plt.savefig(filepath)
    plt.close()

def plot_dimerBox_compare(organisms,domain,func_name1,func_name2):
    
    dimer_words = ["AA",	"TT",	"AT",	"TA",	"AC",	"CA",	"TG",	"GT",	"AG",	"GA",	"TC",	"CT",	"CC",	"GG",	"GC",	"CG"]
    func_values1 = []
    func_values2 = []

    for orga in tqdm(organisms):
        for value_arr,name in zip([func_values1, func_values2],
                                    [func_name1,func_name2]):
            row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,func_name=name,orga_name = orga.name)
            counts_norm,_ = kmer_tools.search_kmer.dimer_count_normalized(row,dimer_words)
            value_arr.append(counts_norm)
    print(len(func_values1), len(func_values2))
    print(len(func_values1[0]), len(func_values2[0]))

    ylabel = "normalized dimer count"
    xlabel = None
    title = func_name1 + "-" + func_name2 + "  comparison - " + domain
    path = box_compare_folder
    filename = domain + "_" + func_name1 + "-" + func_name2 + "dimer_boxplot.png"
    plot_double_boxplot(func_values1, func_values2,[func_name1,func_name2],dimer_words,xlabel,ylabel,title,filepath= path + filename,colors=['blue','red'])

