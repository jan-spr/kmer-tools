from cmath import nan
import sys
from kmer_tools.search_kmer.reorganize_funcs import k_words_expand, kmer_list

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
GC_bar_folder = "../../plots/elem_data/GC_cont/"
dimer_bar_folder = "../../plots/elem_data/dimer_bar/"

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


def plot_word_bar_deviation_mean_categ(organisms,function_names,function_titles,domain,count_func,colors,words=None,error = False,folder=None):
    plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1)
    n=len(function_names)
    if(words is None): 
      #words = ["AA",	"TT",	"AT",	"TA",	"AC",	"CA",	"TG",	"GT",	"AG",	"GA",	"TC",	"CT",	"CC",	"GG",	"GC",	"CG"]
      words = kmer_list(k_value=2)
    width=1
    for i,func_name in tqdm(enumerate(function_names)):
        count_all = []
        err_all = []
        for j,orga in enumerate(organisms):
            domain=kmer_tools.download_genbank.name_to_domain(orga.name)
            counts,err = kmer_tools.search_kmer.read_kdat_orgaFunc(k_data_folder,domain,func_name=func_name,k_value=2,count_func=count_func,orga_name = orga.name)

            #make bar plot: x-axis=words, y-axis = wordcount/expected

            count_all.append(counts)
            err_all.append(err)

        d_count_all = np.std(np.array(count_all),axis=0)
        count_all   = np.mean(np.array(count_all),axis=0)

        #print(np.arange(len(words))+i*width/n)
        ax.bar(np.arange(len(words))-width/2+i*width/n,count_all,yerr=d_count_all,linewidth=0.5,alpha=1,label =function_titles[i] ,width=width/n,error_kw=dict(ecolor='gray', lw=1, capsize=2, capthick=1))
    ax.set_xticks(np.arange(len(words)))
    ax.set_xticklabels(words)
    ax.set_ylabel("data - expected count (mean)")
    ax.set_title(domain+" - "+"mean deviation from expectation in each element")
    plt.tight_layout()
    plt.legend()
    if (folder == None): folder = "deviation_mean_combined"
    plt.savefig((dimer_bar_folder+ folder  +"/"+domain+'_'+func_name+'_absolute.png'))
    plt.close()