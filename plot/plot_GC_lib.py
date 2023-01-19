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
GC_bar_folder = "../../plots/elem_data/GC_cont/"

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

def sort_GC(organisms):
    sorted_orgas = []
    gc_arr = []
    for orga in organisms:
        domain=kmer_tools.download_genbank.name_to_domain(orga.name)
        row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,func_name="intergenics",orga_name = orga.name)  #to be replaced by "total_genome"
        length = row["length"]
        gc_arr.append ((row["G"]+row["C"])/length)

    sorting = np.argsort(gc_arr)
    #print(gc_arr)
    #gc_arr = [gc_arr[i] for i in sorting]
    #print("sorted: ",gc_arr)
    orgas_gc_sorted=[organisms[i] for i in sorting]
    return orgas_gc_sorted

def plot_GC_bar(organisms,func_name="genes",group_name="",sort_orga=True,sort_GC=False):
    gc = []
    at = []
    letters=['A','T','C','G']
    organism_names = [orga.name for orga in organisms]
    if(sort_orga):
        organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
        organisms = kmer_tools.download_genbank.names_to_orgas(organism_names)
    organism_labels = [kmer_tools.download_genbank.orga_name_normal_dict[name] for name in organism_names]
    print(organism_labels)

    for i,orga in enumerate(organisms):
        domain=kmer_tools.download_genbank.name_to_domain(orga.name)
        spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name=func_name,k_value=1,count_func=kmer_tools.search_kmer.word_num,orga_name = orga.name)
        values = []
        for letter in letters:
            values.append ( np.double(spec.get_count(letter)))
            #print(letter," - ",val)
            
        print(values)
        total = np.sum(values)
        print("total: ",total)
        at.append((values[0]+values[1])/total)
        gc.append((values[2]+values[3])/total)
    #for index, row in f.iterrows():
    #    chrom_name = row["chromosome"]
    #    #make bar plot: x-axis=words, y-axis = wordcount/expected
    #
    #    total = row["A"] + row ["T"] + row["G"] + row ["C"]
    #    gc.append((row["G"]+row["C"])/total)
    #    at.append((row["A"]+row["T"])/total)
    print(group_name+" - GC-AT : ",(gc),"-",(at))
    
    plt.figure(figsize=(6,8),dpi=300)
    ax=plt.subplot(1,1,1)

    col     = ["b","orange"]
    labels  = ["G+C","A+T"]
    width  = np.zeros(len(organism_names))
    left  = np.zeros(len(organism_names))
    for i,val in enumerate ([gc,at]):
        width = np.array(val)
        ax.barh(np.arange(len(organism_names)),width=width, height=0.8, left=left, color =col[i], label=labels[i])
        #width += 0.1
        left += width + 0.01
    
    ax.set_yticks(np.arange(len(organism_names)))
    ax.set_yticklabels(organism_names)
    ax.set_title(group_name+" - G+C content")
    #plt.xlim(0,1.4)
    ax.tick_params(axis='both', which='both', labelsize=6)
    plt.legend()
    plt.subplots_adjust(left=0.35,right=0.99)
    #plt.tight_layout()
    plt.savefig(GC_bar_folder+group_name+'_GC_content.png')
    plt.close()