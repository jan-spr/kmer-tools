import sys
from tabnanny import verbose

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

from cmath import nan
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from plot_hm_lib import get_orga_kmer_heatmap

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data/"
#k_data_folder   = "../../data/k_data_strand/"
graph_folder    = "../../plots/kmer_correlation/correlation_graph/k_correlation_orga/"
compare_graph_folder    = "../../plots/kmer_correlation/correlation_graph/k_correlation_func/"
orga_group_folder    = "../../plots/kmer_correlation/correlation_graph/k_correlation_func/orgaGroups/"
elem_plot_folder = "../../plots/elem_data/"


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

def plot_corr_k(domain,organism_names,k_values,function_names,colors):
    #drawers =[]
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    for i,func_name in tqdm(enumerate(function_names)):
        print(func_name)
        mean_values = []
        err_values = []
        for k in k_values:
            matrix = get_orga_kmer_heatmap(organism_names,func_name,function_names[0],k,domain = domain)
            
            #print(matrix.matrix)

            print(np.isnan(matrix.matrix).sum())
            print(matrix.matrix.size)
            if (np.isnan(matrix.matrix).sum()/matrix.matrix.size < 0.5):
                m = np.nanmean(matrix.matrix)
                e = np.nanstd(matrix.matrix)/np.sqrt(len(matrix))
            else:
                m = e = nan
            print("mean, error: ",m,e)
            mean_values.append(m)
            err_values.append(e)
        ax.errorbar(x=np.array(k_values), y=np.array(mean_values), yerr=np.array(err_values), label=func_name, 
            fmt='o--',  linestyle=(linestyles+linestyles)[i],alpha=0.6,capsize=5)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('mean correlation')
    ax.set_xlabel('word size k [bp]')
    ax.set_xticks(k_values)
    #ax.set_ylim((-1,1))
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    plt.savefig(graph_folder+domain+'_v1.png')



def plot_corr_k_categ(domain,organisms,organism_names,k_values,function_pairs,function_pair_titles,colors,file_label=None):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    for i,func_pair in tqdm(enumerate(function_pairs)):
        mean_values = []
        err_values = []
        for k in k_values:
            print("k = ",k,":")
            #orga_names = []
            corr_vals = []
            for orga in organisms:
                genome = Oligo.File.read_genome(orga.name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
                
                spectra_funcs =[]
                for j,func_name in enumerate(func_pair):
                    spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k,orga_name=orga.name)
                    
                    if (True and (spec.k is None)):
                        kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(genome[0]),"genes",k)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
                        spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k)
                    spectra_funcs.append( spec)
                    
                try:
                    #print(Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix,Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix.item(1))
                    corr_val = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs,verbose=0).matrix.item(1)
                    #if k==1 or k== 2:   
                        #print(orga.name,func_pair)
                        #print([spec.get_bins() for spec in spectra_funcs])
                    if ((kmer_tools.search_kmer.get_spec_coverage(spectra_funcs [0])<5000 or kmer_tools.search_kmer.get_spec_coverage(spectra_funcs [1])<5000)):
                        continue
                    print("corr value: ",corr_val)
                    corr_vals.append(corr_val)
                except:
                    continue
            
        
            print("correlations:",func_pair," - ",corr_vals)
            corr_vals = [val for val in corr_vals if val == val]
            m = np.nanmean(corr_vals)
            e = np.nanstd(corr_vals)/np.sqrt(len(corr_vals))
            mean_values.append(m)
            err_values.append(e)
        print(func_pair,mean_values)
        print(len(k_values),len(mean_values),len(err_values))
        print(k_values,mean_values,err_values)
        ax.errorbar(x=np.array(k_values), y=np.array(mean_values), yerr=np.array(err_values), label=function_pair_titles[i][0]+" - "+function_pair_titles[i][1], 
        fmt='o--',  linestyle=linestyles[i],alpha=0.6,capsize=5)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('mean correlation')
    ax.set_xlabel('word size k [bp]')
    ax.set_xticks(k_values)
    #ax.set_ylim((-1,1))
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "_" + domain + '_v1.pdf'
    else:
        filename = domain + '_v1.pdf'
    plt.savefig(compare_graph_folder+filename)
    plt.savefig(compare_graph_folder+filename+".png")


def plot_corr_OrgaGroup(organism_names, domain_list, group_names, orga_to_group_list, k_value, function_pairs, function_pair_titles, colors, file_label=None):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    for i,func_pair in tqdm(enumerate(function_pairs)):
        mean_values = []
        err_values = []
        for group_name in group_names:
            print("group_name = ",group_name,":")
            #orga_names = []
            corr_vals = []
            group_orgas = [name for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            group_idxs = [i for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            for j,orga in zip(group_idxs,group_orgas):
                domain = domain_list[j]
                print(domain)
                genome = Oligo.File.read_genome(orga,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
                func_spec_pair =[]
                for k,func_name in enumerate(func_pair):
                    spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga)
                    if (spec.k is None):
                        spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,"genes",k_value,orga_name=orga)
                        kmers = spec.get_kmers()
                        spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
                    func_spec_pair.append(spec)
                    
                    try:
                        #print(Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix,Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix.item(1))
                        corr_val = Oligo.Kmer.KmerSpectrum.correlate_spectra(func_spec_pair,verbose=0).matrix.item(1)
                        #if k==1 or k== 2:   
                            #print(orga.name,func_pair)
                            #print([spec.get_bins() for spec in spectra_funcs])
                        #print("corr value: ",corr_val)
                        corr_vals.append(corr_val)
                    except:
                        continue
            
        
            print("correlations:",func_pair," - ",corr_vals)
            corr_vals = [val for val in corr_vals if val == val]
            m = np.nanmean(corr_vals)
            e = np.nanstd(corr_vals)
            mean_values.append(m)
            err_values.append(e)
        print(func_pair,mean_values)
        print(len(group_names),len(mean_values),len(err_values))
        print(group_names,mean_values,err_values)
        ax.errorbar(x=np.arange(len(group_names)), y=np.array(mean_values), yerr=np.array(err_values), label=function_pair_titles[i][0]+" - "+function_pair_titles[i][1], 
        fmt='o--',  linestyle=linestyles[i],alpha=0.6,capsize=5)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('mean correlation')
    ax.set_xticks(np.arange(len(group_names)))
    ax.set_xticklabels(group_names)
    #ax.set_ylim((-1,1))
    plt.title(file_label + " k="+str(k_value))
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "_k="+str(k_value) + '_v1.pdf'
    else:
        filename = domain + "_k="+str(k_value) + '_v1.pdf'
    plt.savefig(orga_group_folder+filename)
    plt.savefig(orga_group_folder+filename+'.png')

def plot_corr_OrgaGroup_k(organism_names, domain_list, group_names, orga_to_group_list, k_values, function_pairs, function_pair_titles, colors, file_label=None):
    fig, ax = plt.subplots(figsize=(9,5.5),dpi=600)
    for k_value in k_values:
        for i,func_pair in tqdm(enumerate(function_pairs)):
            mean_values = []
            err_values = []
            for group_name in group_names:
                print("group_name = ",group_name,":")
                #orga_names = []
                corr_vals = []
                group_orgas = [name for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
                group_idxs = [i for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
                for j,orga in zip(group_idxs,group_orgas):
                    domain = domain_list[j]
                    genome = Oligo.File.read_genome(orga,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
                    func_spec_pair =[]
                    for k,func_name in enumerate(func_pair):
                        spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga)
                        
                        #row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,func_name,orga_name=orga)
                        #elem_len = row["A"] + row["T"]+ row["C"]+ row["G"]
                        mark_entry =False
                        
                        if (spec.k is None or kmer_tools.search_kmer.get_spec_coverage(spec) < 5000):
                            #spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,"genes",k_value,orga_name=orga)
                            #kmers = spec.get_kmers()
                            #spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
                            mark_entry = True
                        func_spec_pair.append(spec)
  
                    try:
                        #print(Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix,Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix.item(1))
                        if mark_entry:
                            corr_val = np.nan
                        else:
                            corr_val = Oligo.Kmer.KmerSpectrum.correlate_spectra(func_spec_pair,verbose=0).matrix.item(1)
                        #if k==1 or k== 2:   
                            #print(orga.name,func_pair)
                            #print([spec.get_bins() for spec in spectra_funcs])
                        #print("corr value: ",corr_val)
                        corr_vals.append(corr_val)
                    except:
                        continue
                
            
                print("correlations:",func_pair," - ",corr_vals)
                corr_vals = [val for val in corr_vals if val == val]
                m = np.nanmean(corr_vals)
                e = np.nanstd(corr_vals)
                mean_values.append(m)
                err_values.append(e)
            print(func_pair,mean_values)
            print(len(group_names),len(mean_values),len(err_values))
            print(group_names,mean_values,err_values)
            ax.errorbar(x=np.arange(len(group_names)), y=np.array(mean_values), yerr=np.array(err_values), label=function_pair_titles[i][0]+" - "+function_pair_titles[i][1]+" k="+str(k_value), 
            fmt='o--',  linestyle=linestyles[i],alpha=0.6,capsize=5)
    #drawer = Oligo.Plot.MultiDrawer(drawers)
    ax.set_ylabel('mean correlation')
    ax.set_xticks(np.arange(len(group_names)))
    ax.set_xticklabels(group_names)
    #ax.set_ylim((0.1,1.1))
    plt.title(file_label)
    plt.legend()
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + "_k="+''.join(str(k) for k in k_values) + '_v1.pdf'
    else:
        filename = domain + "_k="+''.join(str(k) for k in k_values) + '_v1.pdf'
    plt.savefig(orga_group_folder+filename)
    plt.savefig(orga_group_folder+filename+".png")


def plot_corr_OrgaGroup_box(organism_names, domain_list, group_names, orga_to_group_list, k_value, function_pairs, function_pair_titles, colors, file_label=None):
    box_values_all = []
    for i,func_pair in tqdm(enumerate(function_pairs)):
        value_array = []
        mean_values = []
        err_values = []
        for group_name in group_names:
            print("group_name = ",group_name,":")
            #orga_names = []
            corr_vals = []
            group_orgas = [name for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            group_idxs = [i for i,name in enumerate(organism_names) if (orga_to_group_list[i] == group_name)]
            for j,orga in zip(group_idxs,group_orgas):
                domain = domain_list[j]
                genome = Oligo.File.read_genome(orga,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
                func_spec_pair =[]
                for k,func_name in enumerate(func_pair):
                    spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga)
                    if (spec.k is None):
                        spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,"genes",k_value,orga_name=orga)
                        kmers = spec.get_kmers()
                        spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
                    func_spec_pair.append(spec)
                    
                    try:
                        #print(Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix,Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs).matrix.item(1))
                        corr_val = Oligo.Kmer.KmerSpectrum.correlate_spectra(func_spec_pair,verbose=0).matrix.item(1)
                        #if k==1 or k== 2:   
                            #print(orga.name,func_pair)
                            #print([spec.get_bins() for spec in spectra_funcs])
                        #print("corr value: ",corr_val)
                        corr_vals.append(corr_val)
                    except:
                        continue
            
        
            print("correlations:",func_pair," - ",corr_vals)
            corr_vals = [val for val in corr_vals if val == val]
            value_array.append(corr_vals)
            m = np.nanmean(corr_vals)
            e = np.nanstd(corr_vals)
            mean_values.append(m)
            err_values.append(e)
        print(func_pair,mean_values)
        print(len(group_names),len(mean_values),len(err_values))
        print(group_names,mean_values,err_values)
        box_values_all.append(value_array)
        
    #drawer.plot('./plots/correlation/chromosomal/correlation_vs_k/'+domain+'_.png', 
    #  xticks=k_values, xlabel='word size k [bp]', ylabel='mean correlation', ylim=(-1,1),figsize=(10,10))
    if (file_label is not None):
        filename = file_label + '_bv1.pdf'
    else:
        filename = domain + '_bv1.pdf'
    labels = [function_pair_titles[i][0]+" - "+function_pair_titles[i][1] for i in range(len(function_pair_titles))]
    kmer_tools.plot.plot_box_lib.plot_n_boxplot(value_arrays=box_values_all, labels=labels,xlabel_array=group_names,x_label=None,y_label='mean correlation',title=None,filepath=orga_group_folder+filename)

def plot_bar_length(organism_names, domain_list,title):
    fig, ax = plt.subplots(figsize=(10,5),dpi=600)
    length_arr = []
    for j,orga_name in tqdm(enumerate(organism_names)):
        length = kmer_tools.search_elem_stat.read_prop_orga(orga_name,"length",domain_list[j],"gaps")
        length_arr.append(length)
    
    ax.bar(np.arange(len(organism_names)),length_arr)
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(organism_names)))
    labels = [kmer_tools.download_genbank.orga_name_normal_dict[name] for name in organism_names]
    labels = [label.split(" ")[0] + " " + label.split(" ")[1][0] + "." for label in labels]
    ax.set_xticklabels(labels,rotation="vertical")
    ax.set_title(title)
    ax.set_ylabel("length (log) [bp]")
    ax.set_xlim((-0.5,len(organism_names)-0.5))
    output_folder = elem_plot_folder
    filename = "length barplot "+domain_list[0]
    plt.tight_layout()
    plt.savefig(output_folder+filename+".png")
    plt.savefig(output_folder+filename+".pdf")