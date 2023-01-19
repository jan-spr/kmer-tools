#Functions for plotting various types of heatmaps:
# - symmetric and asymmetric correlation heatmaps
# - density heatmaps



from cmath import nan
import sys

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transf
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FixedLocator
import seaborn as sns
from tqdm import tqdm

from kmer_tools.plot.get_spec_data_lib import *

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data_strand/"
heatmap_folder  = "../../plots/kmer_correlation/heatmaps/domain_heatmaps/"
heatmap_combined_folder  = "../../plots/kmer_correlation/heatmaps/domain_heatmaps/all_orgas/"
heatmap_repeats_folder  = "../../plots/kmer_correlation/heatmaps/domain_heatmaps/repeats/"
graph_folder    = "../../plots/kmer_correlation/k_correlation_orga/"
coverage_folder = "../../plots/elem_data/coverage/"
coverage_hm_folder = "../../plots/elem_data/coverage_hm/"
func_hm_folder = "../../plots/kmer_correlation/heatmaps/all_orgas/"
elem_plot_folder = "../../plots/elem_data/"

separator_width=4

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

def lerp(a,b,x):
    return a + (b-a)*x

def check_spec_empty(spec):
    kmers = spec.get_kmers()
    counts = [spec.get_count(kmer) for kmer in kmers]
    if np.sum(counts)== 0:  return False
    else:                   return True

def empty_filter(spectra):
    empty_arr = [check_spec_empty(spec) for spec in spectra]
    empty_ind = [i for i,bool in zip(np.arange(len(empty_arr)), empty_arr) if bool]
    spectra_filtered = [spec for spec, bool in zip(spectra, empty_arr) if bool]
    return empty_ind, spectra_filtered

def fit_cbar_height(hm_ax,cbar_ax,outline=True):
    pos_main = hm_ax.get_position().get_points() #[[x0, y0], [x1, y1]]
    pos_cbar = cbar_ax.get_position().get_points()
    bottom = pos_main[0][1]
    top = pos_main[1][1]
    left   = pos_cbar[0][0]
    right  = pos_cbar[1][0]
    cbar_ax.set_position(transf.Bbox([[left,bottom],[right,top]]))
    if outline:
        for spine in cbar_ax.spines.values():
            spine.set(visible=True, edgecolor="black")


def add_domain_split(orga_names_x,orga_names_y,ax,split_x=True,split_y=True,labels=False):
    xticks = []
    yticks = []
    for labels,ticks in zip([orga_names_x,orga_names_y],[xticks,yticks]):
        prev_domain = None
        for i, domain in enumerate([kmer_tools.download_genbank.name_to_domain(name) for name in labels]):
            #print(domain, " - ",labels[i])
            if domain != prev_domain:
                #print("Switch")
                ticks.append(i)
                prev_domain = domain
    if(split_x):
        ax.xaxis.set_minor_locator(FixedLocator(xticks[1:]))
        ax.vlines(xticks[1:],0,len(orga_names_y),color="white",lw=separator_width)
    if(split_y):
        ax.yaxis.set_minor_locator(FixedLocator(yticks[1:]))
        ax.hlines(yticks[1:],0,len(orga_names_x),color="white",lw=separator_width)

def add_sublabel_split(orga_names_x,orga_names_y,fig,ax,split_x=True,split_y=True,show_labels=False,class_split=False,mark_wrap=False):
    if class_split:     label_func = kmer_tools.download_genbank.name_to_class
    else:               label_func = kmer_tools.download_genbank.name_to_domain
    #label_func = kmer_tools.download_genbank.name_to_shortlabel # edit for nucleomorph plots
    xticks = []
    yticks = []
    xlabels = []
    ylabels = []

    for orga_labels,ticks,ax_labels in zip([orga_names_x,orga_names_y],[xticks,yticks],[xlabels,ylabels]):
        prev_domain = None
        for i, domain in enumerate([label_func(name) for name in orga_labels]):
            #print(domain, " - ",orga_labels[i])
            if domain != prev_domain:
                #print("Switch")
                ticks.append(i)
                prev_domain = domain
                ax_labels.append(domain)
        ticks.append(i + 1)
    print(xticks,"/n",yticks)
    xticks_n=[]
    xlabels_n=[]
    for i in range(len(xticks)-1):
        if np.abs(xticks[i]- xticks[i+1]) >2: 
            xticks_n.append(xticks[i])
            xlabels_n.append(xlabels[i])
    xticks_n.append(xticks[-1])
    xticks=xticks_n
    xlabels=xlabels_n
    print(xticks_n,xlabels)
    if(split_x):
        if not mark_wrap: ax.xaxis.set_minor_locator(FixedLocator(xticks[1:-1]))
        ax.vlines(xticks[1:-1],0,len(orga_names_y),color="white",lw=separator_width)
        if show_labels:
            ax.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(xticks[:-1], xticks[1:])]))
            print(xlabels)
            ax.set_xticklabels(xlabels, rotation=0)
    if(split_y):
        if not mark_wrap: ax.yaxis.set_minor_locator(FixedLocator(yticks[1:-1]))
        ax.hlines(yticks[1:-1],0,len(orga_names_x),color="white",lw=separator_width)
        if show_labels:
            ax.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(xticks[:-1], xticks[1:])]))
            ax.set_yticklabels(ylabels, rotation=90)

    lenx = len(orga_names_x)
    leny = len(orga_names_y)
    add_title=True
    if add_title:
        pos_main = ax.get_position().get_points() #[[x0, y0], [x1, y1]]
        print("pos_main: ",pos_main)
        x_lims = [pos_main[0][0],pos_main[1][0]]
        y_lims = [pos_main[1][1],pos_main[0][1]]
        x_pos = pos_main[0][0]
        y_pos = pos_main[1][1]
        x_stops = xticks
        y_stops = yticks
        print("x_lims,y_pos,x_stops,func_names1: ", x_lims,y_pos,x_stops,xlabels)
        print("xticks: ", xticks)
        
        offset=0.11
        for i in range(len(xticks)-1):
            left_n = x_stops[i]
            right_n = x_stops[i+1]
            right= lerp(x_lims[0],x_lims[1],right_n/lenx)
            left= lerp(x_lims[0],x_lims[1],left_n/lenx)
            bottom = top = y_pos
            width=right-left
            height=0.0
            
            rect=[left,bottom,width,height]
            print("left,bottom,width,height: ", rect)
            #title_ax = fig.add_axes(rect)
            #print("title_ax pos: ", title_ax.get_position().get_points())
            #title_ax.set_title(func_names1[i])
            if (np.abs(left_n-right_n)>3 and split_x):
                print ("adding text: ",xlabels[i] )
                ax.text(x=(right+left)/2,y=y_pos+offset,s=xlabels[i],transform = fig.transFigure,fontsize=17,ha='center',va='center')
                pass
        for i in range(len(yticks)-1):
            top_n = y_stops[i]
            bottom_n = y_stops[i+1]
            top= lerp(y_lims[0],y_lims[1],top_n/leny)
            bottom= lerp(y_lims[0],y_lims[1],bottom_n/leny)
            left = right = x_pos
            height=bottom-top
            width=right-left
            height=0.0
            rect=[left,bottom,width,height]
            print("left,bottom,width,height: ", rect)
            #title_ax = fig.add_axes(rect)
            #print("title_ax pos: ", title_ax.get_position().get_points())
            #title_ax.set_title(func_names1[i])
            if (np.abs(top_n-bottom_n)>3 and split_y):
                ax.text(y=(top+bottom)/2,x=x_pos-offset,s=ylabels[i],transform = fig.transFigure,fontsize=17,rotation='vertical',ha='center',va='center')

    if(split_x):
        ax.xaxis.set_minor_locator(FixedLocator(xticks[:-1])) 
    if(split_y):
        ax.yaxis.set_minor_locator(FixedLocator(yticks[:-1]))

def add_wrap_split(fig,ax,orga_names_x,orga_names_y,func_names1,func_names2,split_x=True,split_y=True,show_labels=False,class_split=False,mark_wrap=False,add_title=True):

    xticks = []
    yticks = []
    xlabels = []
    ylabels = []
    lenx = len(orga_names_x)
    leny = len(orga_names_y)

    for orga_labels,ticks,ax_labels in zip([orga_names_x,orga_names_y],[xticks,yticks],[xlabels,ylabels]):
        prev_name_location = 0
        name_location=0
        for i, name_location in enumerate([
                kmer_tools.download_genbank.orga_names_sorted.Sorted_Names_All_normal.index(name) 
                for name in orga_labels
                ]):
            #print(domain, " - ",orga_labels[i])
            if name_location < prev_name_location:
                #print("Switch")
                ticks.append(i)
            prev_name_location=name_location
        ticks.append(i + 1)
    print(xticks,"/n",yticks)
    if add_title:
        pos_main = ax.get_position().get_points() #[[x0, y0], [x1, y1]]
        print("pos_main: ",pos_main)
        x_lims = [pos_main[0][0],pos_main[1][0]]
        y_lims = [pos_main[1][1],pos_main[0][1]]
        x_pos = pos_main[0][0]
        y_pos = pos_main[1][1]
        x_stops = [0] + xticks
        y_stops = [0] + yticks
        print("x_lims,y_pos,x_stops,func_names1: ", x_lims,y_pos,x_stops,func_names1)
        offset=0.085
        for i in range(len(xticks)):
            left_n = x_stops[i]
            right_n = x_stops[i+1]
            right= lerp(x_lims[0],x_lims[1],right_n/lenx)
            left= lerp(x_lims[0],x_lims[1],left_n/lenx)
            bottom = top = y_pos
            width=right-left
            height=0.0
            
            rect=[left,bottom,width,height]
            print("left,bottom,width,height: ", rect)
            #title_ax = fig.add_axes(rect)
            #print("title_ax pos: ", title_ax.get_position().get_points())
            #title_ax.set_title(func_names1[i])
            ax.text(x=(right+left)/2,y=y_pos+offset,s=func_names1[i],transform = fig.transFigure,fontsize='xx-large')
        for i in range(len(yticks)):
            top_n = y_stops[i]
            bottom_n = y_stops[i+1]
            top= lerp(y_lims[0],y_lims[1],top_n/leny)
            bottom= lerp(y_lims[0],y_lims[1],bottom_n/leny)
            left = right = x_pos
            height=bottom-top
            height=0.0
            rect=[left,bottom,width,height]
            print("left,bottom,width,height: ", rect)
            #title_ax = fig.add_axes(rect)
            #print("title_ax pos: ", title_ax.get_position().get_points())
            #title_ax.set_title(func_names1[i])
            ax.text(y=(top+bottom)/2,x=x_pos-offset,s=func_names2[i],transform = fig.transFigure,fontsize='xx-large',rotation='vertical')

    if(split_x):
        ax.xaxis.set_minor_locator(FixedLocator(xticks[:-1])) 
    if(split_y):
        ax.yaxis.set_minor_locator(FixedLocator(yticks[:-1]))


def plot_orga_heatmap(matrix, organism_names, output_file, vmin, vmax, figsize=(15,15), x_label=None, y_label=None):
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.heatmap(data=matrix,  vmin=-1., vmax=1., cbar=False, cmap='rainbow', square=True)
    ax.xaxis.tick_top() # x axis on top
    ax.xaxis.set_label_position('top')
    #print(orga_names)
    ax.set_yticks(np.arange(len(organism_names))+0.5)
    ax.set_xticks(np.arange(len(organism_names))+0.5)
    ax.set_xticklabels(organism_names,rotation=90)
    ax.set_yticklabels(organism_names,rotation=0)

    if(x_label is not None): ax.set_xlabel(x_label)
    if(y_label is not None): ax.set_ylabel(y_label)

    ax.tick_params(axis='both', which='minor', labelsize=2)
    #plt.yticks(rotation=90)
    im = ax.imshow(matrix,cmap='rainbow',vmin=-1., vmax=1.,)
    cbar = fig.colorbar(im,fraction=0.046, pad=0.04,drawedges=False)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_visible(True)
        ax.spines[axis].set_color('black')
        ax.spines[axis].set_visible(True)
        ax.spines[axis].set_color('black')    

    ax.figure.tight_layout()

    #plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)
    
    print("saving plot: "+output_file)
    plt.savefig(output_file)

def plot_heatmap(matrix,figsize,labels_y,labels_x,output_file=None,cmap='rainbow',v_min=None,v_max=None,square=True,log=False, x_label=None, y_label=None,split_x=False,split_y=False,split_class=False,mark_wrap=True):
    fig, axs = plt.subplots(1,2,figsize=figsize,gridspec_kw={'width_ratios' :[1, 0.05]})
    #ax_1 = axs[0]
    cbar_ax = axs[1]
    sns.set(rc = {'figure.figsize':figsize})
    norm = Normalize(vmin=v_min,vmax=v_max)
    if (log):   norm = LogNorm()
    ax = sns.heatmap(data=matrix, ax=axs[0],vmin=v_min, vmax=v_max, cbar=True, cbar_ax=axs[1], square=square, cmap=cmap, norm=norm,edgecolor="black")
    ax.xaxis.tick_top() # x axis on top
    ax.xaxis.set_label_position('top')
    ax.yaxis.tick_left() # x axis on top
    ax.yaxis.set_label_position('left')
    #print(orga_names)
    #print(matrix.shape)
    ax.set_yticks(np.arange(len(labels_y))+0.5)
    ax.set_xticks(np.arange(len(labels_x))+0.5)
    if False:
        #labels_x = [label.split(" ")[0]+" "+label.split(" ")[1][0] for label in  labels_x]
        ax.set_xticklabels([label.split(" ")[0]+" "+label.split(" ")[1][0]+"." for label in  labels_x],rotation=90)
        ax.set_yticklabels([label.split(" ")[0]+" "+label.split(" ")[1][0]+"." for label in  labels_y],rotation=0)
    else:
        
        ax.set_xticklabels(labels_x,rotation=90)
        ax.set_yticklabels(labels_y,rotation=0)
    

    if (False):
        if(x_label is not None): ax.set_xlabel(x_label,fontsize=20)
        if(y_label is not None): ax.set_ylabel(y_label,fontsize=20)

    
    plt.tight_layout()
    orga_short=True
    if (orga_short):
        ax.set_yticklabels([name.split(" ")[0]+" "+name.split(" ")[1][0]+"." for name in labels_y],rotation=0)
        #ax.set_xticklabels([name.split(" ")[0]+" "+name.split(" ")[1][0]+"." for name in labels_x],rotation=90)
    
    add_sublabel_split(labels_x,labels_y,fig,axs[0],split_x=split_x,split_y=split_y,class_split=split_class,mark_wrap=mark_wrap)
    
    


    ax.tick_params(axis='both', which='minor', length=80,labelsize=2,color="purple",width=separator_width)
    #plt.yticks(rotation=90)
    #im = ax.imshow(matrix,cmap='rainbow',vmin=-1., vmax=1.,)
    #cbar = fig.colorbar(im,fraction=0.046, pad=0.04,drawedges=False)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_visible(True)
        ax.spines[axis].set_color('black')
        ax.spines[axis].set_visible(True)
        ax.spines[axis].set_color('black')    
    #plt.tight_layout()
    

    if False: 
        add_wrap_split(labels_x,labels_y,ax,split_x=split_x,split_y=split_y,class_split=split_class,mark_wrap=True)
        #print(func_names1,func_names2)

    # changing height of cbar to match heatmap
    fit_cbar_height(ax,cbar_ax)
    cbar_ax.set_frame_on(True)

    #plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)
    return fig,axs
    

def plot_orga_heatmap_asym(matrix,names_left,names_top,output_file,v_min=None,v_max=None,cmap="rainbow",square=True,x_label=None, y_label=None,figsize=(13,13),split_class=False,mark_wrap=False):
    label_y=names_left
    label_x=names_top
    return plot_heatmap(matrix,figsize,label_y,label_x,output_file,cmap,v_min=v_min,v_max=v_max,square=square,log=False,x_label=x_label, y_label=y_label,split_x=True,split_y=True,split_class=split_class,mark_wrap=mark_wrap)



def plot_corr_chromosomal(domain,organism_names,k_values,function_names):
    for i,func_name in tqdm(enumerate(function_names)):
        print(func_name)
        for k in k_values:
            matrix = get_orga_kmer_heatmap(organism_names,func_name,function_names[0],k,domain=domain)
            
            # Save Matrix
            # matrix.save('./data/matrix/chromosomal/'+domain+'_'+func_name+'_k='+str(k)+'_correlation.matrix')
            # Display Matrix
            print(matrix.matrix)
            output_file = heatmap_folder+domain+'/'+func_name+' '+domain+' k='+str(k)+'.kmer.hm.png'
            plot_orga_heatmap(matrix.matrix,organism_names,output_file,vmin=-1,vmax=1)

def plot_corr_all_orgas(domain_dict,organism_names,k_values,function_names,figsize=(25,25),repeatmask=False,sorted=True,filter_empty=False):
    if (sorted): organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
    for i,func_name in tqdm(enumerate(function_names)):
        print(func_name)
        for k in k_values:
            matrix = get_orga_kmer_heatmap(organism_names,func_name,function_names[0],k,domain_dict=domain_dict)
            
            # Save Matrix
            # matrix.save('./data/matrix/chromosomal/'+domain+'_'+func_name+'_k='+str(k)+'_correlation.matrix')
            # Display Matrix
            print(matrix.matrix)
            output_file = heatmap_combined_folder+func_name+' k='+str(k)+'.kmer.hm.png'
            if (repeatmask): output_file = heatmap_repeats_folder+func_name+' k='+str(k)+'.kmer.hm.png'
            #print (organism_names)
            #print([kmer_tools.download_genbank.orga_name_normal_dict[name] for name in organism_names])
            plot_orga_heatmap(matrix.matrix,[kmer_tools.download_genbank.orga_name_normal_dict[name] for name in organism_names],output_file,vmin=-1,vmax=1,figsize=figsize)

def plot_coverage_bar(organis_names,elem_names,domain,log = True):
    sorted_idx=[0,6,1,3,5,2,7,8,9,10,11,12,16,17,13,4,14,15,18,19]
    sorted_idx=np.arange(len(elem_names))
    num=len(elem_names)
    j=0

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(10,7))
    
    xpos = np.arange(len(sorted_idx))+0.5
    
    values = []
    for elem in tqdm(elem_names):
        cov = []
        for orga_name in organis_names:
            cov.append(kmer_tools.search_elem_stat.read_prop_orga(orga_name,"coverage_density",domain,elem))
        values.append(np.mean(cov))
    print(len(values),len(sorted_idx))
    values = np.array(values)[sorted_idx]

    bars = ax.bar(xpos, values,align='center')
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    ax.set_ylabel("coverage")
    #ax[j].set_title(property_titles[j]+' of '+function_titles[i]+' in archaea')

    # Label with specially formatted floats
    #ax.bar_label(hbars, fmt='%.2f')
    ax.set_xlim((0,len(sorted_idx)))  # adjust xlim to fit labels

    #adjustments to first and last plots
    ax.set_title('average genome coverage in '+domain)
    fig.subplots_adjust(bottom=0.29)
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True) # labels along the bottom edge are off
    ax.set_xticks(xpos)
    ax.set_xticklabels(np.array(elem_names)[sorted_idx], rotation='vertical')
    if (log): ax.set_yscale('log')
    #ax[num-1].set_ylim((0,1))
    plt.savefig(coverage_folder+domain+'_coverage_repeat.png')

def plot_coverage_heatmap(organism_names,elem_names,domain,log = True,square=True,sorted=True, label = None):
    matrix =[]
    if (sorted): organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
    print(organism_names)
    for elem in tqdm(elem_names):
        cov = []
        for orga_name in organism_names:
            cov.append(kmer_tools.search_elem_stat.read_prop_orga(orga_name,"coverage_density",domain,elem))
        #print(len(cov))
        matrix.append(cov)
    #print (matrix)
    if label is not None: label = "_" + label
    else: label = ""
    output_file = coverage_hm_folder +domain+'_coverage' + label +'.png'
    orga_labels = [kmer_tools.download_genbank.orga_name_normal(name) for name in organism_names]
    fig,axs = plot_heatmap(np.transpose(matrix)+1e-5,figsize=(7,4),labels_y=orga_labels,labels_x=elem_names,output_file=output_file,cmap='jet',square=square,log=log,split_y=True)
    plt.savefig(output_file)

def plot_coverage_heatmap_list(organism_names,elem_names,element_labels,domain_list,log = True,square=True,sorted=True, label = None):
    matrix =[]
    if (sorted): organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
    print(organism_names)
    for elem in tqdm(elem_names):
        cov = []
        for i,orga_name in enumerate(organism_names):
            cov.append(kmer_tools.search_elem_stat.read_prop_orga(orga_name,"coverage_density",domain_list[i],elem))
        #print(len(cov))
        matrix.append(cov)
    #print (matrix)
    matrix = np.transpose(matrix)
    #filter empty
    filter_inds=[]
    new_mat=[]
    for i,row in enumerate(matrix):
        if np.sum(row)>10e-11:
            filter_inds.append(i)
            new_mat.append(row)
    matrix=new_mat
    
    if label is not None: label = "_" + label
    else: label = ""
    output_file = coverage_hm_folder +'_coverage' + label +'.png'
    organism_names=[organism_names[i] for i in filter_inds]
    orga_labels = [kmer_tools.download_genbank.orga_name_normal(name) for name in organism_names]
    print(np.array(matrix).shape)
    print(len(orga_labels),len(elem_names))
    fig,axs = plot_heatmap(np.array(matrix)+1e-5,figsize=(4,6),labels_y=orga_labels,labels_x=element_labels,output_file=output_file,cmap='jet',square=square,log=log,split_y=True,split_class=True)
    plt.savefig(output_file)

def plot_prop_heatmap(organism_names,elem_names,domain,prpty,folder,log = True,square=True,sorted=True):
    matrix =[]
    if (sorted): organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
    print(organism_names)
    for elem in tqdm(elem_names):
        cov = []
        for orga_name in organism_names:
            cov.append(kmer_tools.search_elem_stat.read_prop_orga(orga_name,prpty,domain,elem))
        #print(len(cov))
        matrix.append(cov)
    #print (matrix)
    output_file = elem_plot_folder + folder + '/' + domain + '_coverage.png'
    orga_labels = [kmer_tools.download_genbank.orga_name_normal_dict[name] for name in organism_names]
    plot_heatmap(np.transpose(matrix)+1e-5,figsize=(8,10),label_y=orga_labels,label_x=elem_names,output_file=output_file,cmap='jet',square=square)

def plot_dimer_hm(organisms,domain=None,domain_dict=None,func_name="intergenics",sorted_gc=False):
    print(func_name)
    organisms_sorted = organisms
    if(sorted_gc): organisms_sorted = kmer_tools.plot.sort_GC(organisms)
    dimer_words = ["AA",	"TT",	"AT",	"TA",	"AC",	"CA",	"TG",	"GT",	"AG",	"GA",	"TC",	"CT",	"CC",	"GG",	"GC",	"CG"]
    dimer_rows = []
    for orga in tqdm(organisms_sorted):
        row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,func_name=func_name,orga_name = orga.name)
        counts_norm,_ = kmer_tools.search_kmer.dimer_count_normalized(row,dimer_words)
        dimer_rows.append(counts_norm)
    orga_labels = [kmer_tools.download_genbank.orga_name_normal_dict[orga.name] for orga in organisms_sorted]
    output_folder = elem_plot_folder + "dimer_hm" + "/"
    output_name = domain + "_" + func_name + "_" + "dimer_hm.png"
    #print (dimer_rows)
    #print (dimer_rows[0])
    plot_heatmap(matrix=dimer_rows,figsize=(7,10),label_y=orga_labels,label_x=dimer_words,output_file = output_folder + output_name,cmap="rainbow",v_min=None,v_max=None,square=True)

def plot_func_heatmap(orga_names,func_name1,func_name2,k_value,domain_list=None,domain_dict=None,sorted=True,func_label="",figsize=(15,15),filter_empty=False,filer_Nan=False,v_min=None,v_max=None,split_class=False,domain=None):
    domain_degree = np.count_nonzero([(entry is not None) for entry in [domain_list,domain_dict,domain]])
    if(domain_degree != 1): print("!!domains over- or under-defined!! degree: ",domain_degree)
    if (sorted): orga_names = kmer_tools.download_genbank.sort_orga_names(orga_names,dicto=True)

    spectra1= get_orga_kmer_spectra(orga_names, func_name1, "genes", k_value, domain=domain, domain_dict=domain_dict, domain_list=domain_list, fix_empty = True)
    if (func_name1 == func_name2):  
            spectra2 = spectra1
    else:   spectra2= get_orga_kmer_spectra(orga_names, func_name2, "genes", k_value, domain=domain, domain_dict=domain_dict, domain_list=domain_list, fix_empty = True)

    filter_ind1=np.arange(len(spectra1))
    filter_ind2=np.arange(len(spectra2))
    if filter_empty:
        filter_ind1,spectra1 = empty_filter(spectra1)
        filter_ind2,spectra2 = empty_filter(spectra2)
    matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra1=spectra1,spectra2=spectra2)
    normal_names = [kmer_tools.download_genbank.orga_name_normal(name) for name in orga_names]

    names1 = [normal_names[i] for i in filter_ind1]
    names2 = [normal_names[i] for i in filter_ind2]

    otuput_folder = func_hm_folder
    if func_label=="": func_label = func_name1 +" - "+func_name2+" "+"k="+str(k_value)
    filepath=otuput_folder+func_label+"_hm.png"
    if filer_Nan:
        nan_x=~np.isnan(matrix.matrix).all(axis=0)
        nan_y=~np.isnan(matrix.matrix).all(axis=1)
        print(nan_x,nan_y)
        print(np.array(matrix.matrix)[:,nan_x])
        print(np.array(matrix.matrix)[nan_y])
        print(np.array(matrix.matrix)[nan_y,nan_y])
        filter_both = [x_val and y_val for x_val,y_val in zip(nan_x,nan_y)]
        print(filter_both)
    fig,axs =plot_orga_heatmap_asym(matrix.matrix,names1,names2,output_file=filepath,square=True,x_label=func_name2,y_label=func_name1,figsize=figsize,v_min=v_min,v_max=v_max,split_class=split_class)
    plt.savefig(otuput_folder+func_label+"_hm.png")
    plt.savefig(otuput_folder+func_label+"_hm.pdf")

def plot_funcs_mul_heatmap(orga_names,func_names1,func_names2,k_value,domain_list=None,domain_dict=None,sorted=True,func_label="",figsize=(15,15),filter_empty=False,filer_Nan=False,v_min=None,v_max=None,split_class=False,domain=None,mark_wrap=False):
    domain_degree = np.count_nonzero([(entry is not None) for entry in [domain_list,domain_dict,domain]])
    if(domain_degree != 1): print("!!domains over- or under-defined!! degree: ",domain_degree)
    if (sorted): orga_names = kmer_tools.download_genbank.sort_orga_names(orga_names,dicto=True)
    spectra1 = []
    spectra2 = []
    
    for func_name in func_names1:
        spectra1 += get_orga_kmer_spectra(orga_names, func_name, "genes", k_value, domain=domain, domain_dict=domain_dict, domain_list=domain_list, fix_empty = True)
    if func_names2 == func_names1:
        spectra2 = spectra1
    else:
        for func_name in func_names2:    
            spectra2 += get_orga_kmer_spectra(orga_names, func_name, "genes", k_value, domain=domain, domain_dict=domain_dict, domain_list=domain_list, fix_empty = True)

    filter_ind1=np.arange(len(spectra1))
    filter_ind2=np.arange(len(spectra2))
    if filter_empty:
        filter_ind1,spectra1 = empty_filter(spectra1)
        filter_ind2,spectra2 = empty_filter(spectra2)
    matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra1=spectra1,spectra2=spectra2)
    normal_names = [kmer_tools.download_genbank.orga_name_normal(name) for name in orga_names]

    names1 = [normal_names[i % len(normal_names)] for i in filter_ind1]
    names2 = [normal_names[i % len(normal_names)] for i in filter_ind2]

    otuput_folder = func_hm_folder
    if func_label=="": func_label = " - ".join(func_names1) + " " + "k="+str(k_value)
    filepath=otuput_folder+func_label+"_hm.png"
    if filer_Nan:
        nan_x=~np.isnan(matrix.matrix).all(axis=0)
        nan_y=~np.isnan(matrix.matrix).all(axis=1)
        print(nan_x,nan_y)
        print(np.array(matrix.matrix)[:,nan_x])
        print(np.array(matrix.matrix)[nan_y])
        print(np.array(matrix.matrix)[nan_y,nan_y])
        filter_both = [x_val and y_val for x_val,y_val in zip(nan_x,nan_y)]
        print(filter_both)
    fig,axs = plot_orga_heatmap_asym(matrix.matrix,names1,names2,output_file=filepath,square=True,figsize=figsize,v_min=v_min,v_max=v_max,split_class=split_class,mark_wrap=mark_wrap)

    if mark_wrap: 
        add_wrap_split(fig,axs[0],names1,names2,func_names1,func_names2,split_x=True,split_y=True,class_split=split_class,mark_wrap=True,add_title=True)

    print("saving plot: "+filepath)
    plt.savefig(filepath,dpi=400)
    plt.savefig(filepath[0:-3]+"pdf",dpi=400)


