from posixpath import split
import sys

from numpy import square
from joblib import Parallel, delayed
import scipy as sp

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('../..')   
import kmer_tools

from kmer_tools.download_genbank.orga_names_sorted import *
from kmer_tools.download_genbank.orga_names_labeling import *
from kmer_tools.search_kmer.repeat_ids import *
from plot_hm_lib import *
from plot_GC_lib import *
from plot_dimer_bar_lib import *
from plot_dimer_hm_lib import *
from plot_dimer_graph_lib import *
from kmer_tools.search_kmer.count_func_lib import *
from plot_graph_lib import *
from plot_box_lib import *

# ------------------------------------------------------------------------------------
# -------------------------------- PREREQUISITES    ----------------------------------
# ------------------------------------------------------------------------------------

functions=[
    Oligo.File.read_genes,
    Oligo.File.read_exons,
    Oligo.File.read_introns,
    Oligo.File.read_cds,
    Oligo.File.read_mRNAs,
    Oligo.File.read_intergenics,
    Oligo.File.read_protein_coding_genes,
    Oligo.File.read_pseudo_genes,
    Oligo.File.read_pseudo_genes_exons,
    Oligo.File.read_pseudo_genes_introns,
    Oligo.File.read_non_protein_coding_genes,
    Oligo.File.read_ncRNAs,
    Oligo.File.read_lncRNAs,
    Oligo.File.read_miRNAs,
    Oligo.File.read_snoRNAs,
    Oligo.File.read_snRNAs,
    Oligo.File.read_tRNAs,
    Oligo.File.read_rRNAs,
    Oligo.File.read_precursor_RNAs,
    Oligo.File.read_gaps
]

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


merge_func_names = [
    "pc_exons",
    "pc_introns",
    "npc_exons",
    "npc_introns",
    "pseudo_exons",
    "pseudo_introns",
]

function_name_lim1=[
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
    #'precursor_RNAs',
    #'gaps',
]

function_name_lim_merge=[
    'genes',
    'exons',
    'introns',
    'cds',
    'mRNAs',
    'protein_coding_genes',
    "pc_exons",
    "pc_introns",
    'pseudo_genes',
    "pseudo_exons",
    "pseudo_introns",
    'non_protein_coding_genes',
    "npc_exons",
    "npc_introns",
    'ncRNAs',
    'tRNAs',
    'rRNAs',
    'intergenics',
]

function_name_lim=[
    'genes',
    'exons',
    'introns',
    'cds',
    #'mRNAs',
    'intergenics',
    'protein_coding_genes',
    'pseudo_genes',
    #'pseudo_genes_exons',
    #'pseudo_genes_introns',
    'non_protein_coding_genes',
    'ncRNAs',
    #'lncRNAs',
    #'miRNAs',
    #'snoRNAs',
    #'snRNAs',
    'tRNAs',
    'rRNAs',
    #'precursor_RNAs',
    #'gaps',
]

function_name_4=[
    'intergenics',
    'protein_coding_genes',
    'pseudo_genes',
    'non_protein_coding_genes',
]
function_name_title_4=[
    'intergenics',
    'pc genes',
    'pseudo genes',
    'npc genes',
]
function_name_exon_4=[
    'intergenics',
    'pc_exons',
    'pseudo_exons',
    'npc_exons',
]

function_name_RNA=[
    #'ncRNAs',
    'lncRNAs',
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'tRNAs',
    'rRNAs',
    #'precursor_RNAs',
]
function_name_RNA2=[
    #'ncRNAs',
    'lncRNAs',
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'tRNAs',
    'rRNAs',
    #'precursor_RNAs',
]

function_pairs_RNA=[
    #['ncRNAs','mRNAs'],
    ['ncRNAs','lncRNAs'],
    ['ncRNAs','miRNAs'],
    ['ncRNAs','snoRNAs'],
    ['ncRNAs','snRNAs'],
    ['ncRNAs','rRNAs'],
    ['ncRNAs','tRNAs'],
]
function_pairs_RNA2=[
    ['tRNAs','rRNAs'],
    ['tRNAs','snRNAs'],
    ['tRNAs','snoRNAs'],
    ['rRNAs','snRNAs'],
    ['rRNAs','snoRNAs'],
    ['snRNAs','snoRNAs'],
]

function_pairs_miRNA=[
    ['miRNAs','lncRNAs'],
    ['miRNAs','snoRNAs'],
    ['miRNAs','snRNAs'],
    ['miRNAs','rRNAs'],
    ['miRNAs','tRNAs'],
]


function_names_repeat=[
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

function_pairs = [
  ['protein_coding_genes','pseudo_genes'],
  ['intergenics','protein_coding_genes'],
  ['intergenics','pseudo_genes'],
  ['protein_coding_genes','non_protein_coding_genes'],
  ['non_protein_coding_genes','pseudo_genes'],
  ['intergenics','non_protein_coding_genes'],
]

function_pair_titles = [
  ['pc','pseudo'],  
  ['intergenics','pc'],
  ['intergenics','pseudo'],
  ['npc','pc'],
  ['npc','pseudo'],
  ['npc','intergenics'],
]

function_pairs = [
  ['intergenics','protein_coding_genes'],
  ['intergenics','non_protein_coding_genes'],
  ['intergenics','pseudo_genes'],
  ['protein_coding_genes','non_protein_coding_genes'],
  ['protein_coding_genes','pseudo_genes'],
  ['non_protein_coding_genes','pseudo_genes'],
]

function_pair_titles = [
  ['intergenics','pc'],
  ['intergenics','npc'], 
  ['intergenics','pseudo'],
  ['pc','npc'],
  ['pc','pseudo'],
  ['npc','pseudo'],
]

function_pairs_exons = [
  ['intergenics','pc_exons'],
  ['intergenics','npc_exons'],
  ['intergenics','pseudo_exons'],
  ['pc_exons','npc_exons'],
  ['pc_exons','pseudo_exons'],
  ['npc_exons','pseudo_exons'],
]

function_pairs_introns = [
  ['intergenics','pc_introns'],
  ['intergenics','npc_introns'],
  ['intergenics','pseudo_introns'],
  ['pc_introns','npc_introns'],
  ['pc_introns','pseudo_introns'],
  ['npc_introns','pseudo_introns'],
]
function_pairs_intergenics = [
  ['intergenics','pc_introns'],
  ['intergenics','npc_introns'],
  ['intergenics','pseudo_introns'],
  ['intergenics','pc_exons'],
  ['intergenics','npc_exons'],
  ['intergenics','pseudo_exons'],
]

covariance_funcs = [
    "intergenics",
    "pc_introns",
    "npc_introns",
    "pseudo_introns",
    "npc_exons",
    "pseudo_exons",
    "pc_exons",
    "cds"
    ]

#domains = ["animalia","animalia2","embryophyta","fungi","protista"]
#for domain in domains:
#    print(domain+": ",len(orgas[domain]),[orga.name for orga in orgas[domain]])
#orgas["animalia"] += orgas["animalia2"]

domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids","nucleomorph"]
orgas = {}
domain_dict={}
for domain in domains:
    orgas[domain] = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")

orga_names_all=[]
for domain in domains[0:4]:
    orgas[domain] = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga in orgas[domain]:
        domain_dict[orga.name]=domain
    orga_names_all += [orga.name for orga in orgas[domain] ]
    

repeat_names = [
    "Homo sapiens",
    "Pan troglodytes",
    "Mus musculus",
    "Rattus norvegicus",
    "Oryctolagus cuniculus",
    "Canis lupus familiaris",
    "Equus caballus",
    "Capra hircus",
    "Bos taurus",
    "Sus scrofa",
    "Monodelphis domestica",

    "Gorilla gorilla gorilla",
    "Pongo abelii",
    "Macaca mulatta",
    "Microtus ochrogaster",
    "Ochotona princeps",
    "Microcebus murinus",
    
    "Ficedula albicollis",
    "Strigops habroptila",
    "Falco rusticolus",
    "Calypte anna",
]

function_pairs_var = [
  ['pc_exons','npc_exons'],
  ['pc_exons','pseudo_exons'],
  ['npc_exons','pseudo_exons'],
]



function_pair_titles_var = [
  ['pc-exons','npc-exons'],
  ['pc-exons','pseudo-exons'],
  ['npc-exons','pseudo-exons'],
]

# ------------------------------------------------------------------------------------
# ----------------------------- ORGANISMS NORMAL    ----------------------------------
# ------------------------------------------------------------------------------------

function_name_cds=[
    'intergenics',
    'introns',
    'exons',
    'cds',
]

function_name_cds_sub=[
    'cds',
    'pc_exons',
    'pseudo_exons',
    'npc_exons',
]

function_names_hm=[
    'genes',
    'protein_coding_genes',
    #'pseudo_genes',
    'non_protein_coding_genes',
    'exons',
    'pc_exons',
    'npc_exons',
    #'pseudo_exons',
    'cds',
    'mRNAs',
    'introns',
    'pc_introns',
    'npc_introns',
    #'pseudo_introns',
    'intergenics',
    'ncRNAs',
    #'lncRNAs',
    #'miRNAs',
    #'snoRNAs',
    #'snRNAs',
    'tRNAs',
    'rRNAs',
    #'precursor_RNAs',
    #'gaps',
]

for domain in domains[4:5]:
    #print(domain)
    #print(domain+": ",len(orgas[domain]),[orga.name for orga in orgas[domain]])
    organism_names = [orga.name for orga in orgas[domain]]
    organism_names = kmer_tools.download_genbank.sort_mito_names(organism_names)
    #print(organism_names)
    #organisms = kmer_tools.download_genbank.names_to_orgas(organism_names,domain=domain)
    #plot_corr_chromosomal(domain,organism_names,[2],function_names[1:2])
    #plot_corr_k(domain,organism_names,[1,2,3,4,5,6],function_name_lim,colors)
    #plot_corr_k_categ(domain,orgas[domain],organism_names,[1,2,3,4,5,6],function_pairs_var,function_pair_titles_var,colors,file_label="ex_in")
    #plot_corr_k_categ(domain,orgas[domain],organism_names,[1,2,3,4,5,6],function_pairs_RNA,function_pairs_RNA,colors,file_label="ncRNAs")
    #plot_corr_k_categ(domain,orgas[domain],organism_names,[1,2,3,4,5,6],function_pairs_intergenics,function_pairs_intergenics,colors,file_label="intergenic - pc npc pseudo")
    #plot_corr_k_categ(domain,orgas[domain],organism_names,[1,2,3,4,5,6],function_pairs_RNA2,function_pairs_RNA2,colors,file_label="ncRNAs")
    #plot_corr_k_categ(domain,orgas[domain],organism_names,[1,2,3,4,5,6],function_pairs_miRNA,function_pairs_miRNA,colors,file_label="miRNA")

    #plot_coverage_chromosomal(organism_names,function_names,domain)
    #plot_corr_all_orgas(domain_dict,orga_names_all,[2,5],function_name_lim,sorted=True)
    #print(orgas[domain])
    #plot_func_heatmap(orgas[domain],domain,func_name1,func_name2,k_value=3,sorted=True,func_label="",figsize=(10,10))

    #plot_coverage_heatmap(organism_names,function_name_lim1,domain,log = True,square=False,sorted=False, label= "standart")
    #plot_coverage_heatmap(organism_names,function_name_lim1,domain,log = True,square=False,sorted=False, label= "standart")
    #plot_dimer_compare_model("Homo sapiens","animalia","genes")
    #plot_dimer_compare_elems("Homo sapiens","animalia",function_name_cds)
    #plot_dimer_compare_elems(organism_names[0],domain,function_name_cds)
    #plot_prop_heatmap(organism_names,function_name_lim1,domain,"length",'length_hm',log = True,square=False,sorted=True)
    #plot_GC_bar(orgas[domain],func_name="genes",group_name=domain,sort_orga=True,sort_GC=False)
    #plot_word_bar_deviation_mean_categ(orgas[domain],function_name_RNA,function_name_RNA,domain,count_func=dimer_minus_expectation_norm,colors=colors,)
    #plot_word_bar_deviation_mean_categ(orgas[domain],function_name_RNA,function_name_RNA,domain,count_func=dimer_minus_expectation_norm,colors=colors,folder="deviation_mean_RNA")
    #plot_dimerBox_compare(organisms,domain,'exons','introns')
    #for func_name in function_name_lim:
        #plot_dimer_hm(organisms,domain=domain,domain_dict=None,func_name=func_name)
    
    #plot_dimer_mean_devs(organism_names,function_name_RNA2,function_name_RNA2,domain=domain,label=domain)

    func_names = ['exons','introns','intergenics','lncRNAs','tRNAs','rRNAs']+function_name_RNA
    func_labels = ['exons','introns','intergenics','lncRNAs','tRNAs','rRNAs']+function_name_RNA
    #func_names = ['pc_exons','pseudo_exons','npc_exons','pc_introns','pseudo_introns','npc_introns',]
    #func_labels = ['pc-exons','pseudo-exons','npc-exons','pc-introns','pseudo-introns','npc-introns',]

    for func_name,func_label in zip(func_names[0:3],func_labels[0:3]):
        #plot_dimer_word_count(organism_names,func_name,func_label,domain=domain,filter_nan=False,dimer_bias=True)
        pass


    #plot_dimer_sum_diffs(organism_names,function_pairs,domain=domain,label=None)
    #plot_dimer_cum_bias(organism_names,function_names_hm,function_names_hm,domain=domain,label=None,filter_nan=False,total_exp=True)
    #plot_dimer_cum_bias_diff(organism_names,function_pairs_introns,function_pairs_introns,domain=domain,label="pearson")
    for orga_name in organism_names[0:1]:
        #plot_dimer_orga_categ_bias(orga_name,function_name_cds,function_name_cds,domain=None,label=None)
        #plot_dimer_orga_categ_bias(orga_name,function_name_cds,function_name_cds,domain=None,label=None)
        #plot_dimer_compare_model(orga_name,domain,function_name_cds[1])
        pass
    #plot_dimer_categ_bias_box(organism_names,["intergenics","exons"],["intergenics","exons"],domain=None,label=None)

    for func_name in function_name_lim:
        #plot_func_heatmap(organism_names,func_name,func_name,k_value=5,sorted=False,filter_empty=True,figsize=(9,9),domain=domain)
        pass

    #plot_dimer_covariance_heatmap(organism_names,["intergenics","introns","npc_exons","pseudo_exons","pc_exons","cds"],domain=domain)
    #plot_dimer_covariance_heatmap_combined(organism_names,covariance_funcs,organisms_label=domain,label="multiple genome regions",domain=domain)
    #plot_bar_length(organism_names, domain_list=[domain for name in organism_names],title="genome length of organisms in "+domain)


orga_names_sorted = kmer_tools.download_genbank.sort_orga_names(orga_names_all,dicto=True)
domain_list=[kmer_tools.download_genbank.name_to_domain(orga_name) for orga_name in orga_names_sorted]
for Word in dimer_dist.dimer_words:
    #plot_dimer_word_bias_orgas(orga_names_sorted,"exons","exons",Word,label=None,domain_list=domain_list)
    pass
#plot_func_heatmap(orga_names_all,domain_dict,"exons","intergenics",k_value=5,sorted=True,filter_empty=True,figsize=(23,23))
print([(name,domain) for name,domain in zip(orga_names_sorted,domain_list)])
#plot_func_heatmap(orga_names_sorted,"exons","exons",domain_list=domain_list,k_value=5,sorted=False,filter_empty=True,figsize=(20,20),split_class=False)


orga_names_select=[]
for domain in [
    "animalia",
    "embryophyta",
    "fungi",
    "protista",
    #"mitochondria",
    #"plastids"
    ]:
    orgas[domain] = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga in orgas[domain]:
        domain_dict[orga.name]=domain
    orga_names_select += [orga.name for orga in orgas[domain]]

funcs_normal = [
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
    #'gaps',
    "pc_exons",
    "pc_introns",
    "npc_exons",
    "npc_introns",
    "pseudo_exons",
    "pseudo_introns",
]

func_rna = [
    'non_protein_coding_genes',
    'ncRNAs',
    'lncRNAs',
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'tRNAs',
    'rRNAs',
    'precursor_RNAs',
]

function_names_solo=[
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'rRNAs',
    'precursor_RNAs',
    'gaps',
]


#plot_func_heatmap(orga_names_all,"introns","intergenics",domain_dict=domain_dict, v_min=-1,v_max=1,k_value=3,sorted=True,filter_empty=True,figsize=(22,22))
def plot_step(i):
    plot_func_heatmap(orga_names_all,function_names_solo[i],function_names_solo[i],domain_dict=domain_dict, v_min=-1,v_max=1,k_value=3,sorted=True,filter_empty=True,figsize=(22,22))
    plot_func_heatmap(orga_names_all,function_names_solo[i],function_names_solo[i],domain_dict=domain_dict, v_min=-1,v_max=1,k_value=5,sorted=True,filter_empty=True,figsize=(22,22))
    pass
func_nam = function_names+merge_func_names
n=len(function_names_solo)
#Parallel(n_jobs=n)(delayed(plot_step)(i) for i in range(n))



#plot_func_heatmap(orga_names_all,"exons","exons",domain_dict=domain_dict,k_value=5,sorted=True,filter_empty=True,figsize=(22,22))

def plot_step(i):
    func1=func_nam[i]
    for func2 in func_nam[0:i+1]:
        for k_value in [
            3,
            5,
            ]:
            plot_func_heatmap(orga_names_all,func1,func2,domain_dict=domain_dict, v_min=-1,v_max=1,k_value=k_value,sorted=True,filter_empty=True,figsize=(22,22))
            pass
n=len(func_nam)
#Parallel(n_jobs=n)(delayed(plot_step)(i) for i in range(n))


def plot_step_mito(i):
    func=func_nam[i]
    for k_value in [
        #3,
        5,
        ]:
        plot_func_heatmap(Mitochondria_Names_sorted_All_NCBI,func,func,domain_list=["mitochondria" for name in Mitochondria_Names_sorted_All_NCBI],k_value=k_value,sorted=False,filter_empty=True,figsize=(15,15),split_class=True)
        pass
n=len(func_nam)

#plot_func_heatmap(Mitochondria_Names_sorted_All_NCBI,"exons","exons",domain_list=["mitochondria" for name in Mitochondria_Names_sorted_All_NCBI],k_value=5,sorted=False,filter_empty=True,figsize=(15,15),split_class=True)
#Parallel(n_jobs=n)(delayed(plot_step_mito)(i) for i in range(n))


# ------------------------------------------------------------------------------------
# ----------------------------- ORGANISMS GROUPS    ----------------------------------
# ------------------------------------------------------------------------------------

function_pairs = [
  ['intergenics','exons'],
  ['intergenics','introns'],
]
function_pairs_pc = [
  ['intergenics','pc_exons'],
  ['intergenics','pc_introns'],
]
function_pairs_npc = [
  ['intergenics','npc_exons'],
  ['intergenics','npc_introns'],
]
function_pairs_pseudo = [
  ['intergenics','pseudo_exons'],
  ['intergenics','pseudo_introns'],
]

function_pairs_RNA=[
    #['ncRNAs','mRNAs'],
    ['ncRNAs','lncRNAs'],
    ['ncRNAs','miRNAs'],
    ['ncRNAs','snoRNAs'],
    ['ncRNAs','snRNAs'],
    ['ncRNAs','rRNAs'],
    ['ncRNAs','tRNAs'],
]

function_pairs_RNA_total=[
    #['ncRNAs','mRNAs'],
    ['total_genome','lncRNAs'],
    ['total_genome','miRNAs'],
    ['total_genome','snoRNAs'],
    ['total_genome','snRNAs'],
    ['total_genome','rRNAs'],
    ['total_genome','tRNAs'],
]
function_pairs_RNA_inter=[
    #['ncRNAs','mRNAs'],
    ['intergenics','lncRNAs'],
    ['intergenics','miRNAs'],
    ['intergenics','snoRNAs'],
    ['intergenics','snRNAs'],
    ['intergenics','rRNAs'],
    ['intergenics','tRNAs'],
]

function_pairs_RNA_pc=[
    #['ncRNAs','mRNAs'],
    ['pc_exons','lncRNAs'],
    ['pc_exons','miRNAs'],
    ['pc_exons','snoRNAs'],
    ['pc_exons','snRNAs'],
    ['pc_exons','rRNAs'],
    ['pc_exons','tRNAs'],
]

function_pair_titles = [
  ['pc','pseudo'],  
  ['intergenics','pc'],
  ['intergenics','pseudo'],
  ['npc','pc'],
  ['npc','pseudo'],
  ['npc','intergenics'],
]

function_pairs_pc_comp=[
    ['pc_introns','npc_introns'],
    ['pc_exons','npc_exons'],
    ['pc_introns','pseudo_introns'],
    ['pc_exons','pseudo_exons'],
]

function_pairs_npc_comp=[
    ['npc_introns','pseudo_introns'],
    ['npc_exons','pseudo_exons'],
]

function_pairs_ncRNAs=[
    ['ncRNAs','lncRNAs'],
    ['ncRNAs','miRNAs'],
    ['ncRNAs','snoRNAs'],
    ['ncRNAs','snRNAs'],
    #['ncRNAs','rRNAs'],
    #['ncRNAs','tRNAs'],
]

function_pairs_trRNAs=[
    ['ncRNAs','rRNAs'],
    ['ncRNAs','tRNAs'],
    ['tRNAs','rRNAs'],
]
function_pairs_trRNAs=[
    ['total_genome','rRNAs'],
    ['total_genome','tRNAs'],
    ['tRNAs','rRNAs'],
]

orga_names_select = []
domain_list = []
group_list = []

vertebrate_names = [name for name in Animalia_Names_NCBI if 
            Animalia_organism_class[name] == "Mammalia" or 
            Animalia_organism_class[name] == "Aves" or 
            Animalia_organism_class[name] == "Rep., Amphib. and Fish" or
            Animalia_organism_class[name] == "Rep., Am., Fish"]

invertebrate_names = [name for name in Animalia_Names_NCBI if 
            Animalia_organism_class[name] == "Invertebrates" or 
            Animalia_organism_class[name] == "Insecta"]


orga_names_select += vertebrate_names
orga_names_select += invertebrate_names
domain_list += ["animalia" for orga in vertebrate_names]
group_list  += ["vertebrates" for orga in vertebrate_names]
domain_list += ["animalia" for orga in invertebrate_names]
group_list  += ["invertebrates" for orga in invertebrate_names]
#print([x for x in zip(orga_names_select,domain_list,group_list)])

eudicot_names = Embryophyta_Names_NCBI[1:10]
#print (eudicot_names)
monocotyledon_names = Embryophyta_Names_NCBI[11:]
#print (monocotyledon_names)
protist_pure_names = Protista_Names_NCBI[5:-6]
#print (protist_pure_names)

for domain in [
    #"animalia",
    "embryophyta",
    "protista",
    "fungi",
    #"mitochondria",
    #"plastids"
    ]:
    orgas[domain] = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga in orgas[domain]:
        domain_dict[orga.name]=domain
    orga_names_select += [orga.name for orga in orgas[domain]]
    domain_list += [domain for orga in orgas[domain]]
    group_list += [domain for orga in orgas[domain]]

group_names = ["vertebrates", "invertebrates", "embryophyta", "protista", "fungi"]
#plot_corr_OrgaGroup(orga_names_select, domain_list, group_names, group_list, 5, function_pairs_RNA[0:4], function_pairs_RNA[0:4], colors, file_label="ncRNA correlations")
#plot_corr_OrgaGroup(orga_names_select, domain_list, group_names, group_list, 5, function_pairs_RNA_total[0:4], function_pairs_RNA_total[0:4], colors, file_label="ncRNA - total genome correlations")
#plot_corr_OrgaGroup(orga_names_select, domain_list, group_names, group_list, 5, function_pairs_RNA_inter, function_pairs_RNA_inter, colors, file_label="ncRNA - intergenics correlations")
#plot_corr_OrgaGroup(orga_names_select, domain_list, group_names, group_list, 5, function_pairs_RNA_pc, function_pairs_RNA_pc, colors, file_label="ncRNA - pc exons correlations")

function_names_bias=[
    "exons",
    "introns",
    "cds",
    "intergenics",
    #"tRNAs",
    #"rRNAs"
    ]

for func_pairs,label in zip([
    function_pairs_pc,
    #function_pairs_npc,
    function_pairs_pseudo
    ],[
        "pc genes",
        #"npc genes",
        "pseudo genes"
        ]):
    label = "intergenics - exons and introns in " + label
    #plot_corr_OrgaGroup(orga_names_select, domain_list, group_names, group_list, 5, func_pairs, func_pairs, colors, file_label=label)
    #plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [2,3,5], func_pairs, func_pairs, colors, file_label=label)
label = "pc genes - npc genes correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [2,3,5], function_pairs_pc_comp[0:2], function_pairs_pc_comp[0:2], colors, file_label=label)
label = "pc genes - pseudo genes correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [2,3,5], function_pairs_pc_comp[2:4], function_pairs_pc_comp[2:4], colors, file_label=label)
label = "npc genes - pseudo genes correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [2,3,5], function_pairs_npc_comp, function_pairs_npc_comp, colors, file_label=label)
label = "ncRNAs correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [5], function_pairs_ncRNAs, function_pairs_ncRNAs, colors, file_label=label)
label = "tRNAs rRNAs total correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [5], function_pairs_trRNAs, function_pairs_trRNAs, colors, file_label=label)
label = "ncRNAs internal correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [5], function_pairs_RNA2, function_pairs_RNA2, colors, file_label=label)
label = "miRNAs internal correlations"
#plot_corr_OrgaGroup_k(orga_names_select, domain_list, group_names, group_list, [5], function_pairs_miRNA, function_pairs_miRNA, colors, file_label=label)

#plot_dimer_covariance_heatmap(vertebrate_names,["intergenics","introns","npc_exons","pseudo_exons","pc_exons","cds"],organisms_label="vertebreates",domain="animalia")
#plot_dimer_covariance_heatmap(invertebrate_names,["intergenics","introns","npc_exons","pseudo_exons","pc_exons","cds"],organisms_label="invertebreates",domain="animalia")
#plot_dimer_covariance_heatmap(vertebrate_names,["intergenics","introns","npc_exons","pseudo_exons","pc_exons","cds"],organisms_label="vertebrates",domain="animalia")
#plot_dimer_covariance_heatmap(invertebrate_names,["intergenics","introns","npc_exons","pseudo_exons","pc_exons","cds"],organisms_label="invertebrates",domain="animalia")
covariance_funcs = [
    "intergenics",
    "pc_introns",
    "npc_introns",
    "pseudo_introns",
    "npc_exons",
    "pseudo_exons",
    "pc_exons",
    "cds"
    ]

#plot_dimer_covariance_heatmap_combined(invertebrate_names,covariance_funcs,organisms_label="invertebrates",label="multiple genome regions",domain="animalia")
#plot_dimer_covariance_heatmap_combined(vertebrate_names,covariance_funcs,organisms_label="vertebrates",label="multiple genome regions",domain="animalia")
#plot_dimer_covariance_heatmap_combined(orga_names_select,covariance_funcs,organisms_label="all organisms",label="multiple genome regions",domain_list=domain_list)
for func_name in covariance_funcs:
    #plot_dimer_covariance_heatmap_combined(orga_names_select,[func_name],organisms_label="all organisms",label=func_name,domain_list=domain_list)
    pass

for func_name,func_label in zip (["cds", "exons", "introns", "intergenics"], ["coding sequences", "exons", "introns", "intergenics"]):
    domains =  [
        "animalia", 
        "embryophyta", 
        "protista", 
        "fungi",
        ]
    #plot_dimer_TpA_CpG_scatter(orga_names_select, domains, domain_list, group_names, func_name, func_label, colors=colors, file_label="")
    #plot_dimer_CpA_CpG_scatter(orga_names_select, domains, domain_list, group_names, func_name, func_label, colors=colors, file_label="")
    #plot_dimer_CpG_vs_GC_scatter(orga_names_select, domains, domain_list, group_names, func_name, func_label, colors=colors, file_label="",TpA=True)

funcs = ["pc_exons","pseudo_exons","npc_exons"]
for domain in domains:
    #plot_dimer_TpA_CpG_scatter_func(orga_names_select, "animalia", domain_list, group_names, funcs, funcs, colors=colors, file_label="")
    #plot_dimer_word_scatter_func("GA","AT", orga_names_select,domain, domain_list, group_names, funcs, funcs, colors=colors, file_label="")
    #plot_dimer_word_scatter_func_chromo("CA","CG", orga_names_select,domain, domain_list, group_names, funcs, funcs, colors=colors, file_label="")
    #plot_dimer_word_scatter_func_chromo("TA","CG", orga_names_select,domain, domain_list, group_names, funcs, funcs, colors=colors, file_label="")
    #plot_dimer_word_scatter_func_chromo("GA","AT", orga_names_select,domain, domain_list, group_names, funcs, funcs, colors=colors, file_label="")
    pass


func_sequence=function_name_cds
for classification in ["Monocots","Eudicots","Euglenozoa","Protists","Algae","higher Fungi","Mammalia","Aves","Rep., Am., Fish","Invertebrates","Insecta"]:
    #plot_dimer_mean_devs([name for name in Sorted_Names_All_NCBI if all_orga_classes[name] == classification  ],func_sequence,label=classification)
    pass

for domain in domains[0:4]:
    organism_names = [orga.name for orga in orgas[domain]]
    organism_names = kmer_tools.download_genbank.sort_orga_names(organism_names)
    #plot_dimer_mean_devs(organism_names,func_sequence,domain=domain,label=domain)

#plot_dimer_mean_devs([name for name in Animalia_Names_NCBI if (Animalia_organism_class[name] == "Mammalia" or Animalia_organism_class[name] =="Aves") ],func_sequence,label="mammalia and aves")
#plot_dimer_mean_devs(vertebrate_names,func_sequence,domain=None,label="vertebrates")
#plot_dimer_mean_devs(invertebrate_names,func_sequence,domain=None,label=" all invertebrates")



#plot_dimer_mean_devs([name for name in Animalia_Names_NCBI if Animalia_organism_class[name] == "Mammalia"],function_name_RNA2,domain="animalia",label="mammalia RNAs full")
#plot_dimer_mean_devs(invertebrate_names,function_name_RNA2,domain="animalia",label="invertebrate RNAs full")
#plot_dimer_mean_devs(vertebrate_names,function_name_RNA2,domain="animalia",label="vertebrate RNAs full")
#plot_dimer_mean_devs(monocotyledon_names,function_name_RNA2,domain="embryophyta",label="monocotyledon RNAs full")
#plot_dimer_mean_devs(eudicot_names,function_name_RNA2,domain="embryophyta",label="eudicots RNAs full")
#plot_dimer_mean_devs(eudicot_names,function_name_RNA2,domain="embryophyta",label="eudicots RNAs full")
#plot_dimer_mean_devs([name for name in Protista_Names_NCBI if Protista_organism_classv2[name] == "Algen" or Protista_organism_classv2[name] == "Euglenozoa"],function_name_RNA2,domain="protista",label="algae and euglenozoa RNAs full")



#plot_dimer_cum_bia_OrgaGroups(orga_names_select, domain_list, group_names, group_list, 
#    function_names=function_names_bias, 
#    function_titles=function_names_bias
#    )

#------------------------------- Chlorop., Algae, Embryophyta ---------------------

for domain in [
    #"animalia",
    "embryophyta",
    #"fungi",
    "protista",
    #"mitochondria",
    "plastids",
    ]:
    orgas[domain] = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga in orgas[domain]:
        domain_dict[orga.name]=domain
    orga_names_select += [orga.name for orga in orgas[domain]]

orga_names_select = []
domain_list = []
label_list = []

orga_names_select += sort_nucleo_names(Nucleomorph_names_sorted_NCBI)
domain_list += ["nucleomorph" for i in range(len(Nucleomorph_names_sorted_NCBI))]
label_list  += ["nucleo" for i in range(len(Nucleomorph_names_sorted_NCBI))]

orga_names_select += sort_plas_names(Plastid_Names_sorted_All_NCBI[4:])
domain_list += ["plastids" for i in range(len(Plastid_Names_sorted_All_NCBI[4:]))]
label_list  += ["chloro" for i in range(len(Plastid_Names_sorted_All_NCBI[4:]))]

algae_names = [name for name in Protista_Names_NCBI if Protista_organism_class[name] == "Algae"]
orga_names_select += sort_orga_names(algae_names)
domain_list += ["protista" for i in range(len(algae_names))]
label_list  += ["algae" for i in range(len(algae_names))]

orga_names_select += sort_orga_names(Embryophyta_Names_NCBI)
domain_list += ["embryophyta" for i in range(len(Embryophyta_Names_NCBI))]
label_list  += ["embryophyta" for i in range(len(Embryophyta_Names_NCBI))]

#print(orga_names_select)
#print(domain_list)

func_nam = ["exons","introns"]
#(plot_func_heatmap)(orga_names_select,func_nam[0],func_nam[0],domain_list=domain_list,v_min=-1,v_max=1,k_value=5,sorted=False,filter_empty=True,figsize=(15,15)) 
for func in func_nam:
    #(plot_func_heatmap)(orga_names_select,func,func,domain_list=domain_list,v_min=-1,v_max=1,k_value=5,sorted=False,filter_empty=True,figsize=(13,13)) 
    pass

#Parallel(n_jobs=len(func_nam))(
#    delayed(plot_func_heatmap)(orga_names_select,func,func,domain_list=domain_list,k_value=5,sorted=False,filter_empty=True,figsize=(15,15),split_class=True) 
#    for func in func_nam
#    )

# ---------------------------------------------------------------------------------
# ----------------------------------- REPEATS    ----------------------------------
# ---------------------------------------------------------------------------------

repeat_funcs_labels = kmer_tools.search_kmer.repeat_funcs_labels
repeat_funcs_labels_SINE = kmer_tools.search_kmer.repeat_funcs_labels[1:13]
repeat_funcs_labels_LINE = kmer_tools.search_kmer.repeat_funcs_labels[13:]

repeat_funcs_labels_SINE2 = [
    "SINE",
    #"SINE_Alu",
    "SINE_Alu_pure",
    #"SINE_Alu_FAM",
    #"SINE_Alu_FRAM",
    #"SINE_Alu_FLAM",
    #"SINE_Alu_PB1",
    #"SINE_Alu_B1",
    "SINE_B2",
    "SINE_B4",
    #"SINE_ID",
    "CSINE",
]
repeat_func_names_SINE2=[
    "SINE (all)",
    #"SINE_Alu",
    "SINE Alu",
    #"SINE_Alu_FAM",
    #"SINE_Alu_FRAM",
    #"SINE_Alu_FLAM",
    #"SINE_Alu_PB1",
    #"SINE_Alu_B1",
    "SINE B2",
    "SINE B4",
    #"SINE_ID",
    "C SINE",
]
#print(repeat_funcs_labels_SINE)
#print(repeat_funcs_labels_LINE)

func_name1 = "SINE_Alu_pure"
func_name2 = "SINE_Alu_B1"

dict_repeats = {}
orga_names_all = []
for dict,domain in [
    (animalia_ids,"animalia"),
    (embryophyta_ids,"embryophyta"),
    (fungi_ids,"fungi")
    ]:  
    dict_repeats.update(dict)
    orga_names = [name for name in dict] 
    orga_names_all += orga_names
    for name in orga_names:
        domain_dict[name]=domain
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    domain_names = [orga.name for orga in organisms]

    #plot_coverage_heatmap(orga_names,repeat_funcs_labels_SINE,domain,log = True,square=False,sorted=True, label= "repeat_SINE")
    #plot_coverage_heatmap(orga_names,repeat_funcs_labels_LINE,domain,log = True,square=False,sorted=True, label= "repeat_LINE")
    #plot_func_heatmap(organisms,domain_dict,func_name1,func_name2,k_value=5,sorted=True,filter_empty=False)

orga_names = [name for name in dict_repeats] 
orga_names_all=kmer_tools.download_genbank.sort_orga_names(orga_names)
animalia_names = [name for name in animalia_ids]
domain_names = [orga.name for orga in organisms]
domain_list=[kmer_tools.download_genbank.name_to_domain(name) for name  in orga_names_all]

#plot_coverage_heatmap(orga_names,repeat_funcs_labels,domain,log = True,square=False,sorted=True, label= "repeat_full")
#plot_func_heatmap(orga_names,domain_dict,"SINE","LINE",k_value=5,sorted=True,filter_empty=True,figsize=(13,13))
plot_coverage_heatmap_list(orga_names_all,repeat_funcs_labels_SINE2,repeat_func_names_SINE2,domain_list,log = True,square=False,sorted=True, label= "repeat_SINE")
#plot_func_heatmap(orga_names_all,"SINE","SINE",domain_list=domain_list,k_value=5,sorted=False,filter_empty=True,figsize=(12,12),split_class=True)
funcs_repeat = [
    "SINE",
    "LINE",
    "LTR",
    ]

funcs_LINE = [
    "LINE",
    "LINE_L1",
    "LINE_L2",
    "LINE_L3",
]

for i,func1 in enumerate(funcs_repeat):
    for func2 in funcs_repeat[0:i+1]:
        for k in [
            3,
            5
            ]:
            #plot_func_heatmap(orga_names,domain_dict,func1,func2,v_min=-1,v_max=1,k_value=k,sorted=True,filter_empty=True,figsize=(12,12))
            pass

def plot_step(i):
    func1=funcs_repeat[i]
    for func2 in funcs_repeat[0:i+1]:
        for k in [
            3,
            5,
            ]:
            #plot_func_heatmap(orga_names,func1,func2,k_value=k,domain_dict=domain_dict,v_min=-1,v_max=1,sorted=True,filter_empty=True,figsize=(13,13),split_class=True)
            pass

#plot_funcs_mul_heatmap(orga_names,funcs_repeat,funcs_repeat,k_value=5,domain_dict=domain_dict,v_min=-1,v_max=1,sorted=True,filter_empty=True,figsize=(20,20),split_class=False,mark_wrap=True)
n=len(funcs_repeat)
#Parallel(n_jobs=n)(delayed(plot_step)(i) for i in range(n))

#plot_corr_all_orgas(domian_dict,animalia_repeat_names,[2,3],function_names_repeat,repeatmask=True)
#plot_corr_chromosomal(domain,repeat_names,[2],repeat_functions)
#plot_corr_k(domain,organism_names,[1,2,3],repeat_functions,colors)
#plot_coverage_chromosomal(repeat_names,repeat_functions,domain)
