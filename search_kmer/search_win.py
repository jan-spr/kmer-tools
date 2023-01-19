import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
sys.path.append('../..')                                    #add parent folder

import Oligo
import kmer_tools

from kmer_search_funcs import *
from reorganize_funcs import *
from repeat_search_funcs import *
from loci_save_lib import *

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
    'gaps'
]

functions_merge = kmer_tools.search_kmer.merge_funcs
function_names_merge = kmer_tools.search_kmer.merge_func_names
functions_merge_dict={}
for func,name in zip(functions_merge,function_names_merge):
    functions_merge_dict[name] = func

functions_merge_loci = kmer_tools.search_kmer.merge_funcs_loci
function_names_merge = kmer_tools.search_kmer.merge_func_names
functions_merge_dict_loci={}
for func,name in zip(functions_merge,function_names_merge):
    functions_merge_dict[name] = func

functions_merge_Oligo = kmer_tools.search_kmer.merge_funcs_Oligo
functions_merge_dict_Oligo={}
for func,name in zip(functions_merge_Oligo,function_names_merge):
    functions_merge_dict_Oligo[name] = func

search_kmer_funcs= functions + functions_merge_Oligo
search_kmer_funcs_names= function_names + function_names_merge

functions_merge_single = [
    Oligo.File.read_exons,
    Oligo.File.read_introns,
    Oligo.File.read_protein_coding_genes,
    Oligo.File.read_pseudo_genes,
    Oligo.File.read_non_protein_coding_genes,
]

function_names_merge_single = [
    'exons',
    'introns',
    'protein_coding_genes',
    'pseudo_genes',
    'non_protein_coding_genes',
]

functions_merge_single_dict={}
for func,name in zip(functions_merge_single,function_names_merge_single):
    functions_merge_single_dict[name] = func


repeat_function_names=[
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

domain = "embryophyta"
domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids","nucleomorph"]

#organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
#for orga_group in [Fungi_names_Basidiomycota_p,Fungi_names_Ascomycota]:
for domain in domains:
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    domain_names = [name for name in fungi_ids]
    orga_names = [orga.name for orga in organisms]
    #for names in [domain_names]:       
    #    func_loop_kmerOpt(names,domain=domain,functions = functions,function_names = function_names,k_values=[1,2,3],Windows=False,stop_overwrite=True)
    print(domain)
    store_kmer_chromosomal(organisms,search_kmer_funcs_names+["total_genome"],domain,k_values=[1,2])
    #store_kmer_orga(organisms,function_names=search_kmer_funcs_names+["total_genome"],domain=domain,k_values=[1,2])
    fix_kmer_filename_windows("../../data/k_spectra_strand/"+domain+"/")
    #func_loop(orga_names,domain,save_loci,function_names_merge_single,function_dict=functions_merge_single_dict,Windows=True,stop_overwrite=False)


for dict,domain in [
    (animalia_ids,"animalia"),
    (embryophyta_ids,"embryophyta"),
    (fungi_ids,"fungi")
    ]:
    print(domain)      
    orga_names = [name for name in dict] 
    organisms = kmer_tools.download_genbank.names_to_orgas(orga_names)
    #func_loop(orga_names[20:60], domain, save_loci, function_names=[repeat_funcs_labels[i] for i in repeat_v2_indeces], functions=[repeat_funcs_save[i] for i in repeat_v2_indeces], Windows=True, stop_overwrite=False)
    store_kmer_chromosomal(organisms,repeat_funcs_labels,domain,k_values=[1,2])
    #store_kmer_orga(organisms,function_names=repeat_funcs_labels,domain=domain,k_values=[1,2,3,4,5,6])


print([repeat_funcs_labels[i] for i in [0,13,14,15,16]])
print([repeat_funcs_labels[i] for i in [1,2,3,4,5,6,7,8,9,10,11,12]])