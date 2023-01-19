import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                     #Windows Oligo
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')           #Linux Oligo
sys.path.append('../..')                                    #add parent folder

import Oligo
import kmer_tools
from kmer_tools.search_kmer.repeat_search_funcs import *
from elem_search_funcs import *

functions=[
    Oligo.File.read_genes,
    Oligo.File.read_exons,
    Oligo.File.read_introns,
    Oligo.File.read_cds,
    Oligo.File.read_mRNAs,
    Oligo.File.read_intergenics,
    Oligo.File.read_protein_coding_genes,
    Oligo.File.read_pseudo_genes,
    #Oligo.File.read_pseudo_genes_exons,
    #Oligo.File.read_pseudo_genes_introns,
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


functions_repeat=[
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
    Oligo.File.read_gaps,
    #read_LTR,
    #read_SINE,
    #read_SINE_Alu,
    #read_LINE,
    #read_LINE_L1,
    #read_LINE_L2,
]

function_names_repeat=[
    'genes',
    'exons',
    'introns',
    'cds',
    'mRNAs',
    'intergenics',
    'protein_coding_genes',
    'pseudo_genes',
    'pseudo_genes_exons',
    'pseudo_genes_introns',
    'non_protein_coding_genes',
    'ncRNAs',
    'lncRNAs',
    'miRNAs',
    'snoRNAs',
    'snRNAs',
    'tRNAs',
    'rRNAs',
    'precursor_RNAs',
    'gaps',
    #'LTR',
    #'SINE',
    #'SINE_Alu',
    #'LINE',
    #'LINE_L1',
    #'LINE_L2',
]

functions_merge_Oligo = kmer_tools.search_kmer.merge_funcs_Oligo
function_names_merge = kmer_tools.search_kmer.merge_func_names

search_kmer_funcs= functions + functions_merge_Oligo
search_kmer_funcs_names= function_names + function_names_merge


functions_merge_loci = kmer_tools.search_kmer.merge_funcs_loci
function_names_merge = kmer_tools.search_kmer.merge_func_names
"""
functions_repeat_dict={}
for func,name in zip(functions_repeat,function_names_repeat):
    functions_repeat_dict[name] = func

function_names = [    
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
    ]

functions = [functions_repeat_dict[name] for name in function_names]
"""


domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids","nucleomorph"]
for domain in ["nucleomorph"]:
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    #func_loop_elem(organisms,domain,functions_repeat,function_names_repeat)
    names = [orga.name for orga in organisms]
    func_loop_elem(names, domain, functions=search_kmer_funcs, function_names=search_kmer_funcs_names)

for dict,domain in [
    (animalia_ids,"animalia"),
    (embryophyta_ids,"embryophyta"),
    (fungi_ids,"fungi"),
    ]:
    orga_names = [name for name in dict] 
    #func_loop_elem(orga_names,domain,functions,function_names)
    #func_loop_elem(orga_names, domain, functions=functions_merge_loci, function_names=function_names_merge)
