#functions for reading certain gene-elements from chromosome
#mainly for repeatmasker data as gb-data is read with Oligo-functions

import sys
from tabnanny import verbose

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('../..')   
import kmer_tools

def merge_loci_func(func1,func2,chromosome):
    print(func1)
    loci1=func1(chromosome=chromosome)
    print(func2)
    loci2=func2(chromosome=chromosome)

    if (len(loci1) == 0 or len(loci2) == 0):
        return []
    else:
        return Oligo.Loci.substract(loci1,Oligo.Loci.invert(loci2,verbose=0),verbose=0)

def read_protein_coding_genes_loci(chromosome):
    domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    func_name = "protein_coding_genes"
    loci = kmer_tools.search_kmer.read_loci(chromosome,domain,func_name)


def read_protein_coding_genes_exons(chromosome):
    return merge_loci_func(
        Oligo.File.read_protein_coding_genes,
        Oligo.File.read_exons,
        chromosome)

def read_protein_coding_genes_introns(chromosome):
    return merge_loci_func(
        Oligo.File.read_protein_coding_genes,
        Oligo.File.read_introns,
        chromosome)

def read_non_protein_coding_genes_exons(chromosome):
    return merge_loci_func(
        Oligo.File.read_non_protein_coding_genes,
        Oligo.File.read_exons,
        chromosome)

def read_non_protein_coding_genes_introns(chromosome):
    return merge_loci_func(
        Oligo.File.read_non_protein_coding_genes,
        Oligo.File.read_introns,
        chromosome)

def read_pseudo_genes_exons(chromosome):
    return merge_loci_func(
        Oligo.File.read_pseudo_genes,
        Oligo.File.read_exons,
        chromosome)

def read_pseudo_genes_introns(chromosome):
    return merge_loci_func(
        Oligo.File.read_pseudo_genes,
        Oligo.File.read_introns,
        chromosome)

def merge_loci_func_l(func_name1,func_name2,chromosome):
    domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    loci1 = kmer_tools.search_kmer.read_loci(chromosome,domain,func_name1)
    loci2 = kmer_tools.search_kmer.read_loci(chromosome,domain,func_name2)
    if (len(loci1) == 0 or len(loci2) == 0):
        return []
    else:
        return Oligo.Loci.substract(loci1,Oligo.Loci.invert(loci2,verbose=0),verbose=0)

def read_protein_coding_genes_exons_l(chromosome):
    return merge_loci_func_l(
        "protein_coding_genes",
        "exons",
        chromosome)

def read_protein_coding_genes_introns_l(chromosome):
    return merge_loci_func_l(
        "protein_coding_genes",
        "introns",
        chromosome)

def read_non_protein_coding_genes_exons_l(chromosome):
    return merge_loci_func_l(
        "non_protein_coding_genes",
        "exons",
        chromosome)

def read_non_protein_coding_genes_introns_l(chromosome):
    return merge_loci_func_l(
        "non_protein_coding_genes",
        "introns",
        chromosome)

def read_pseudo_genes_exons_l(chromosome):
    return merge_loci_func_l(
        "pseudo_genes",
        "exons",
        chromosome)

def read_pseudo_genes_introns_l(chromosome):
    return merge_loci_func_l(
        "pseudo_genes",
        "introns",
        chromosome)

"""
def merge_loci_func(func1,func2,chromosome):
    print(func1)
    loci1=func1(chromosome=chromosome)
    print(func2)
    loci2=func2(chromosome=chromosome)

    if (len(loci1) == 0 and len(loci2) == 0):
        print("len0",loci1,loci2)
        return []
    else:
        print("loci_loci")
        return Oligo.Loci.loci_in_loci(loci1,loci2)
"""

"""
def def_get_func_exons(func,chromosome,verbose=0):
    loci=func(chromosome=chromosome, add_ex_in=True)
    exons = []
    for gene in loci:
        exons += gene.exons
    if chromosome is not None:
        Oligo.Loci.add_target(exons, chromosome)
    if verbose:
        if not exons:
            print('No valid exons found in '+str(chromosome))
    return exons

def def_get_func_introns(func,chromosome,verbose=0):
    loci=func(chromosome=chromosome)
    introns = []
    for gene in loci:
        introns += gene.introns
    if chromosome is not None:
        Oligo.Loci.add_target(introns, chromosome)
    if verbose:
        if not introns:
            print('No valid exons found in '+str(chromosome))
    return introns

def read_protein_coding_genes_exons(chromosome):
    return def_get_func_exons(
        Oligo.File.read_protein_coding_genes,
        chromosome)

def read_protein_coding_genes_introns(chromosome):
    return def_get_func_introns(
        Oligo.File.read_protein_coding_genes,
        chromosome)

def read_non_protein_coding_genes_exons(chromosome):
    return def_get_func_exons(
        Oligo.File.read_non_protein_coding_genes,
        chromosome)

def read_non_protein_coding_genes_introns(chromosome):
    return def_get_func_introns(
        Oligo.File.read_non_protein_coding_genes,
        chromosome)

def read_pseudo_genes_exons(chromosome):
    return def_get_func_exons(
        Oligo.File.read_pseudo_genes,
        chromosome)

def read_pseudo_genes_introns(chromosome):
    return def_get_func_introns(
        Oligo.File.read_pseudo_genes,
        chromosome)
"""

merge_funcs = [
    read_protein_coding_genes_exons,
    read_protein_coding_genes_introns,
    read_non_protein_coding_genes_exons,
    read_non_protein_coding_genes_introns,
    read_pseudo_genes_exons,
    read_pseudo_genes_introns,
]

merge_func_names = [
    "pc_exons",
    "pc_introns",
    "npc_exons",
    "npc_introns",
    "pseudo_exons",
    "pseudo_introns",
]

merge_funcs_loci = [
    read_protein_coding_genes_exons_l,
    read_protein_coding_genes_introns_l,
    read_non_protein_coding_genes_exons_l,
    read_non_protein_coding_genes_introns_l,
    read_pseudo_genes_exons_l,
    read_pseudo_genes_introns_l,
]

merge_funcs_Oligo = [
    Oligo.File.read_protein_coding_exons,
    Oligo.File.read_protein_coding_genes_introns,
    Oligo.File.read_non_protein_coding_genes_exons,
    Oligo.File.read_non_protein_coding_genes_introns,
    Oligo.File.read_pseudo_genes_exons,
    Oligo.File.read_pseudo_genes_introns,
]



#----------- TOTAL GENOME ------------


def read_total_genome(chromosome):
    return None