import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
sys.path.append('../..')                                    #add parent folder

import Oligo
import kmer_tools

from kmer_tools.download_genbank.orga_names_sorted import *
from kmer_tools.search_kmer.kmer_search_funcs import *
from kmer_tools.search_kmer.reorganize_funcs import *
from kmer_tools.search_kmer.repeat_search_funcs import *
from kmer_tools.search_kmer.loci_save_lib import *
 









Mamalia2_names = [
    "Homo sapiens",
    "Gorilla gorilla",
    "Pongo abelii",
    #"Pan troglodyte",
    "Macaca mulatta",
    "Rattus norvegicus",
    "Mus musculus",
    "Oryctolagus cuniculus",
    "Microtus ochrogaster",
    "Ochotona princeps",
    "Microcebus murinus",
]

seqs_folder = "../download_genbank/"

domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids"]
domain = "animalia"

organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
domain_names = [orga.name for orga in organisms]

genome_HS = Oligo.File.read_genome("Homo sapiens",seqs_filename=Oligo.File.search(seqs_folder+"animalia"+".seqs"))
chromosome = genome_HS[0]
loci = loci=Oligo.File.read_genes(chromosome=chromosome)
output_filepath = ""
output_filename =str(chromosome)+"_"+"genes k=2 test strand"+".kmer"
k_value = 2
Oligo.Search.search_kmer(
        data_seq=chromosome.get_seq(), 
        k=k_value,
        output_filename=output_filepath+output_filename,
        loci=loci,
        clear_file=True,
        regard_strand=True
        )

for locus in loci[0:100]:
    print(locus,locus.strand)
quit()
