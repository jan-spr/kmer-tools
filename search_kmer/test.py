import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
sys.path.append('../..')                                    #add parent folder

import Oligo
import kmer_tools

from kmer_tools.download_genbank.orga_names_sorted import *
from kmer_search_funcs import *
from reorganize_funcs import *
from repeat_search_funcs import *
from loci_save_lib import *
from kmer_tools.plot.dimer_func_lib import *
 
functions=[
    #Oligo.File.read_genes,
    #Oligo.File.read_exons,
    #Oligo.File.read_introns,
    #Oligo.File.read_cds,
    Oligo.File.read_mRNAs,
    #Oligo.File.read_intergenics,
    #Oligo.File.read_protein_coding_genes,
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
    #'genes',
    #'exons',
    #'introns',
    #'cds',
    'mRNAs',
    #'intergenics',
    #'protein_coding_genes',
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
    read_LTR,
    read_SINE,
    read_SINE_Alu,
    read_LINE,
    read_LINE_L1,
    read_LINE_L2,
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
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

functions_repeat_dict={}
for func,name in zip(functions_repeat,function_names_repeat):
    functions_repeat_dict[name] = func

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


domain = "animalia"

seqs_folder = "../download_genbank/"
entries = ["chromosome","length","G+C"]

domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids"]
domain = "animalia"

#organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
#for orga_group in [Fungi_names_Basidiomycota_p,Fungi_names_Ascomycota]:

organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
domain_names = [orga.name for orga in organisms]




animalia_names_repeat = [
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

animalia_names_fix = [
    "Caenorhabditis briggsae",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Strongyloides ratti",
]

repeat_function_names=[
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

for dict,domain in [
    (animalia_ids,"animalia"),
    (embryophyta_ids,"embryophyta"),
    (fungi_ids,"fungi")
    ]:      
    orga_names = [name for name in dict] 
    orga_names = [orga_name_NCBI_dict[name] for name in orga_names]
    #func_loop_kmerOpt(orga_names,domain=domain,functions = functions,function_names = function_names,k_values=[1,2,3],Windows=False,stop_overwrite=True)
    #func_loop(orga_names,domain,save_loci,functions_repeat_dict,function_names,Windows=False,stop_overwrite=True)
print(domain)
#store_kmer_chromosomal(name_to_orga(domain_names,organisms),function_names,domain,k_values=[1,2,3])
#fix_kmer_filename_windows("../../data/k_spectra/"+domain+"/")

orga_names = [
    "Homo sapiens"
]
domain = "animalia"
file=seqs_folder+domain+".seqs"
organisms = kmer_tools.download_genbank.names_to_orgas(orga_names)
seq = "GTTGCCCAGGGCTCACCTCTGATGGGTGGGCAGGTGAGCGCTTCCAACAGCTTCTCGAGGCTGCACTGCAGAAATGCCAACGAGGACTGGATGTCGGCACTGTGTCCCCGGCTCTGGGATGTGCCCCTCCACCACCTCTCCATCCCAG"
#print("length:",len(seq))
#Oligo.Search.search_kmer(seq, k=2, output_filename='test.kmer')
orga_names = [name for name in animalia_ids] 
domain = "animalia"
for name in orga_names[0::3]:
    #organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    file=seqs_folder+domain+".seqs"
    genome = Oligo.File.read_genome(name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
    for chromosome in genome[0:1]:
        #print(str(chromosome))
        #print("exon")
        #loci = read_SINE_f(chromosome = chromosome)
        #data_seq=chromosome.get_seq()
        #for i,locus in enumerate(loci[0:10]):
        #    print(i+1,data_seq[locus.start:locus.get_end()])
        #    pass
        pass
k_data_folder    = "../../data/k_data_strand/"        
func_pair=["npc_exons","intergenics"]
for domain in ["embryophyta"]:
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    orga_names = [orga.name for orga in organisms]
    corr_vals=[]
    for orga_name in orga_names:
        genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
        func_spec_pair =[]
        for k,func_name in enumerate(func_pair):
            spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value=3,orga_name=orga_name)
            
            row = kmer_tools.search_kmer.read_kdat_orgaRow(k_data_folder,domain,func_name,orga_name=orga_name)
            elem_len = row["A"] + row["T"]+ row["C"]+ row["G"]
            mark_entry =False
            
            if (spec.k is None or elem_len < 0.00001*row["length"]):
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
    print("mean value, error:",m,e)

for dict,domain in [
    #(animalia_ids,"animalia"),
    #(embryophyta_ids,"embryophyta"),
    #(fungi_ids,"fungi")
    ]:
    orga_names = [name for name in dict] 
    orga_name = "Pongo abelii"
    orga_names = ["Pongo abelii"]
    organisms = kmer_tools.download_genbank.names_to_orgas(orga_names)
    genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
    #read_repeat_loci(genome[0],repeat_class="SINE/Alu",verbose=1)

files = [
    r"E:\Genomik\bacteria and archaea\scripts\data\k_mer\k_data\archaea\kdata_archaea_intergenics_k=1-2.txt",
    r"E:\Genomik\bacteria and archaea\scripts\data\k_mer\k_data\archaea\kdata_archaea_protein_coding_genes_k=1-2.txt",
]


"""
for file in files:
    orga_name = "Acidilobus saccharovorans 345-15 ?"
    orga_row = {}
    f=pd.read_csv(file, sep="\t",header =0)
    for index, row in f.iterrows():
        #print(row["chromosome"])
        if (orga_name ==row["chromosome"]):
            orga_row = row

    di_dist=dimer_dist("Homo sapiens","cds")
    di_dist.row = orga_row
    di_dist.length = orga_row["A"] + orga_row["C"] + orga_row["G"] + orga_row["T"]
    print("results.")
    print (di_dist.get_func_diff("count","mono_expectation"))
"""
