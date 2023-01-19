#Functions to search kmer Spectra in different Domains, Organisms and genome elemts/areas

#interesting categories:
# pc, npc, pseudo, integenic
# subgroups:
# intron, exon

#exotic interesting elements:
#a) RNAs: lncRNA, miRNA, snoRNA, snRNA, precursorRNA
#b) LTR, L1
#c) npc: tRNA,rRNA

#organisms:
#

import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                     #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')          #Linux
import Oligo
from tqdm import tqdm
import os

seqs_folder = "../download_genbank/"

def search_kmer_chromosomal(organisms,functions,function_names,domain,k_values):
    for i,func in tqdm(enumerate(functions)):
        for orga in tqdm(organisms):
                genome = Oligo.File.read_genome(orga.name,seqs_filename=Oligo.File.search(domain+".seqs"))
                for k in k_values:
                    for j,chromo in enumerate(genome):
                        chromo.load_seq()

                        Oligo.Search.search_kmer(
                            data_seq=chromo.get_seq(), k=k,
                            data_seq_name=str(chromo)+"\t"+function_names[i], 
                            output_filename=('./data/k_mer/k_spectra/'+domain+'/%s_%s_k=%d.kmer' % (str(chromo),function_names[i],k)).replace('?','_'),
                            #loci=fix_loci_circular(loci=func(organism=orga,seqs_filename=Oligo.File.search(domain+".seqs"),verbose=0)[j],cutoff=chromo.length*0.8),
                            loci=None,
                            clear_file=True
                            )
                        chromo.unload_seq()

def chromo_gb_filename(chromosome,domain):
    return "../../data/genbank/"+domain+"/"+chromosome.orga.name+"_"+chromosome.name+".gb"

def delete_file_windows(file,Windows=False):          # helper function to circumvent Windows filename bug
    if (Windows):
        try:
            os.remove(file+"")
        except:
            print("file delete failed")
        print("deleted file: ",file)

def fix_filename_windows(filename,Windows=False):         # helper function to circumvent Windows filename bug
    if (Windows): 
        try:
            os.replace(filename+"",filename)
        except:
            print("fix filename worked?")
        print("fix filename borked?")

def kmer_single_search(domain,chromosome,k_value,function,function_name,output_filepath,output_filename,data_seq_name=None):
    data_seq_name=str(chromosome)+"\t"+function_name
    loci=function(chromosome=chromosome)
    print(len(loci))
    #print("loci =",loci)
    chromosome.load_seq()
    Oligo.Search.search_kmer(
        data_seq=chromosome.get_seq(), 
        k=k_value,
        data_seq_name=data_seq_name, 
        output_filename=output_filepath+output_filename,
        loci=loci,
        clear_file=True
        )
    chromosome.unload_seq()

def kmer_func_search(domain,chromosome,k_value,functions,function_names,output_filepath,Windows=False):
    chromosome.load_seq()
    for i,func in tqdm(enumerate(functions)):
        output_filename =str(chromosome)+"_"+function_names[i]+"_k="+str(k_value)+".kmer"
        
        delete_file_windows(output_filepath+output_filename,Windows)
            
        data_seq_name=str(chromosome)+"\t"+function_names[i]
        loci=func(chromosome=chromosome)
        print(len(loci))
        Oligo.Search.search_kmer(
            data_seq=chromosome.get_seq(), 
            k=k_value,
            data_seq_name=data_seq_name, 
            output_filename=output_filepath+output_filename,
            loci=loci,
            clear_file=True
            )
        
        fix_filename_windows(output_filepath+output_filename,Windows)
            
    chromosome.unload_seq()

def kmer_func_search2(domain,chromosome,k_value,loci_func,function_name,output_filepath,Windows=False):
    output_filename =str(chromosome)+"_"+function_name+"_k="+str(k_value)+".kmer"
    delete_file_windows(output_filepath+output_filename,Windows)
    data_seq_name=str(chromosome)+"\t"+function_name

    Oligo.Search.search_kmer(
        data_seq=chromosome.get_seq(), 
        k=k_value,
        data_seq_name=data_seq_name, 
        output_filename=output_filepath+output_filename,
        loci=loci_func,
        clear_file=True,
        merge_loci=False,
        regard_strand = True
        )
    fix_filename_windows(output_filepath+output_filename,Windows)


def func_chromo_loop(organism_names,domain,functions,function_names,k_values,Windows=False): #loop over all orgas, k, chromos
    output_filepath="../../data/k_spectra/"+domain+"/"
    #organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga_name in tqdm(organism_names):
        genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
        for k in k_values:
            for j,chromo in enumerate(genome):
                kmer_func_search(domain,chromo,k,functions,function_names,output_filepath = output_filepath,Windows=Windows)
                


def func_loop_kmerOpt(organism_names,domain,functions,function_names,k_values,Windows=False,stop_overwrite=False): #loop optimized for laoding chromosomes and reading genome elements
    output_filepath="../../data/k_spectra/"+domain+"/"
    #organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    for orga_name in tqdm(organism_names):
        genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
        if(len(genome)==0):
            print("Error: orga not found")
            return None
        for j,chromo in enumerate(genome):
            chromo.load_seq()
            for i,func in tqdm(enumerate(functions)):
                print("looking for "+function_names[i])
                if (stop_overwrite): #skip loop if files already exist
                    file_count = 0
                    for k in k_values:
                        output_filename = str(chromo)+"_"+function_names[i]+"_k="+str(k)+".kmer"
                        if (os.path.isfile(output_filepath+output_filename)): file_count += 1
                    if (file_count == len(k_values)): continue
                print(func)
                loci = func(chromosome=chromo)
                #print(loci)
                if loci is not None:
                    loci = Oligo.Loci.merge(loci)
                    print(len(loci))
                for k in k_values:
                    kmer_func_search2(domain,chromo,k,loci,function_names[i],output_filepath = output_filepath,Windows=Windows)
            chromo.unload_seq()

def fix_kmer_filename_windows(folder):
    filenames = next(os.walk(folder), (None, None, []))[2]  # [] if no file
    print(filenames[0:30])
    for filename in tqdm(filenames):
        os.replace(folder+filename,folder+(filename.replace('','_').replace('','_')))
