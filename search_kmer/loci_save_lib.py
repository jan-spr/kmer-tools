#useful functions for reorganizing data in kmer analysis
#2) functions to read and write loci to files

import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                     #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')          #Linux
import Oligo

sys.path.append('../..')                                    #add parent folder
import kmer_tools

import os
from tqdm import tqdm

seqs_folder = "../download_genbank/"
loci_folder = "../../data/loci/"

def chromo_fix_filename(filename):
    return filename.replace('','_').replace('','_').replace('?','_')

def read_loci(chromosome,domain,func_name):
    loci = Oligo.Locus.read(chromo_fix_filename(loci_folder+domain+"/"+str(chromosome)+"_"+func_name+".loci"))
    return loci

def save_loci(chromosome,domain,function,func_name,stop_overwrite=False):
    #if (domain is None): domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    output_filename = chromo_fix_filename(loci_folder+domain+"/"+str(chromosome)+"_"+func_name+".loci")
    if (stop_overwrite): 
        if (os.path.isfile(output_filename)): 
            print ("file already exists")
            return None
    loci = function(chromosome=chromosome)
    Oligo.Loci.Locus.save(loci,output_filename)

def read_loci(chromosome,domain,func_name):
    input_filename = chromo_fix_filename(loci_folder+domain+"/"+str(chromosome)+"_"+func_name+".loci")
    loci=Oligo.Locus.read(input_filename)
    return loci

def func_loop(organism_names,domain,action_func,function_names,function_dict=None,functions=None,Windows=False,stop_overwrite=False): #loop over chromosomes+functions and do given unspecific task
    if function_dict == None:
        if functions == None:
            print("Error: fucntions not defined")
            return None
        function_dict=dict(zip(function_names,functions))
    print(function_dict)
    for orga_name in tqdm(organism_names):
        genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
        if(len(genome)==0):
            print("Error: orga not found")
            return None
        for j,chromo in enumerate(genome):
            for i,func_name in enumerate(function_names):
                print("function: "+func_name)
                action_func(chromo,domain,function_dict[func_name],func_name,stop_overwrite=stop_overwrite)
