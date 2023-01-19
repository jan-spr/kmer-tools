# few helper functions for reading kmer data from k_dat files

import sys

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

import numpy as np
from tqdm import tqdm

seqs_folder     = "../download_genbank/"
spec_folder     = "../../data/k_spectra/"
k_data_folder   = "../../data/k_data_strand/"

def get_orga_kmer_spectra(orga_names,func_name,func_name_bak,k_value,domain=None,domain_dict=None,domain_list=None,fix_empty = True):
    spectra_orgas =[]
    for i,orga_name in enumerate(tqdm(orga_names)):
        if (domain_list is not None):
            domain = domain_list[i]
        elif (domain_dict is not None):
            domain = domain_dict[orga_name]
        #print(orga_name,domain)
        genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
    
        #spectra_chrom.append(Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k_value)).replace('?','_').replace(':','_'),verbose=0))
        
        spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga_name)
        
        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        if (fix_empty and (spec.k is None)):
            kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(genome[0]),func_name_bak,k_value)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
            spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        #print("test: k= ",spec.get_kmers(),[spec.get_count(kmer) for kmer in spec.get_kmers()])
        spectra_orgas.append(spec)

    # Correlate
    #print(spectra_orgas)
    #print([[spec.get_count(kmer) for kmer in spec.get_kmers()] for spec in spectra_orgas])
    return spectra_orgas

def get_orga_kmer_heatmap(orga_names,func_name,func_name_bak,k_value,domain=None,domain_dict=None,fix_empty = True): # symmertic heatmap only (1 func)
    """
    spectra_orgas =[]
    for orga_name in tqdm(orga_names):
        if(domain_dict is not None):
            #print(domain_dict[orga_name])
            genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain_dict[orga_name]+".seqs"),verbose=0)
        else:
            genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
    
        #spectra_chrom.append(Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k_value)).replace('?','_').replace(':','_'),verbose=0))
        if(domain_dict is not None):
            spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain_dict[orga_name],func_name,k_value,orga_name=orga_name)
        else:
            spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga_name)

        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        if (fix_empty and (spec.k is None)):
            
            if(domain_dict is not None):
                kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain_dict[orga_name]+'/%s_%s_k=%d.kmer' % (str(genome[0]),func_name_bak,k_value)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
            else:
                kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(genome[0]),func_name_bak,k_value)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
            spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        #print("test: k= ",spec.get_kmers(),[spec.get_count(kmer) for kmer in spec.get_kmers()])
        spectra_orgas.append(spec)
    """
    spectra_orgas = get_orga_kmer_spectra(orga_names,func_name,func_name_bak,k_value,domain=domain,domain_dict=domain_dict,fix_empty = fix_empty)
    # Correlate
    #print(spectra_orgas)
    #print([[spec.get_count(kmer) for kmer in spec.get_kmers()] for spec in spectra_orgas])
    matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_orgas)
    return matrix

def get_func_kmer_heatmap(orga_name,func_names,func_name_bak,k_value,domain=None,domain_dict=None,fix_empty = True):
    spectra_funcs =[]
    for func_name in tqdm(func_names):
        if(domain_dict is not None):
            #print(domain_dict[orga_name])
            genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain_dict[orga_name]+".seqs"),verbose=0)
        else:
            genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"),verbose=0)
    
        #spectra_chrom.append(Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(chromo),func_name,k_value)).replace('?','_').replace(':','_'),verbose=0))
        if(domain_dict is not None):
            spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain_dict[orga_name],func_name,k_value,orga_name=orga_name)
        else:
            spec = kmer_tools.search_kmer.read_kdat_orgaSpec(k_data_folder,domain,func_name,k_value,orga_name=orga_name)

        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        if (fix_empty and (spec.k is None)):
            
            if(domain_dict is not None):
                kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain_dict[orga_name]+'/%s_%s_k=%d.kmer' % (str(genome[0]),func_name_bak,k_value)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
            else:
                kmers = Oligo.Kmer.KmerSpectrum.read((spec_folder+domain+'/%s_%s_k=%d.kmer' % (str(genome[0]),func_name_bak,k_value)).replace('?','_').replace(':','_'),verbose=0).get_kmers()
            spec = Oligo.Kmer.KmerSpectrum(kmers=kmers,values=np.zeros(len(kmers)),k=k_value)
        #print("test: spec= ",spec)
        #print("test: k= ",spec.k)
        #print("test: k= ",spec.get_kmers(),[spec.get_count(kmer) for kmer in spec.get_kmers()])
        spectra_funcs.append(spec)

    # Correlate
    #print(spectra_orgas)
    #print([[spec.get_count(kmer) for kmer in spec.get_kmers()] for spec in spectra_orgas])
    matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra_funcs)
    return matrix