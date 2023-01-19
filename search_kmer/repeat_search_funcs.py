#Functions to search kmer Spectra in elements of repeatmasker Files

#interesting categories:

#interesting elements:

#b) LTR, L1

#organisms:
#

import enum
import sys
sys.path.insert(1, r'E:\Genomik\Oligo')                     #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')          #Linux
import Oligo

sys.path.append('../..')                                    #add parent folder
import kmer_tools

from tqdm import tqdm
import os
try:
    from repeat_ids import *
except:
    from .repeat_ids import *
try:
    from loci_save_lib import *
except:
    from .loci_save_lib import *

    
repeat_folder = "../../data/repeatmasker/"
loci_folder = "../../data/loci/"

def filter_seq_id(row,id):
    return row[4] == id

def filter_repeat_class(row,label):
    return row[10][0:len(label)] == label

def filter_repeat_subclass_include(row,labels_include = None):
    if labels_include is None: return True
    sub_label = row[9]
    for label_in in labels_include:
        if (len(sub_label) >= len(label_in) and sub_label[0:len(label_in)] == label_in): return True
    return False

def filter_repeat_subclass_exclude(row,labels_exclude = None):
    if labels_exclude is None: return True
    sub_label = row[9]
    for label_ex in labels_exclude:
        if (len(sub_label) >= len(label_ex) and sub_label[0:len(label_ex)] == label_ex): return False
    return True

def filter_class(row,id,label):
    #ret = filter_seq_id(row,id) and filter_elemen_type(row,label)
    #if (ret): print([(i,entry) for i,entry in enumerate(row)])
    return filter_seq_id(row,id) and filter_repeat_class(row,label)

def filter_subclass(row, id, class_label, sub_label_in, sub_label_ex):
    #ret = (filter_seq_id(row,id) 
    #    and filter_repeat_class(row,class_label) 
    #    and filter_repeat_subclass_include(row,sub_label_in) 
    #    and filter_repeat_subclass_exclude(row,sub_label_ex))
    #if (ret): print([(i,entry) for i,entry in enumerate(row)])
    return (
        filter_seq_id(row,id) 
        and filter_repeat_class(row,class_label) 
        and filter_repeat_subclass_include(row,sub_label_in) 
        and filter_repeat_subclass_exclude(row,sub_label_ex)
    )


def read_repeat_loci_class(chromosome,repeat_class,verbose=1):
    if(verbose):
        print("reading "+repeat_class+" from "+str(chromosome))
    domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    if  (domain == "animalia"):      domain_ids = animalia_ids
    elif(domain == "embryophyta"):  domain_ids = embryophyta_ids
    elif(domain == "fungi"):        domain_ids = fungi_ids
    else: return None
    loci = Oligo.Locus.read_repeatmasker(repeat_folder+domain+"/"+domain_ids[chromosome.orga.name]+".out", filter_func=filter_class, filter_args=(chromosome.id,repeat_class,))
    if(verbose):
        print("found "+str(len(loci))+" "+repeat_class+" in "+str(chromosome))
    return loci

def read_repeat_loci_subclass(chromosome,repeat_class,subclasses_in,subclasses_ex,verbose=1):
    if(verbose):
        print("reading "+repeat_class + " - " , subclasses_in , " excluding " , subclasses_ex , " from "+str(chromosome))
    domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    if  (domain == "animalia"):      domain_ids = animalia_ids
    elif(domain == "embryophyta"):  domain_ids = embryophyta_ids
    elif(domain == "fungi"):        domain_ids = fungi_ids
    else: return None
    loci = Oligo.Locus.read_repeatmasker(repeat_folder+domain+"/"+domain_ids[chromosome.orga.name]+".out", filter_func=filter_subclass, filter_args=(chromosome.id,repeat_class,subclasses_in,subclasses_ex,))
    if(verbose):
        print("found "+str(len(loci))+" "+repeat_class+" in "+str(chromosome))
    return loci

def read_loci_file(chromosome,label,verbose=1):
    if(verbose):
        print("reading "+label+" from "+str(chromosome))
    domain = kmer_tools.download_genbank.name_to_domain(chromosome.orga)
    loci = read_loci(chromosome,domain,label)
    if(verbose):
        print("found "+str(len(loci))+" "+label+" in "+str(chromosome))
    return loci

#reading reapeat elements from saved loci file
#functions:
#   a) LTR:     LTR
#   b) SINE:    SINE
#               SINE/Alu
#               (SINE/B1)
#               SINE/B2
#               SINE/B4
#               SINE/ID
#               SINE/tRNA-C
#   c) LINE:    LINE
#               LINE/L1
#               LINE/L2

def read_LTR_f(chromosome):
    return read_loci_file(chromosome,"LTR")

def read_SINE_f(chromosome):
    return read_loci_file(chromosome,"SINE")

def read_SINE_Alu_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu")

def read_SINE_Alu_pure_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_pure")

def read_SINE_Alu_FAM_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_FAM")

def read_SINE_Alu_FRAM_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_FRAM")

def read_SINE_Alu_FLAM_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_FLAM")

def read_SINE_Alu_PB1_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_PB1")

def read_SINE_Alu_B1_f(chromosome):
    return read_loci_file(chromosome,"SINE_Alu_B1")

def read_SINE_B2_f(chromosome):
    return read_loci_file(chromosome,"SINE_B2")

def read_SINE_B4_f(chromosome):
    return read_loci_file(chromosome,"SINE_B4")

def read_SINE_ID_f(chromosome):
    return read_loci_file(chromosome,"SINE_ID")

def read_CSINE_f(chromosome):
    return read_loci_file(chromosome,"CSINE")

def read_LINE_f(chromosome):
    return read_loci_file(chromosome,"LINE")

def read_LINE_L1_f(chromosome):
    return read_loci_file(chromosome,"LINE_L1")

def read_LINE_L2_f(chromosome):
    return read_loci_file(chromosome,"LINE_L2")

def read_LINE_L3_f(chromosome):
    return read_loci_file(chromosome,"LINE_L3")

repeat_funcs_labels = [
    "LTR",
    "SINE",
    "SINE_Alu",
    "SINE_Alu_pure",
    "SINE_Alu_FAM",
    "SINE_Alu_FRAM",
    "SINE_Alu_FLAM",
    "SINE_Alu_PB1",
    "SINE_Alu_B1",
    "SINE_B2",
    "SINE_B4",
    "SINE_ID",
    "CSINE",
    "LINE",
    "LINE_L1",
    "LINE_L2",
    "LINE_L3",
]

repeat_funcs_read = [
    read_LTR_f,
    read_SINE_f,
    read_SINE_Alu_f,
    read_SINE_Alu_pure_f,
    read_SINE_Alu_FAM_f,
    read_SINE_Alu_FRAM_f,
    read_SINE_Alu_FLAM_f,
    read_SINE_Alu_PB1_f,
    read_SINE_Alu_B1_f,
    read_SINE_B2_f,
    read_SINE_B4_f,
    read_SINE_ID_f,
    read_CSINE_f,
    read_LINE_f,
    read_LINE_L1_f,
    read_LINE_L2_f,
    read_LINE_L3_f,
]

#reading reapeat elements from Repeatmasker file

def read_LTR(chromosome):
    return read_repeat_loci_class(chromosome,"LTR",verbose=1)

def read_SINE(chromosome):
    return read_repeat_loci_class(chromosome,"SINE")

def read_SINE_Alu(chromosome):
    return read_repeat_loci_class(chromosome,"SINE/Alu")


def read_SINE_Alu_pure(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = None
    subclass_ex = [
        "FAM",
        "FLAM",
        "FRAM",
        "PB1",
        "B1",
    ]
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)

def read_SINE_Alu_FAM(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = ["FAM"]
    subclass_ex = None
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)

def read_SINE_Alu_FRAM(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = ["FRAM"]
    subclass_ex = None
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)

def read_SINE_Alu_FLAM(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = ["FLAM"]
    subclass_ex = None
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)

def read_SINE_Alu_PB1(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = ["PB1"]
    subclass_ex = None
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)


def read_SINE_Alu_B1(chromosome):
    class_lb = "SINE/Alu"
    subclass_in = ["B1"]
    subclass_ex = None
    return read_repeat_loci_subclass(chromosome,class_lb,subclass_in,subclass_ex)

def read_SINE_B2(chromosome):
    return read_repeat_loci_class(chromosome,"SINE/B2")

def read_SINE_B4(chromosome):
    return read_repeat_loci_class(chromosome,"SINE/B4")

def read_SINE_ID(chromosome):
    return read_repeat_loci_class(chromosome,"SINE/ID")

def read_CSINE(chromosome):
    return read_repeat_loci_class(chromosome,"SINE/tRNA-C")

def read_LINE(chromosome):
    return read_repeat_loci_class(chromosome,"LINE")

def read_LINE_L1(chromosome):
    return read_repeat_loci_class(chromosome,"LINE/L1")

def read_LINE_L2(chromosome):
    return read_repeat_loci_class(chromosome,"LINE/L2")

def read_LINE_L3(chromosome):
    return read_repeat_loci_class(chromosome,"LINE/CR1")

repeat_funcs_save = [
    read_LTR,
    read_SINE,
    read_SINE_Alu,
    read_SINE_Alu_pure,
    read_SINE_Alu_FAM,
    read_SINE_Alu_FRAM,
    read_SINE_Alu_FLAM,
    read_SINE_Alu_PB1,
    read_SINE_Alu_B1,
    read_SINE_B2,
    read_SINE_B4,
    read_SINE_ID,
    read_CSINE,
    read_LINE,
    read_LINE_L1,
    read_LINE_L2,
    read_LINE_L3,
]

repeat_v1_indeces = [0,1,2,9,10,11,12,13,14,15]
repeat_v2_indeces = [3,4,5,6,7,8,16]

domain = "animalia"
organism_names = ["Homo sapiens"]
seqs_folder = "../download_genbank/"
k_values = [2]
stop_overwrite = True
labels = ["LTR","LINE"]



#output_filepath="../../data/k_spectra/"+domain+"/"
#organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
#for orga_name in tqdm(organism_names):
#    genome = Oligo.File.read_genome(orga_name,seqs_filename=Oligo.File.search(seqs_folder+domain+".seqs"))
#    for j,chromo in enumerate(genome):
        #chromo.load_seq()
#        print(chromo.id)
#        for i,label in tqdm(enumerate(labels)):
#            if (stop_overwrite): #skip loop if files already exist
#                file_count = 0
                #for k in k_values:
                    #output_filename = str(chromo)+"_"+function_names[i]+"_k="+str(k)+".kmer"
                    #if (os.path.isfile(output_filepath+output_filename)): file_count += 1
#                if (file_count == len(k_values)): continue
#            loci = Oligo.Locus.read_repeatmasker(repeat_folder+domain+"/"+animalia_ids["Homo sapiens"]+".out", filter_func=filter_combined, filter_args=(chromo.id,label,))
#            loci = Oligo.Loci.merge(loci)
#            print(len(loci))
            #for k in k_values:
                #kmer_func_search2(domain,chromo,k,loci,function_names[i],output_filepath = output_filepath,Windows=Windows)
        #chromo.unload_seq()

#loci = Oligo.Locus.read_repeatmasker(repeat_folder+domain+"/"+animalia_ids["Homo sapiens"]+".out", filter_func=lambda row,c: row[4] == c, filter_args=(chromo.id,"LTR"))