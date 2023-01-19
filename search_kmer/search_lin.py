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

#--------FUNCTIONS-----------

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
    read_LTR,
    read_SINE,
    read_SINE_Alu,
    read_LINE,
    read_LINE_L1,
    read_LINE_L2,
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

function_names_repeat=[
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

functions_merge_Oligo = kmer_tools.search_kmer.merge_funcs_Oligo
functions_merge_dict_Oligo={}
for func,name in zip(functions_merge_Oligo,function_names_merge):
    functions_merge_dict_Oligo[name] = func

functions_repeat_dict={}
for func,name in zip(functions_repeat,function_names_repeat):
    functions_repeat_dict[name] = func

function_names_repeat=[
    'LTR',
    'SINE',
    'SINE_Alu',
    'LINE',
    'LINE_L1',
    'LINE_L2',
]

search_kmer_funcs= functions + functions_merge_Oligo
search_kmer_funcs_names= function_names + function_names_merge


#--------ORGANISMS-----------


Mammalia1_names = [
    #"Homo sapiens",
    #"Pan troglodytes",
    #"Mus musculus",
    #"Rattus norvegicus",
    #"Oryctolagus cuniculus",
    "Canis lupus familiaris",
    #"Felis catus",
    #"Equus caballus",
    #"Capra hircus",
    #"Bos taurus",
    #"Sus scrofa",
    #"Monodelphis domestica",
    "Ornithorhynchus anatinus",
]

Mammalia2_names = [
    "Homo sapiens",
    "Gorilla gorilla",
    "Pongo abelii",
    "Pan troglodyte",
    "Macaca mulatta",
    "Rattus norvegicus",
    "Mus musculus",
    "Oryctolagus cuniculus",
    "Microtus ochrogaster",
    "Ochotona princeps",
    "Microcebus murinus",
]

Aves_names = [
    "Gallus gallus",
    "Taeniopygia guttata",
    "Ficedula albicollis",
    "Strigops habroptila",
    "Falco rusticolus",
    "Aquila chrysaetos chrysaetos",
    "Calypte anna",
]

Fungi_names_Basidiomycota_p = [
    "Encephalitozoon cuniculi",
    "Encephalitozoon intestinalis",
    "Ustilago maydis",
    "Sporisorium reilianum",
    "Cryptococcus gattii",
    "Pyrrhoderma noxium",
]
Fungi_names_Ascomycota = [
    "Schizosaccharomyces pombe",
    "Saccharomyces cerevisiae",
    "Tetrapisispora phaffii",
    "Tetrapisispora blattae",
    "Zygosaccharomyces rouxii",
    "Naumovozyma castellii",
    "Lachancea thermotolerans",
    "Eremothecium gossypii",
    "Candida glabrata",
    "Candida albicans",
    "Debaryomyces hansenii",
    "Komagataella phaffii",
    "Yarrowia lipolytica",
    "Pichia kudriavzevii",
    "Ogataea parapolymorpha",
    "Zymoseptoria tritici ",
    "Neurospora crassa",
    "Thermothelomyces thermophila",
    "Thielavia terrestris",
    "Fusarium oxysporum",
    "Aspergillus nidulans",
    "Aspergillus fumigatus",
    "Penicillium chrysogenum",
]

reptile_fish_names = [
    "Chrysemys picta",
    "Anolis carolinensis",
    "Lacerta agilis",
    "Xenopus tropicalis",
    "Rana temporaria",
    "Geotrypetes seraphini",
    "Takifugu rubripes",
    "Cynoglossus semilaevis",
    "Oreochromis niloticus",
    "Poecilia reticulata",
    "Danio rerio",
    "Salmo salar",
    "Lepisosteus oculatus",
]
crustacean_names = [
    "Asterias rubens",
    "Ciona intestinalis",
    "Schistosoma mansoni",
    "Strongyloides ratti",
    "Caenorhabditis brigsae",
    "Caenorhabditis elegans",
    "Gigantopelta aegis",
    "Pomacea canaliculata",
    "Pecten maximus",
    "Crassostrea virginica",
    "Octopus sinensis",
    "Lepeophtheirus salmonis",
    "Penaeus monodon",
    "Rhipicephalus sanguineus",
]
insect_names = [
    "Rhopalosiphum maidis",
    "Apis mellifera",
    "Bombus terrestris",
    "Nasonia vitripennis",
    "Anopheles gambiae",
    "Drosphila pseudoobscura",
    "Drosophila melanogaster",
    "Drosophila simulans",
    "Aricia agestis",
]

protista_names = [
    "Dictyostelium discoideum",
    "Leishmania donovani",
    "Leishmania infantum",
    "Leishmania major",
    "Trypanosoma brucei",
    "Neospora caninum",
    "Toxoplasma gondii",
    "Cryptosporidium parvum",
    "Plasmodium vivax",
    "Plasmodium falciparum",
    "Plasmodium berghei",
    "Plasmodium cynomolgi",
    "Plasmodium coatneyi",
    "Plasmodium gaboni",
    "Theileria annulata",
    "Theileria parva",
    "Theileria orientalis",
    "Babesia microti",
    "Babesia bigemina",
    "Paramecium tetraurelia",
    "Phaeodactylum tricornutum",
    "Ectocarpus siliculosus",
    "Thalassiosira pseudonana",
    "Nannochloropsis gaditana",
    "Cyanidioschyzon merolae",
    "Chlamydomonas reinhardtii",
    "Micromonas commoda",
    "Bathycoccus prasinos",
    "Ostreococcus lucimarinus",
    "Ostreococcus tauri",
]

embryophyta_names = [
    "Physcomitrium patens",
    "Brachypodium distachyon",
    "Oryza sativa",
    "Zea mays",
    "Sorghum bicolor",
    "Setaria italic",
    "Ananas comosus",
    "Elaeis guineensis",
    "Musa acuminate",
    "Asparagus officinalis",
    "Papaver somniferum",
    "Vitis vinifera",
    "Arabidopsis thaliana",
    "Brassica napus",
    "Camelina sativa",
    "Gossypium raimondii",
    "Theobroma cacao",
    "Citrus sinensis",
    "Populus trichocarpa",
    "Manihot esculenta",
    "Arachis duranensis",
    "Phaseolus vulgaris",
    "Cajanus cajan",
    "Vigna angularis",
    "Medicago truncatula",
    "Lupinus angustifolius",
    "Fragaria vesca",
    "Malus domestica",
    "Prunus mume",
    "Cucumis sativus",
    "Beta vulgaris",
    "Nicotiana attenuata",
    "Sesamum indicum",
    "Helianthus annuus",
]

domain = "embryophyta"

domains_fix=["animalia","embryophyta"]
names_fix = {}
names_fix["animalia"]= [
    "Caenorhabditis briggsae",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Strongyloides ratti",  
]

names_fix["embryophyta"]= [
    "Phaseolus vulgaris",
]


domains = ["animalia","embryophyta","fungi","protista","mitochondria","plastids","nucleomorph"]


#---------CODE---------


for domain in ["nucleomorph"]:
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    orga_names = [orga.name for orga in organisms]
    func_loop_kmerOpt(orga_names, domain=domain, functions = [kmer_tools.search_kmer.read_total_genome], function_names = ["total_genome"], k_values=[1,2,3,4,5,6], Windows=False, stop_overwrite=False)
    func_loop_kmerOpt(orga_names, domain=domain, functions = search_kmer_funcs, function_names = search_kmer_funcs_names, k_values=[1,2,3,4,5,6], Windows=False, stop_overwrite=False)
    #func_loop_kmerOpt(orga_names, domain=domain, functions = functions_merge_Oligo, function_names = function_names_merge, k_values=[1,2,3,4,5,6], Windows=False, stop_overwrite=False)
    #func_loop(orga_names[0:40],domain,save_loci,functions_merge_single_dict,function_names_merge_single,Windows=False,stop_overwrite=False)

print( [repeat_funcs_labels[i] for i in repeat_v2_indeces])
for dict,domain in [
    (animalia_ids,"animalia"),
    (embryophyta_ids,"embryophyta"),
    (fungi_ids,"fungi")
    ]:      
    orga_names = [name for name in dict] 
    #print(orga_names)
    #functions = [functions_repeat_dict[name] for name in function_names]
    organisms = Oligo.File.read_all_orgas(seqs_filename=seqs_folder+domain+".seqs")
    domain_names = [orga.name for orga in organisms]
    #func_loop(orga_names[30:40], domain, save_loci, function_names=repeat_funcs_labels, functions=repeat_funcs_save, Windows=False, stop_overwrite=False)
    #print(domain,": ",domain_names)
    #for names in [domain_names]:       
    #func_loop_kmerOpt(orga_names, domain=domain, functions = repeat_funcs_read, function_names = repeat_funcs_labels, k_values=[1,2,3,4,5,6], Windows=False, stop_overwrite=False)
    #store_kmer_chromosomal(organisms,function_names,domain,k_values=[1,2,3])
