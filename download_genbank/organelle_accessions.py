#acession ids for all organelle of analyzed organisms

#1st: Animalia - Mitochondria:

#Mammalia1:
#Homo sapiens, Pan troglodytes, Mus musculus, Rattus norvegicus, Oryctolagus cuniculus, Canis lupus, Felis catus, Equus caballus, Capra hircus, Bos taurus, Sus scrofa, Monodelphis domestica, Ornithorhynchus anaticus,
#Aves:
#Gallus gallus, Taeniopygia guttata, Ficedula albicollis, Strigops habroptila, Falco rusticolus, Aquila chrysaetos chrysaetos, Calypte anna, 
#Archelosauria:
#Chrysemys picta
#Lepidosauria:
#Anolis carolinensis, Lacerta agilis
#Amphibia:
#Xenopus tropicalis, Rana temporaria, Geotrypetes seraphini
#Actinopterygii:
#Takifugu rubripes, Cynoglossus semilaevis, Oreochromis niloticus, Poecilia reticulata, Danio rerio, Salmo salar, Lepisosteus oculatus, 
#Echinodermata:
#Asterias rubens
#Tunicata:
#Ciona intestinalis
#Spiralia:
#Schistosoma mansoni
#Ecdysozoa:
#Strongyloides ratti, Caenorhabditis brigsae, Caenorhabditis elegans, 
#Mollusca:
#Gigantopelta aegis, Pomacea canaliculata, Pecten maximus, Crassostrea virginica, Octopus sinensis, 
#Crustacea:
#Lepeophtheirus salmonis, Penaeus monodon
#Chelicerata:
#Rhipicephalus sanguineus
#Insecta:
#Rhopalosiphum maidis, Apis mellifera, Bombus terrestris, Nasonia vitripennis, Anopheles gambiae, Drosphila pseudoobscura, Drosophila melanogaster, Drosohpila simulans, Aricia agestis, 


""""
Missing Mitochondria data: 
Canis lupus
Strigops habroptila
Aquila chrysaetos chrysaetos
Calypte anna
Anolis carolinensis
Rana temporaria
Geotrypetes seraphini
Takifugu rubripes
Salmo salar
Strongyloides ratti
Caenorhabditis brigsae

Gigantopelta aegis
Pecten maximus
Crassostrea virginica
Lepeophtheirus salmonis

Rhopalosiphum maidis
Bombus terrestris
Nasonia vitripennis
Drosohpila simulans
Aricia agestis

Microtus ochrogaster

Dictyostelium discoideum???

manually added: 


"""
mammalia1_mitochon_accessions = [
    "NC_012920.1",
    "NC_001643.1",
    "NC_005089.1",
    "NC_001665.2",
    "NC_001913.1",
    "NC_001700.1",
    "NC_001640.1",
    "NC_005044.2",
    "NC_006853.1",
    "NC_000845.1",
    "NC_006299.1",
    "NC_000891.1",
]

aves_mitochon_accessions = [
    "NC_053523.1",
    "NC_007897.1",
    "NC_021621.1",
    "NC_029359.1",
]

fish_mitochon_accessions = [
    "NC_023890.1",
    "NC_021766.1",
    "NC_006839.1",
    "NC_012825.1",
    "NC_013663.1",
    "NC_024238.1",
    "NC_002333.2",
    "NC_004744.1",
    "NC_017929.1",
    "NC_002545.1",
]

mollusc_mitochon_accessions = [
    "NC_001328.1",
    "NC_024586.1",
    "NC_006353.1",
    "NC_002184.1",
    "NC_002074.1",
]

insecta_mitochon_accessions = [
	"NC_001566.1",
    "NC_002084.1",
    "NC_046603.1",
    "NC_024511.2",
]

#Mamalia2:
#Homo sapiens, Gorilla gorilla, Pongo abelii, Pan troglodyte, Macaca mulatta, Rattus norvegicus, Mus musculus, Oryctolagus cuniculus, Microtus ochrogaster, Ochotona princeps, Microcebus murinus,

mammalia2_mitochon_accessions = [
    "NC_012920.1",
    "NC_011120.1",
    "NC_002083.1",
    "NC_001643.1",
    "NC_005943.1",
    "NC_001665.2",
    "NC_005089.1",
    "NC_001913.1",
    "NC_005358.1",
    "NC_028718.1",
]



animalia_mitochondria_accessions = [
    mammalia1_mitochon_accessions,
    aves_mitochon_accessions,
    fish_mitochon_accessions,
    mollusc_mitochon_accessions,
    insecta_mitochon_accessions,
    mammalia2_mitochon_accessions,
]

#2nd: Embryophyta

#Bryophyta:
#Physcomitrella patens
#Tracheophyta—Liliopsida:
#Brachypodium distachyon, Oryza sativa, Zea mays, Sorghum bicolor, Setaria italic, Ananas comosus, Elaeis guineensis, Musa acuminate, Asparagus officinalis, 
#Tracheophyta—Ranunculales:
#Papaver somniferum
#Tracheophyta—Eudicotyledons:
#Vitis vinifera, Arabidopsis thaliana, Brassica napus, Camelina sativa, Gossypium raimondii, Theobroma cacao, Citrus sinensis, Populus trichocarpa, Manihot esculenta, Arachis duranensis, Phaseolus vulgaris, Cajanus cajan, Vigna angularis, Medicago truncatula, Lupinus angustifolius, Fragaria vesca, Malus domestica, Prunus mume, Cucumis sativus, Beta vulgaris, Nicotiana attenuata, Sesamum indicum, Helianthus annuus, 

Embryophyta_mitochon_accessions = [
    "NC_007945.1",
    "NC_011033.1",
    "NC_007982.1",
    "NC_008360.1",
	"NC_012119.1",
    "NC_037304.1",
    "NC_008285.1",
    "NC_045136.1",
    "NC_021092.1",
    "NC_029641.1",
    "NC_018554.1",
    "NC_023337.1",
    "NC_004946.1",
    
]

Embryophyta_plastid_accessions = [
    "NC_005087.1",
    "NC_011032.1",
    "NC_001320.1",
	"NC_001666.2",
    "NC_008602.1",
    "NC_022850.1",
    "NC_026220.1",
    "NC_017602.1",
    "NC_029434.1",
    "NC_007957.1",
    "NC_000932.1",
        
    "NC_014676.2",
    "NC_008334.1",
    "NC_009143.1",
    "NC_010433.1",
    "NC_031429.1",
    "NC_021091.1",
    "NC_003119.8",
	"NC_015206.1",
    "NC_023798.1",
    "NC_007144.1",
    "NC_016433.2",
    "NC_007977.1",

    "NC_016734.1",

]

Embryophyta_etc_accessions = [
    #"NC_016734.1",
]

#3rd: Protista

#Amoebozoa:
#Dictyostelium discoideum
#Euglenozoa:
#Leishmania donovani, Leishmania infantum, Leishmania major, Trypanosoma brucei
#Alveolata:
#Neospora caninum, Toxoplasma gondii, Cryptosporidium parvum, Plasmodium vivax, Plasmodium falciparum, Plasmodium berghei, Plasmodium cynomolgi, Plasmodium coatneyi, Plasmodium gaboni, Theileria annulata, Theileria parva, Theileria orientalis, Babesia microti, Babesia bigemina, Paramecium tetraurelia, 
#Stramenopiles:
#Phaeodactylum tricornutum, Ectocarpus siliculosus, Thalassiosira pseudonana, Nannochloropsis gaditana
#Rhodophyta:
#Cyanidioschyzon merolae
#Chlorophyta:
#Chlamydomonas reinhardtii, Micromonas commoda, Bathycoccus prasinos, Ostreococcus lucimarinus, Ostreococcus tauri

Protista_mitochon_accessions = [
    "NC_000895.1",
    "NC_007243.1",
    "NC_015303.1",
    "NC_034637.1",
    "NC_012643.1",
    "NC_023273.1",
    "NC_008290.1",
]

Protista_plastid_accessions = [
    "NC_001799.1",
    "NW_017385168.1",
    "NC_007758.1",
    "NC_034636.1",
	"NC_012575.1",
    "NC_024811.1",
    "NC_008289.1",
]


#3th: Fungi

#Microsporidia:
#Encephalitozoon cuniculi, Encephalitozoon intestinalis
#Basidiomycota—Ustilaginomycetes:
#Ustilago maydis, Sporisorium reilianum,
#Basidiomycota—Tremellomycetes:
#Cryptococcus gattii
#Basidiomycota—Agaricomycetes:
#Pyrrhoderma noxium
#Ascomycota—Schizosaccharomycetes:
#Schizosaccharomyces pombe, Saccharomyces cerevisiae, Tetrapisispora phaffii, Tetrapisispora blattae, Zygosaccharomyces rouxii, Naumovozyma castellii, Lachancea thermotolerans, Eremothecium gossypii, Candida glabrata, Candida albicans, Debaryomyces hansenii, Komagataella phaffii, Yarrowia lipolytica, Pichia kudriavzevii, Ogataea parapolymorpha, 
#Ascomycota—Dothideomycetes:
#Zymoseptoria tritici, 
#Ascomycota—Sordariomycetes:
#Neurospora crassa, Thermothelomyces thermophila, Thielavia terrestris, Fusarium oxysporum, 
#Ascomycota—Eurotiomycetes:
#Aspergillus nidulans, Aspergillus fumigatus, Penicillium chrysogenum, 

Fungi_mitochon_accessions = [
    "CM008263.2",
    "NC_001326.1",
    "NC_001224.1",
    "NC_005789.1",
    "NC_004691.1",
	"NC_010166.1",
    "NC_002659.1",
    "NC_026614.1",
    "NC_017896.1",
    "CM002802.1",
]

#print(animalia_mitochondria_accessions)
animalia_mitochondria_accessions = [j for i in animalia_mitochondria_accessions for j in i]
#print(animalia_mitochondria_accessions)

mitochondria_accessions = [
    animalia_mitochondria_accessions,
    Embryophyta_mitochon_accessions,
    Protista_mitochon_accessions,
    Fungi_mitochon_accessions
    ]

plastid_accessions = [
    Embryophyta_plastid_accessions,
    Protista_plastid_accessions
    ]