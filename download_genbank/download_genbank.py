import sys
sys.path.insert(1, r'D:\Users\Jan Sprengel\Desktop\Genomik\Oligo')
import Oligo
from tqdm import tqdm

from animalia_accessions import *
from fungi_accessions import *
from unknown_accessions import *
from organelle_accessions import *
from nucleomorph_acessions import *


#Oligo.File.download_genbanks(TakifuguRubripes_acessions,seqs_filename='animalia.seqs',path='../../data/genbank/animalia/', labels=['Genome'])

#for group_acessions in tqdm(fungi_acessions):
#    for organism_acessions in tqdm(group_acessions):
#        while (True):
#            try:
#                print(organism_acessions)
#                Oligo.File.download_genbanks(organism_acessions,seqs_filename='fungi.seqs',path='../../data/genbank/fungi/', labels=['Genome'])
#            except:
#                continue
#            break

for organism_acessions in tqdm(nucleomorph_accessions):
    print(organism_acessions)
    Oligo.File.download_genbanks(organism_acessions,seqs_filename='nucleomorph.seqs',path='../../data/genbank/nucleomorph/', labels=['Genome'])
    while (True):
        try:
            print(organism_acessions)
        except:
            continue
        break

#Oligo.File.download_genbanks(PanTroglodytes_acessions,seqs_filename='nucleomorph.seqs',path='../../data/genbank/nucleomorph/', labels=['Genome'])