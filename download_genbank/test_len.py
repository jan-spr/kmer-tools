from animalia_accessions import *
from fungi_accessions import *
from unknown_accessions import *

length=0
n=0
for group_acessions in fungi_acessions:
    for organism_acessions in group_acessions:
        length += len(organism_acessions)
        n += 1

print("total length:", length)
print("average length for fungi:", length/n)