import numpy as np

def word_num(row,words):
    counts = []
    for i , word in enumerate(words):
        counts.append(row[word])
    return counts, np.zeros(len(counts))

def dimer_count_normalized(row,words):
    counts= []
    for i,word in enumerate(words):
        counts.append( row[word]/(row["A"]+row["T"]+row["C"]+row["G"]+0.001)*4**len(word))
    return counts, np.zeros(len(counts))

def dimer_portion_normalized(row,words):
    counts= []
    for i,word in enumerate(words):
        counts.append( row[word]/(row["A"]+row["T"]+row["C"]+row["G"]+0.001)*4**len(word)-1 )
    return counts, np.zeros(len(counts))

def dimer_minus_portion(row,words):
    counts= []
    for i,word in enumerate(words):
        counts.append( row[word] - (row["A"]+row["T"]+row["C"]+row["G"])/(4**len(word)) )
    return counts, np.zeros(len(counts))

def dimer_minus_expectation(row,words):
    length = row["A"]+row["T"]+row["C"]+row["G"]+0.001
    counts = []
    err    = []
    for i,word in enumerate(words):
        monomers = [word[0],word[1]]
        p_x = [row[letter]/length for letter in monomers]
        p_xy = np.prod(p_x)
        E_xy = p_xy * length
        sigma2 = length * p_xy * (1- p_xy)
        C = row[word] -E_xy
        counts.append(C)
        err.append(np.sqrt(sigma2))
    return counts,err

def dimer_minus_expectation_norm(row,words):
    length = row["A"]+row["T"]+row["C"]+row["G"]+0.001
    n = (length -1 )/16
    counts = []
    err    = []
    for i,word in enumerate(words):
        monomers = [word[0],word[1]]
        p_x = [row[letter]/length for letter in monomers]
        p_xy = np.prod(p_x)
        #print(monomers,p_x,p_xy)
        E_xy = p_xy * length
        sigma2 = length * p_xy * (1- p_xy)
        C = (row[word] -E_xy)/n
        counts.append(C)
        err.append(np.sqrt(sigma2)/n)
    return counts,err