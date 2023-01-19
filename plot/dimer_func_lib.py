import sys

sys.path.insert(1, r'E:\Genomik\Oligo')                 #Windows
sys.path.insert(1, '/home/chilly/Genomik/Oligo/')      #Linux
import Oligo

sys.path.append('..\..')
import kmer_tools

import numpy as np
from tqdm import tqdm
from scipy import stats as stat
k_data_folder   = "../../data/k_data/"
k_data_folder   = "../../data/k_data_strand/"
k_data_folder   = "../../data/k_data/"
#functions to compare dimer quantities between different Genome-elements and models as well
#functionalities:
#   -normalization:         divide total word-count by length of Genome-element to obtain coverage-independent proportion of each word
#   -deviation:             subtraction of normalized dimer distribution by normalized expectation to get a quantifiable deviation
#   -monomer expectation:   calculate expected quantity of dimers from proporion of monomers as independent combination




class dimer_dist:
    letters = ['A','T','C','G']
    dimer_words = ["AA",	"TT",	"AT",	"TA",	"AC",	"CA",	"TG",	"GT",	"AG",	"GA",	"TC",	"CT",	"CC",	"GG",	"GC",	"CG"]
    dimer_words = ["AA",	"TT",	"AT",	"TA",	"AC",	"GT",	"CA",	"TG",	"AG",	"CT",	"GA",	"TC",	"GC",	"CG",	"CC",	"GG"] #revCom together
    def __init__(self,orga_name,func_name,domain=None,total=True,chromo=False):
        self.orga_name = orga_name
        self.func_name = func_name
        if domain is None:
            domain = kmer_tools.download_genbank.name_to_domain(orga_name)
            print("automatic domain for "+orga_name+" = "+domain)
        self.domain = domain
        row = kmer_tools.search_kmer.read_kdat_dimerRow(k_data_folder,domain,func_name=func_name,orga_name = orga_name,chromo=chromo)
        if total:
            self._total_row = kmer_tools.search_kmer.read_kdat_dimerRow(k_data_folder,domain,func_name="total_genome",orga_name = orga_name,chromo=chromo)
        self.row = row
        self.chromo_length = row["length"]
        self.length = row["A"]+row["T"]+row["C"]+row["G"]
        
    def update_func(self,new_func_name):
        self.func_name = new_func_name
        row = kmer_tools.search_kmer.read_kdat_dimerRow(k_data_folder,self.domain,func_name=new_func_name,orga_name = self.orga_name)
        self.row = row
        self.length = row["A"]+row["T"]+row["C"]+row["G"]

    def get_count(self,word):
        return self.row[word]
    
    def _get_total_count(self,word):
        return self._total_row[word]

    def normalize_count(self,count,err=None):
        if type(count) == list:
            return [self.normalize_count(entry,err) for entry,err in zip(count,err)]
        else:
            count=np.array(count)
            length=self.length
            num_tot=4**2
            if err is None:
                #print("monomers, count, Expactation:",count[0],self.dimer_words,count[0]/(length-1)*num_tot)
                return count/(length-1)*num_tot
            else:
                return count/(length-1)*num_tot, err/(length-1)*num_tot

    def normalize_count_arr(self,count,err=None):
        count=np.array(count)
        length=self.length
        num_tot=4**2
        if err is None:
            #print("monomers, Expactation:",self.dimer_words,count[0]/(length-1)*num_tot)
            return count/(length-1)*num_tot
        else:
            return count/(length-1)*num_tot, np.array(err)/(length-1)*num_tot

    def get_count_normalized(self,word):
        length=self.length
        num_kmer=16
        return self.get_count(word)/(length-1)*num_kmer

    def get_mono_expectation(self,word):
        length=self.length
        if(length==0): 
            print("warning: length is zero")
            return 0
        monomers = [word[0],word[1]]
        p_x = [self.row[letter]/length for letter in monomers]
        p_xy = np.prod(p_x)
        #print(monomers,p_x,p_xy)
        E_xy = p_xy * length
        sigma2 = length * p_xy * (1- p_xy)
        err=np.sqrt(sigma2)
        #print("monomers, Expactation:",monomers,E_xy)
        return E_xy,err
    
    def get_total_expectation(self,word):   #get an estimate for error?
        total_count = self._get_total_count(word)
        total_len   = self.chromo_length
        area_len    = self.length
        if area_len==0: return 0
        return total_count/total_len * area_len,0
    
    def get_dist(self,func=None):
        array = None
        if func in ["count", None]:
            array =  np.array([(self.get_count(word),0)                for word in dimer_dist.dimer_words])
        elif func == "normalized":
            array =  np.array([(self.get_count_normalized(word),0)     for word in dimer_dist.dimer_words])
        elif func == "mono_expectation":
            array =  np.array([self.get_mono_expectation(word)     for word in dimer_dist.dimer_words])
        elif func == "total_expectation":
            array =  np.array([self.get_total_expectation(word)    for word in dimer_dist.dimer_words])
        if array is not None:
            array = np.transpose(array)
            #print(func)
            #print(array[0],array[1])
            return array[0],array[1]
        else:
            print ("returning ditribution of " + func.__name__)
            return [func(word) for word in dimer_dist.dimer_words]
    
    def get_func_diff(self,func1,func2):
        vals1,err1= self.normalize_count(self.get_dist(func1))
        vals2,err2= self.normalize_count(self.get_dist(func2))
        return vals1-vals2,np.sqrt(err1**2+err2*2)

    def get_dist_diff(self,dist1,dist2):
        vals1,err1= dist1.normalize_count(self.get_dist(dist1))
        vals2,err2= dist2.normalize_count(self.get_dist(dist2))
        return vals1-vals2,np.sqrt(err1**2+err2*2)


def dimer_bias_diff(orga_name,func_name1,func_name2,domain=None,NaN=False):
    dist = dimer_dist(orga_name,func_name1,domain=domain)
    if (NaN and np.sum(np.abs(dist.get_dist("count")[0]))==0):
        return np.nan,np.nan
    vals1,errs1 = dist.get_func_diff("count","mono_expectation")
    dist.update_func(func_name2)
    if (NaN and np.sum(np.abs(dist.get_dist("count")[0]))==0):
        return np.nan,np.nan
    vals2,errs2 = dist.get_func_diff("count","mono_expectation")
    return vals1-vals2, np.sqrt(errs1**2 + errs2**2)

def dimer_bias_pearson(orga_name,func_name1,func_name2,domain=None,NaN=False):
    dist = dimer_dist(orga_name,func_name1,domain=domain)
    if (NaN and np.sum(np.abs(dist.get_dist("count")[0]))==0):
        return np.nan,np.nan
    vals1,errs1 = dist.get_func_diff("count","mono_expectation")
    dist.update_func(func_name2)
    if (NaN and np.sum(np.abs(dist.get_dist("count")[0]))==0):
        return np.nan,np.nan
    vals2,errs2 = dist.get_func_diff("count","mono_expectation")
    r, p = stat.pearsonr(vals1,vals2)
    return r

def dimer_bias_dist_cum(orga_name,func_name1,func_name2,domain=None):
    diff_vals, diff_errs = dimer_bias_diff(orga_name,func_name1,func_name2,domain=domain,NaN=True)
    return np.sum(np.abs(diff_vals)), np.sqrt(np.sum(diff_errs**2))


organism="Homo sapiens"
domain="animalia"
dist = dimer_dist(organism,"genes",domain)
#print([(x,y,z,w) for x,y,z,w in zip(dist.dimer_words,dist.get_dist("count")[0],dist.get_dist("mono_expectation")[0],dist.get_func_diff("count","mono_expectation")[0])])