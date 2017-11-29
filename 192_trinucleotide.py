### Initial script to count 192 spectrum from VCF file using reference genome fasta file
### No reference context is provided i.e transcribed strand or direction of replication
### Output needs to be tidied and procedure streamlines (unnecessary intermediate files)
### usage python 192_trinucleotide_freq.py file.vcf reference_genome.fa
### output file is file_trinucleotide_frequencies.out


import requests, sys, time, re
import StringIO, pysam
from collections import Counter

parameters = sys.argv
VCF = parameters[1]
fasta = parameters[2]

outfile = parameters[1].split('.')[0]


sys.stdout = StringIO.StringIO()
# For writing multiple lines to a file 
class cfile(file):
    #subclass file to have a more convienient use of writeline
    def __init__(self, name, mode = 'r'):
        self = file.__init__(self, name, mode)

    def wl(self, string):
        self.writelines(string + '\n')
        return None

## take vcf file and pull 5' and 3' nucleotides. vcf mapped to GRCH37
myfile = open(VCF) # open VCF change to sys parameter to stream line
raw = myfile.readlines()
myfile.close()

genome = pysam.Fastafile(fasta)

req = [] 
for line in range(0,len(raw)): # loop over vcf file
       chr_pos = raw[line].strip("\n").split("\t")[1] # take coordinates out
       chr_no = raw[line].strip("\n").split("\t")[0]   # take chromosome position
       req.append(chr_no + "\t" + chr_pos + "\n")      # create variable containing number and pos
       # search reference genome based on req.append
##      genome = pysam.Fastafile('/data4/Declan/Internship/GRCH37/Homo_sapiens.GRCh37.75.dna.toplevel.fa')
       sequence = genome.fetch(chr_no, int(chr_pos)-2, int(chr_pos) + 1)
###     print(sequence)
       r_mut = list(str(sequence))     # take triplet as list
       r_mut[1] = raw[line].strip("\n").split("\t")[4]         #replace middle value with alternative allele from vcf file
       new_string = ''.join(r_mut)     #join 
       fid = cfile("%s_spectr_pos.out"%(outfile), 'a')
       out = str(sequence) + new_string        # write out
       fid.wl(out +"\t" + chr_no + "\t" + chr_pos)
       fid.close()
f = open("%s_spectr_pos.out"%(outfile),"r")
lines = f.readlines()
f.close()
f = open("%s_tidy_data.out"%(outfile),"w")
for line in lines: # remove lines with warning (can prob be removed), refence N's, subs str will be six char long and that the reference nuc != alt allele. 
        if not line.startswith('Y') and not line.startswith('N') and line.strip("\n").split("\t")[0][2] != 'N' and len(line.strip("\n").split("\t")[0])== 6 and line.strip("\n").split("\t")[0][1] != line.strip("\n").split("\t")[0][4] and line.strip("\n").split("\t")[0][4].isalpha():
                f.write(line)
                print(line)

f.close()
#
#
#### Counts

myfile = open("%s_tidy_data.out"%(outfile)) # open VCF change to sys parameter to stream line
raw = myfile.readlines()
myfile.close()

per_row = []
for line in raw:
    per_row.append(line.strip().split('\t'))

per_column = zip(*per_row)

print(per_column[0])

counts = Counter(per_column[0])
dic_out = []
tri_nuc = []
counts.keys()
for key, value in counts.items():
        temp = key + str(value)
        dic_out.append(temp)
counts_trinucleotide = str(dic_out).replace(",","\n").replace("[","").replace("]","").replace("'","").replace(" ", "")
#
tri = open("%s_untidy_tridata.out"%(outfile), 'w')      # write out
tri.writelines(counts_trinucleotide)
tri.close()

### tidy and calculate frequencies

myfile = open("%s_untidy_tridata.out"%(outfile)) # open VCF change to sys parameter to stream line
raw1 = myfile.readlines()
myfile.close()

counts1 = 0

def my_split(s):        # SPLIT AT FIRST NUMBER quick but creates messy output
    return filter(None, re.split(r'(\d+)', s))

for c in range(0,len(raw1)):
        temp0 = raw1[c].strip("\t").split("\n")[0]
        xy = (my_split(temp0))
        xx = str(xy).replace("[","").replace("]","").replace("'","").replace(" ","")
        xx_split = str(xx).split(",")
        temp_count = xx_split[1]
        tri_nucl = str(tri_nuc).replace(",","\n").replace("'","").replace("[","").replace("]","").replace(" ", "").replace("Q", "\t")
        counts1 = counts1 + float(str(xx_split[1]))

#### BELOW IS VER MESSY DICTIONARY OUTPUT DIDNT HAVE NICE FORMATTING
for k in range(0,len(raw1)):
        temp0 = raw1[k].strip("\t").split("\n")[0]
        xy = (my_split(temp0))
        xx = str(xy).replace("[","").replace("]","").replace("'","").replace(" ","")
        xx_split = str(xx).split(",")
        temp2 = str(xx_split[0]) + "Q" + str(xx_split[1]) + "Q" + str(float(xx_split[1])/float(counts1))
        tri_nuc.append(temp2)
        tri_nucl = str(tri_nuc).replace("'","").replace("[","").replace("]","").replace(" ", "").replace("Q", "\t").replace(",","\n")



#print(tri_nucl)
trin = open("%s_trinucleotide_frequencies.out"%(outfile), 'w')  # write out ref triplet + alternative triplet, counts + frequency
trin.writelines(tri_nucl)
trin.close()
