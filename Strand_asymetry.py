### Script to calculate strand asymmetry associated with transcription.
### Requirements:
### genes to be translated to same strand i.e, mutations occuring in genes on minus strand are reverse transcribed.
### Input file contains:
#### Haplotype predicted expression (high/low)
#### ref allele (strand that is transcribed)
#### alt allele (non-transcribed strand)
#### median expression value from GTEx data or GEUVADIS
#### Number of tissues associated with cis acting eQTL. 
### usage python Strand_asymmetry.py input_file expression_threshold(float) minimun_tissue(int)

import sys, re
from collections import Counter

parameters = sys.argv
max_exp = parameters[2]
min_tissues = parameters[3]


outfile = parameters[1].split('.')[0] + ".spec"

myfile = open(parameters[1]) # open VCF change to sys parameter to stream line
raw = myfile.readlines()
myfile.close()

## loop over highly expressed haplotypes

record = []

for line in range(0,len(raw)): # loop over 
        if raw[line].strip("\n").split(" ")[6] == 'high' and len(raw[line].strip("\n").split(" ")[4]) == 1 and len(raw[line].strip("\n").split(" ")[5]) == 1 and raw[line].strip("\n").split(" ")[8] > max_exp and raw[line].strip("\n").split(" ")[9] > min_tissues : # criteria [high expression and SNV mutation with filter on expression level and number f cis-eqtls 
                alt = raw[line].strip("\n").split(" ")[5].replace("\n","").strip()
                ref = raw[line].strip("\n").split(" ")[4].strip()
                record.append(ref + alt + "\n")
trin = open("%s.high"%(outfile), 'w')       # write out
trin.writelines(record)
trin.close()


record_occur = []
# open file of SNV and read as string
fileName=open(outfile+".high")
lines = [i for i in fileName.readlines()]

#print(str(lines).replace("\n",""))
#Count string occurences 
record_occur = str(Counter(lines))
record_oc = Counter(record)

tot_mut = sum(int(x) for x in re.findall(r'[0-9]+', record_occur.strip("\n")))

outp = []
outp.append("Subst" + "\t" + "Count" + "\t" + "%_total" + "\t" + "r_Subst" + "\t" + "r_Count" + "\t" + "r_%_total" + "\n")

outp.append(\
'CT' + '\t' + str(record_oc.get('CT\n')) + '\t' + str(record_oc.get('CT\n')/float(tot_mut)) + '\t' + 'GA' + '\t' + str(record_oc.get('GA\n')) + '\t' + str(record_oc.get('GA\n')/float(tot_mut)) + '\n' + 'CG' + '\t' + str(record_oc.get('CG\n')) + '\t' + str(record_oc.get('CG\n')/float(tot_mut)) + '\t' + 'GC' + '\t' + str(record_oc.get('GC\n')) + '\t' + str(record_oc.get('GC\n')/float(tot_mut)) + '\n' + 'CA' + '\t' + str(record_oc.get('CA\n')) + '\t' + str(record_oc.get('CA\n')/float(tot_mut)) + '\t' + 'GT' + '\t' + str(record_oc.get('GT\n')) + '\t' + str(record_oc.get('GT\n')/float(tot_mut)) + '\n' + 'AT' + '\t' + str(record_oc.get('AT\n')) + '\t' + str(record_oc.get('AT\n')/float(tot_mut)) + '\t' + 'TA' + '\t' + str(record_oc.get('TA\n')) + '\t' + str(record_oc.get('TA\n')/float(tot_mut)) + '\n' + 'AG' + '\t' + str(record_oc.get('AG\n')) + '\t' + str(record_oc.get('AG\n')/float(tot_mut)) + '\t' + 'TC' + '\t' + str(record_oc.get('TC\n')) + '\t' + str(record_oc.get('TC\n')/float(tot_mut)) + '\n' + 'AC' + '\t' + str(record_oc.get('AC\n')) + '\t' + str(record_oc.get('AC\n')/float(tot_mut)) + '\t' + 'TG' + '\t' + str(record_oc.get('TG\n')) + '\t' + str(record_oc.get('TG\n')/float(tot_mut)) + '\n' + str(tot_mut))

trin = open("mut_counts_high.out", 'w')       # write out
trin.writelines(outp)
trin.close()

## loop through low expression haplotypes

record1 = []
for line in range(0,len(raw)): # loop over 
        if raw[line].strip("\n").split(" ")[6] == 'low' and len(raw[line].strip("\n").split(" ")[4]) == 1 and len(raw[line].strip("\n").split(" ")[5]) == 1 and raw[line].strip("\n").split(" ")[8] > max_exp and raw[line].strip("\n").split(" ")[9] >= min_tissues: # criteria [high expression and SNV mutation with filter on expression level and number f cis-eqtls 
                alt1 = raw[line].strip("\n").split(" ")[5]
                ref1 = raw[line].strip("\n").split(" ")[4]
                record1.append(ref1 + alt1 + "\n")
trin = open("%s.low"%(outfile), 'w')       # write out
trin.writelines(record1)
trin.close()

record_occur1 = []

fileName=open(outfile+".low")
lines = [i for i in fileName.readlines()]

record_occur1 = str(Counter(lines)).strip()
record_oc1 = Counter(lines)
tot_mut1 = sum(int(x) for x in re.findall(r'[0-9]+', record_occur1.strip("\n")))


outp1 = []
outp1.append("Subst" + "\t" + "Count" + "\t" + "%_total" + "\t" + "r_Subst" + "\t" + "r_Count" + "\t" + "r_%_total" + "\n")

outp1.append('CT' + '\t' + str(record_oc1.get('CT\n')) + '\t' + str(record_oc1.get('CT\n')/float(tot_mut1)) + '\t' + 'GA' + '\t' + str(record_oc1.get('GA\n')) + '\t' + str(record_oc1.get('GA\n')/float(tot_mut1)) + '\n' + \
\
'CG' + '\t' + str(record_oc1.get('CG\n')) + '\t' + str(record_oc1.get('CG\n')/float(tot_mut1)) + '\t' + 'GC' + '\t' + str(record_oc1.get('GC\n')) + '\t' + str(record_oc1.get('GC\n')/float(tot_mut1)) + '\n' + \
\
'CA' + '\t' + str(record_oc1.get('CA\n')) + '\t' + str(record_oc1.get('CA\n')/float(tot_mut1)) + '\t' + 'GT' + '\t' + str(record_oc1.get('GT\n')) + '\t' + str(record_oc1.get('GT\n')/float(tot_mut1)) + '\n' + \
\
'AT' + '\t' + str(record_oc1.get('AT\n')) + '\t' + str(record_oc1.get('AT\n')/float(tot_mut1)) + '\t' + 'TA' + '\t' + str(record_oc1.get('TA\n')) + '\t' + str(record_oc1.get('TA\n')/float(tot_mut1)) + '\n' + \
\
'AG' + '\t' + str(record_oc1.get('AG\n')) + '\t' + str(record_oc1.get('AG\n')/float(tot_mut1)) + '\t' + 'TC' + '\t' + str(record_oc1.get('TC\n')) + '\t' + str(record_oc1.get('TC\n')/float(tot_mut1)) + '\n' + \
\
'AC' + '\t' + str(record_oc1.get('AC\n')) + '\t' + str(record_oc1.get('AC\n')/float(tot_mut1)) + '\t' + 'TG' + '\t' + str(record_oc1.get('TG\n')) + '\t' + str(record_oc1.get('TG\n')/float(tot_mut1)) + '\n' + str(tot_mut1))

trin = open("mut_counts_low.out", 'w')       # write out
trin.writelines(outp1)
trin.close()
