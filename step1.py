#!/usr/bin/python3
import os
import sys
import subprocess
import pandas as pd
import re

print('Hello, welcome!\n')
print('This program is designed for users to process a series of sequences related to user-defined protein families and taxonomic group. After screening these sequences by sequence number and sequence length for the user, applying clustalo to perform multiple sequence alignments. Then plot the level of protein sequence conservation across the species within that taxonomic group. Then scan the selected protein sequences of interest with motifs from the PROSITE database. Finally, the hydropathy plot of and calculate statistics of protein properties of selected sequences are offered in the program. \n')

os.system('cp /localdisk/data/BPSM/ICA2/pullseq .')
#step1: Filter the number of sequences for the user and ask the user if they want to proceed to the next step
protein1 = input("Please input your query protein name: \n")
organism = input("Please input your organism: \n")
#ask user to input the query protein and organism
subprocess.call(f"esearch -db protein -query '{protein1}[protein] AND {organism}[organism]' | efetch -format fasta > protein.txt" , shell=True)
#Use esearch then efetch to get the fasta seqs

pro_count = open("protein.txt").read().count('>')
#count how many sequences are in protein.txt
f = open("protein.txt")
org = []
for line in f:
  if line.startswith('>'):
    m = line.split("[")[1]
    n = m.split("]")[0]
    org.append(n)
    set(org)
    
#first create a new list containing species names as org

org_num = len(set(org))
if pro_count > 1000:
  print('WARNING: The number of sequences is too large to proceed to the next step, please exit or re-enter the protein and organism\n')
if pro_count <= 1000:
  if org_num > 1:
      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from " + str(org_num) + " different species\n")
  if org_num <= 1:
      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from the same specie\n")
#Filter the number of sequences for the user and ask the user if they want to proceed to the next step


#step2: this step contain the assignment step2 and some part pf step4
#1:Use sequence length as the criterion to filter sequences for users to determine
#2:do the clustalo
#3:do the plot
#4:use pepstats to calculate statistics of protein properties 

f = open("protein.txt")
acc = []
for line in f:
  if line.startswith('>'):
    a = line.split(" ")[0]
    acc.append(a)
   
#Extract all accession numbers as a list
f = open("protein.txt")
spe = []
for line in f:
  if line.startswith('>'):
    b = line.split("[")[1]
    c = b.split("]")[0]
    spe.append(c)
  
#Extract all species name as a list
f = open("protein.txt")
seqtotal = []
seq = []
for line in f:
    if line.startswith('>'):
        seqtm=''.join(seq)
        seqtotal.append(seqtm)
        seq=[]
        d = line.split("]")[1]
    else:
        e = line
        seq.append(e)
#Extract all sequences as a list. 
seqtm=''.join(seq)
seqtotal.append(seqtm)
seqtotal.pop(0)

#since the first element is 0, so use 'pop' to delete the 0 element and add the last sequence in the seqtotal list 

length = []
for i in seqtotal:
  g = len(i)
  length.append(g)
#count the sequences length and make a list 


s1 = pd.Series(acc)
s2 = pd.Series(spe)
s3 = pd.Series(seqtotal)
s4 = pd.Series(length)
overview = pd.DataFrame( { 'Accession number' : s1, 'species' : s2, 'sequence' : s3, 'length' : s4} )
pd.set_option('display.max_rows', None)
#make a data frame of accession number, species, sequences, length of sequences. And show this to the users
print(overview)

def minlength():
  length = input('What is the minimum length of the sequence you want to filter \n')
  length = int(length)
  os.system(f'./pullseq -i protein.txt -m {length} > selectseq.txt') 
  return min
#write a new function: select the minimum limit of sequences and save as a new file:selectseq.txt
def maxlength():
  length = input('What is the maximum length of the sequence you want to filter \n')
  length = int(length)
  os.system(f'./pullseq -i protein.txt -a {length} > selectseq.txt') 
  return max
#write a new function: select the maximum limit of sequences and save as a new file:selectseq.txt

y1 = input('Now you can see a dataframe of sequence length on your screen. Do you want to filter for sequence length? please input yes/no: \n')
# Users can choose whether to filter sequences by sequence length.
if y1 == 'no':
  x = 'clustalo -i protein.txt -o alignment.txt'
  os.system(x)
  i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
  os.system(i)
  os.system(f'pepstats -sequence protein.txt -outfile final.pepstats')
#If the user chooses not to filter the sequence, the programme will apply protein.txt to do the multiple alignment and plot the alignment result
if y1 == 'yes':
  y2 = input('How do you want to limit the length of the sequence? Please input: max/min/max and min: \n')
  if y2 == 'max':
    maxlength()
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
    os.system(f'pepstats -sequence selectseq.txt -outfile final.pepstats')
  if y2 == 'min':
    minlength()
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
    os.system(f'pepstats -sequence selectseq.txt -outfile final.pepstats')
  if y2 == 'max and min':
    length1 = input('What is the maximum length of the sequence you want to filter(Please enter Arabic numerals): \n')
    length2 = input('What is the minimum length of the sequence you want to filter(Please enter Arabic numerals): \n')
    length1 = int(length1)
    length2 = int(length2)
    os.system(f'./pullseq -i protein.txt -m {length2} -a {length1} > selectseq.txt') 
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
    os.system(f'pepstats -sequence selectseq.txt -outfile final.pepstats')
#If the user wants to limit the length of the sequence, the user has three choices. Users can choose to limit the maximum or minimum length of the sequence, or they can choose to limit the maximum and minimum length of the sequence at the same time. For example, the user can select sequences whose length is less than x and generate all these sequences as a new file  

print('please see the alignment result and plotcon pdf in your folder. \n')


#step3
os.system('mkdir motif')
os.system('mkdir hydropathy_plot')
os.system('touch temp.fasta')
pp = {}
if y1 == 'no':
  f = open('protein.txt')
  count = 0
  for i in acc:
    pp[i] = ''
    a = seqtotal[count]
    count += 1
    pp[i] = a


  for acc1,seq120 in pp.items() :
    print(f"{acc1} {seq120[:5]}")
    u = acc1.replace('>','')
    file = open(f'temp.fasta','w')
    file.write(f'{acc1}\n')
    file.write(f'{seq120}')
    file.close()
    os.system(f'pepwindow -sequence temp.fasta -window 2 -graph pdf -goutfile ./hydropathy_plot/{u}')
    os.system(f'patmatmotifs -sequence temp.fasta -outfile ./motif/{u}.patmatmotifs')
    with open(f'./motif/{u}.patmatmotifs') as motifs :
      for eachline in motifs:
        motifs2 = ''
        if re.search(r'Motif = (.*)',eachline):
          motifs2 = re.search(r'Motif = (.*)',eachline).group(1)
          with open(f'motifsfinal.txt','a') as mymotif:
            mymotif.write(f'{u}: ' + f'the name of motif : {motifs2}\n')


if y1 == 'yes':
  fy = open('selectseq.txt')
  accy = []
  for line in fy:
    if line.startswith('>'):
      a = line.split(" ")[0]
      accy.append(a)

  fy = open('selectseq.txt')
  seqtotaly = []
  seqy = []
  for line in fy:
    if line.startswith('>'):
      seqtmy=''.join(seqy)
      seqtotaly.append(seqtmy)
      seqy=[]
      d = line.split("]")[1]
    else:
      e = line
      seqy.append(e)
  seqtmy=''.join(seqy)
  seqtotaly.append(seqtmy)
  seqtotaly.pop(0)

  count = 0 
  for i in accy:
    pp[i] = ''
    a = seqtotaly[count]
    count += 1
    pp[i] = a

  for acc1,seq120 in pp.items() :
    u = acc1.replace('>','')
    file = open(f'temp.fasta','w')
    file.write(f'{acc1}\n')
    file.write(f'{seq120}')
    file.close()
    os.system(f'pepwindow -sequence temp.fasta -window 2 -graph pdf -goutfile ./hydropathy_plot/{u}')
    os.system(f'patmatmotifs -sequence temp.fasta -outfile ./motif/{u}.patmatmotifs')
    with open(f'./motif/{u}.patmatmotifs') as motifs :
      for eachline in motifs:
        motifs2 = ''
        if re.search(r'Motif = (.*)',eachline):
          motifs2 = re.search(r'Motif = (.*)',eachline).group(1)
          with open(f'motifsfinal.txt','a') as mymotif:
            mymotif.write(f'{u}: ' + f'the name of motif : {motifs2}\n')

print('please see the motifs detail in motif folder\n An overview file of the motif name is available. Please check motifsfinal.txt in your folder\n')









    
      




