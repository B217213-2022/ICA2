#!/usr/bin/python3
import os
import sys
import subprocess
import pandas as pd
import re

#step1 
#protein1 = input("what is your query protein: ")
#organism = input("what is your organism: ")
#ask user to input the query protein and organism
#subprocess.call(f"esearch -db protein -query '{protein1}[protein] AND {organism}[organism]' | efetch -format fasta > protein.txt" , shell=True)
#Use esearch then efetch to get the fasta seqs

#pro_count = open("protein.txt").read().count('>')
#count how many sequences are in protein.txt
#f = open("protein.txt")
#org = []
#for line in f:
#  if line.startswith('>'):
#    m = line.split("[")[1]
#    n = m.split("]")[0]
#    org.append(n)
#    set(org)
#    print(org)
#first create a new 

#org_num = len(set(org))
#if pro_count > 1000:
#  print('The number of sequences is too large to proceed to the next step, please exit or re-enter the protein and organism')
#else:
#  if org_num > 1:
#      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from " + str(org_num) + " different species")
#  else:
#      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from the same specie")
#
#step2
f = open("protein.txt")
acc = []
for line in f:
  if line.startswith('>'):
    a = line.split(" ")[0]
    acc.append(a)
    print(acc)
#Extract all accession numbers as a list
f = open("protein.txt")
spe = []
for line in f:
  if line.startswith('>'):
    b = line.split("[")[1]
    c = b.split("]")[0]
    spe.append(c)
    print(spe)
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
print(seqtotal)
#since the first element is 0, so use 'pop' to delete the 0 element and add the last sequence in the seqtotal list 

length = []
for i in seqtotal:
  g = len(i)
  length.append(g)
#count the sequences length and make a list 
print(length)

s1 = pd.Series(acc)
s2 = pd.Series(spe)
s3 = pd.Series(seqtotal)
s4 = pd.Series(length)
pd.DataFrame( { 'Accession number' : s1, 'species' : s2, 'sequence' : s3, 'length' : s4} )
#make a data frame of accession number, species, sequences, length of sequences.

def minlength():
  length = input('What is the minimum length of the sequence you want to filter \n')
  length = int(length)
  os.system(f'./pullseq -i protein.txt -m {length} > selectseq.txt') 
  return min
#write a new function: select the minimum limit of sequences and save as a new file:selectseq.txt
def maxlength():
  length = input('What is the maximum length of the sequence you want to filter \n')
  length = int(length)
  os.system(f'./pullseq -i protein.txt -m {length} > selectseq.txt') 
  return max
#write a new function: select the maximum limit of sequences and save as a new file:selectseq.txt

y1 = input('Now you can see a dataframe of sequence length in your folder. Do you want to filter for sequence length? please input yes/no: \n')
if y1 == 'no':
  x = 'clustalo -i protein.txt -o alignment.txt'
  os.system(x)
  i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
  os.system(i)
else:
  y2 = input('what is the limitaion you want to choose for sequence length? Please input: max/min/max and min: \n')
  if y2 == 'max':
    maxlength()
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
  if y2 == 'min':
    minlength()
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
  if y2 == 'max and min':
    length1 = input('What is the maximum length of the sequence you want to filter \n')
    length2 = input('What is the minimum length of the sequence you want to filter \n')
    length1 = int(length1)
    length2 = int(length2)
    os.system(f'./pullseq -i protein.txt -m {length2} -a {length1}> selectseq.txt') 
    x = 'clustalo -i selectseq.txt -o alignment.txt'
    os.system(x)
    i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
    os.system(i)
        

print('please see the alignment result and plot pdf in your folder. \n')


#step3
pp = {}
if y1 == 'no':
  f = open('protein.txt')
  count = 0 
  for i in acc:
    pp[i] = ''
    a = seqtotal[count]
    count += 1
    pp[i] = a

    

  os.system('mkdir motif')

  for acc1,seq120 in pp.items() :
    u = acc1.replace('>','')
    file = open(f'temp.fasta','w')
    file.write(f'{acc1}\n')
    file.write(f'{seq120}')
    file.close
      #os.system(f'cd ./fasta')
    os.system(f'patmatmotifs -sequence temp.fasta -outfile ./motif/{u}.patmatmotifs')
    with open(f'./motif/{u}.patmatmotifs') as motifs :
      for eachline in motifs:
        motifs2 = ''
        if re.search(r'Motif = (.*)',eachline):
          motifs2 = re.search(r'Motif = (.*)',eachline).group(1)
          with open(f'motifsfinal.txt','a') as mymotif:
            mymotif.write(f'{u}: ' + f'the name of motif : {motifs2}\n')

if y1 == 'yes':
  f = open('selectseq.txt')
  count = 0 
  for i in acc:
    pp[i] = ''
    a = seqtotal[count]
    count += 1
    pp[i] = a

    

  os.system('mkdir motif')

  for acc1,seq120 in pp.items() :
    u = acc1.replace('>','')
    file = open(f'temp.fasta','w')
    file.write(f'{acc1}\n')
    file.write(f'{seq120}')
    file.close
      #os.system(f'cd ./fasta')
    os.system(f'patmatmotifs -sequence temp.fasta -outfile ./motif/{u}.patmatmotifs')
    with open(f'./motif/{u}.patmatmotifs') as motifs :
      for eachline in motifs:
        motifs2 = ''
        if re.search(r'Motif = (.*)',eachline):
          motifs2 = re.search(r'Motif = (.*)',eachline).group(1)
          with open(f'motifsfinal.txt','a') as mymotif:
            mymotif.write(f'{u}: ' + f'the name of motif : {motifs2}\n')








range(len(acc))

    
      

#x = 'clustalo -i protein.txt -o alignment.txt'
#os.system(x)

#i = 'plotcon -sequences alignment.txt -winsize 4 -graph pdf'
#os.system(i)

#p = 'patmatmotifs alignment.txt -rformat listfile'
#os.system(p)

#for i in seqtotal:
#  p = 'patmatmotifs {i} -rformat listfile'
#  os.system(p)




