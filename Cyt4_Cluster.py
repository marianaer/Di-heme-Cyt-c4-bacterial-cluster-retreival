#!/usr/bin/env python
# coding: utf-8

# Importing libraries
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import re
import subprocess
import os
from path import Path
from collections import defaultdict


# Signing in to NCBI
Entrez.email = "mescobar@gmail.com" 

gen="Acidithiobacillus sp."
query="AVK42812"
protquery=gen+"_"+query


# Running BlastP. Excluiding Archaeobacteria, metagenomes and uncultured samples in Non Redundant database.
blast_handle = NCBIWWW.qblast('blastp', 'nr', sequence=query, entrez_query='all [filter] NOT(environmental samples[organism] OR metagenomes[orgn] OR txid2157[organism] OR txid119977[organism])', hitlist_size=100)


# Opening Blats results
with open("blast_"+protquery+".xml", "w") as out_handle:
    out_handle.write(blast_handle.read())
blast_handle.close()

blast_handle = open("blast_"+protquery+".xml")

blast_record = NCBIXML.read(blast_handle)


# ITerating through blast file
res=[]
for record in NCBIXML.parse(open("blast_"+protquery+".xml")):
    for align in record.alignments:
             for hsp in align.hsps:
                a=(((align.title).split('['))[1])
                res.append({'Taxon': (a.split(']'))[0] ,'Hit ID':  align.accession, 'E-value': hsp.expect, 'Length':hsp.align_length ,'Identity':(hsp.identities/hsp.align_length*100), 'Query':record.query_id, 'Hit Protein name':re.sub(r'\s\[.*', '',align.title.split("|")[2])})

res[0]['Hit Protein name']="Cyc 1"    
df=pd.DataFrame(res)
df['Hit ID'].to_csv('IDblast.txt',  sep=' ', index=False, header=False)


# Shell script for fetching adjacent proteins:
subprocess.call("./flankinggenes.sh")

# Importing adjacent proteins
flank=[]
seen=[]

with open("FlankingProteins.txt") as flanks:
    for line in flanks:
        line=line.rstrip("\n")
        line=re.sub(r'\.\d', '', line)
        line=line.split('\t')  
        if line not in seen:
            seen.append(line)
            flank.append(line)


# Ordeiring proteins with their flanking proteins
def adjacent_proteins(protein):
    currentProt=[]
    for prot in flank:
        if prot[3]==protein:
            currentProt.append(prot)
    return(currentProt)


prots=[]
passed=[]

# getting rid of proteins with no  id
for prot in flank:
    if prot[0]=='--':
        flank.remove(prot)
# Sorting
for prot in flank:
    a=adjacent_proteins(prot[3])
    while len(a)<5:
        b=[a[0][0], "-", "-", a[0][3]]
        a.insert(0, b)
      
    for i, prot in enumerate(a):
        if a[i][3]==a[i][1]:
            temp=a[i]
            a[i]=a[2]
            a[2]=temp
            Left_2=a[2-2]
            Left_1=a[2-1]
            middle=a[2],i
            Right_1=a[2+1]
            Right_2=a[2+2]
    prots.append({'Hit ID': a[i][3] , '-2 Protein name': Left_2[2], '-2 Protein ID':Left_2[1],'-2 Conserved Domain':"-",'-1 Protein name': Left_1[2], '-1 Protein ID':Left_1[1],'-1 Conserved Domain':"-",'+1 Protein name': Right_1[2], '+1 Protein ID':Right_1[1], '+1 Conserved Domain':"-",'+2 Protein name': Right_2[2], '+2 Protein ID': Right_2[1], '+2 Conserved Domain':"-", 'Hit Protein Conserved Domain':"-"})
    
seen = set()
new_prots = []
for d in prots:
    t = tuple(d.items())
    if t not in seen:
        seen.add(t)
        new_prots.append(d)

pd.DataFrame(new_prots)


# Merging both the adjacent protein dict with our initial blast result dict

d = defaultdict(dict)
for l in (new_prots, res):
    for elem in l:
        d[elem['Hit ID']].update(elem)
ProteinList = d.values()
P=pd.DataFrame(list(ProteinList))
P = P.fillna('-')


# Writing IDs of the proteins in order to find their CDs
with open("CDDSearch.txt", "w") as out_handle:
    for prot in flank:
        line=prot[1]
        out_handle.write(line+"\n")
   



# Shell script for fetching adjacent proteins through perl script:
subprocess.run("./bwrpsb.pl < CDDSearch.txt > hitdata.txt", shell=True)


# Parsing the CDD Batch file
with open('hitdata.txt', 'r') as file:
    data = file.read().splitlines(True)
with open('ConservedDomains.txt', 'w') as out:
    out.writelines(data[11:])

cdd=[]
with open('ConservedDomains.txt', 'r') as file:

    for line in file:
        line=line.split()
        if len(line)>0: # Ignore blank lines
            cdd.append({'Hit ID': line[2], 'CD':line[10]})


# Adding each protein its conserved domain
c=[]
p3=list(ProteinList)
for i in p3:
    if '-2 Protein ID' and i['Hit ID'] not in c:
        c.append(i['Hit ID'])
        
cont=0
for i in range(len(p3)):
    if '-1 Protein ID' in p3[i]:
        cont=cont+1

# Assigning to each protein its respective CD based on th ID:
for count in range(cont):
    for j in cdd:
        for i in c:
            if p3[count]['Hit ID']==i:
                if p3[count]['-2 Protein ID']==j['Hit ID']: 
                    p3[count]['-2 Conserved Domain']=j['CD']
                elif p3[count]['-1 Protein ID']==j['Hit ID']:
                    p3[count]['-1 Conserved Domain']=j['CD'] 
                elif p3[count]['+2 Protein ID']==j['Hit ID']:
                    p3[count]['+2 Conserved Domain']=j['CD']
                elif p3[count]['+1 Protein ID']==j['Hit ID']:
                    p3[count]['+1 Conserved Domain']=j['CD'] 
                elif p3[count]['Hit ID']==j['Hit ID']: 
                    p3[count]['Hit Protein Conserved Domain']=j['CD'] 
                    
# Column sorting:
results=pd.DataFrame(p3)
columnsTitles = ['Query', 'Hit ID', 'Hit Protein name', 'Taxon', 'Length', 'Identity', 'E-value', '-2 Protein name', '-2 Conserved Domain', '-2 Protein ID', '-1 Protein name', '-1 Conserved Domain', '-1 Protein ID', 'Hit Protein name', 'Hit Protein Conserved Domain', '+1 Protein name', '+1 Conserved Domain', '+1 Protein ID', '+2 Protein name', '+2 Conserved Domain', '+2 Protein ID' ]
results = results.reindex(columns=columnsTitles)

# Save to a file
results.to_csv(protquery+'_Results.csv',  sep='\t', index=False, header=True)

# Delete all .txt files created during the process
direc = Path("/home/mescobar/Escritorio/Mauro/")
files = direc.walkfiles('*.txt')
files2 = direc.walkfiles('*.xml')

for file in files:
    file.remove()
for file in files2:
    file.remove()   

