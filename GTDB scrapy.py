#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scrapy
from urllib.request import urlopen
from bs4 import BeautifulSoup
import re
import pandas as pd


# In[2]:


def completeness(taxon):
    html = urlopen(taxon)
    soup =BeautifulSoup(html, 'lxml')
    text = soup.get_text()
    a = soup.find_all("a")

    cleantext = BeautifulSoup(str(a), "lxml").get_text()
    gid=cleantext.split()
    gid=str(gid[34])
    gid=gid.replace(",", "")
    gid=gid.replace("]", "")

    url = ("https://gtdb.ecogenomic.org/genomes?gid="+gid+"#3")
    html = urlopen(url)
    soup =BeautifulSoup(html, 'lxml')
    text = soup.get_text()

    c=text.find("Completeness")
    c=c+len("completeness")+1
    comp=text[c:c+4]
    return(comp)


# In[3]:


seen=[]
tax=[]
taxa=[]
results=[]

with open ("/home/mescobar/Escritorio/Mauro/Acidiferrobacter_thyooxidans/Acidiferrobacter_thiooxydans_WP_114282240_Results.csv") as file:
    for line in file:
        line=line.split("\t")
        results.append({'Query':line[0], 'Hit ID':line[1], 'Hit Protein name':line[2], 'Taxon':line[3], 'Length':line[4], 'E-value':line[5], '-2 Protein name':line[6], '-2 Conserved Domain':line[7], '-2 Protein ID':line[8], '-1 Protein name':line[9], '-1 Conserved Domain':line[10], '-1 Protein ID':line[11], 'Hit Protein name':line[12], 'Hit Protein Conserved Domain':line[13], '+1 Protein name':line[14], '+1 Conserved Domain':line[15], '+1 Protein ID':line[16], '+2 Protein name':line[17], '+2 Conserved Domain':line[18], '+2 Protein ID':line[19] })
        if line[3] not in seen:
            tax.append(line[3])
            seen.append(line[3])
seen=[]

for t in tax:
    t=t.replace(" ", "+")
    taxa.append(t)
taxa.pop(0)
tax.pop(0)


# In[5]:


results


# In[6]:


url_List=[]
for taxon in taxa:
    url_List.append("https://gtdb.ecogenomic.org/searches?q=%25"+taxon+"%25&s=al")


# In[7]:


comp_List=[]
cont=0
for url in url_List:
    comp_List.append({'Taxon': tax[cont], 'Completeness':completeness(url)})
    cont=cont+1


# In[8]:


completeness_List=[]
pd.DataFrame(comp_List)


# In[9]:


comp_List=pd.DataFrame(comp_List)
comp_List


# In[10]:


results=pd.DataFrame(results)
#results = results.drop([0], axis=0)


# In[11]:


full_Results=pd.merge(results, comp_List)
full_Results=full_Results.replace('\\n', '', regex=True)
columnsTitles = ['Query', 'Hit ID', 'Hit Protein name', 'Taxon', 'Completeness','Length', 'E-value', '-2 Protein name', '-2 Conserved Domain', '-2 Protein ID', '-1 Protein name', '-1 Conserved Domain', '-1 Protein ID', 'Hit Protein name', 'Hit Protein Conserved Domain', '+1 Protein name', '+1 Conserved Domain', '+1 Protein ID', '+2 Protein name', '+2 Conserved Domain', '+2 Protein ID' ]
full_Results = full_Results.reindex(columns=columnsTitles)
full_Results


# In[61]:


full_Results.to_csv('Acidiferrobacter_thiooxydans_FULL.csv', index=False, sep='\t', header= True, na_rep='NA')

