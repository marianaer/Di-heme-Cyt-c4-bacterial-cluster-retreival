#!/bin/bash

ids='./IDblast.txt'
while read line; #iterando sobre cada línea
do
	
	if [[ $line =~ ^WP ]]
	then 
		efetch -db protein -id $line -format ipg |awk -F "\t" '$2=="RefSeq" {print $3}'| xargs -n 1 sh -c 'efetch -db nuccore -id "$0" -format "gbc"'|xtract -insd CDS protein_id product |grep -A 2 -B 2 $line |sed "s/$/\t$line/"| head -5 >> FlankingProteins.txt


	#si no es no redundante. El sed de al final agrega una última columna con el nombre de la proteína de en medio (para identificar)
	else
		efetch -db protein -id $line -format ipg | awk -F "\t" '$2=="INSDC" {print $3}'| xargs -n 1 sh -c 'efetch -db nuccore -id "$0" -format "gbc"'| xtract -insd CDS  protein_id product| grep -A 2 -B 2 $line |sed "s/$/\t$line/" | head -5 >> FlankingProteins.txt
	fi

done < $ids
