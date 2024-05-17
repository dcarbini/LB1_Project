#!/bin/bash
"""DENYS CARBINI - LM BIOINFORMATICS - LB1 PROJECT"""

#AIM: Implement a HMM representing Kunitz Domain based on MSA based on structural information.

"""1. COLLECTION OF DATA"""
#PDB [https://www.rcsb.org/search/advanced] -> collection of structures#
'''
PFAM:
( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Data Collection Resolution <= 3 AND Polymer Entity Sequence Length = [ 50 - 80 ] AND Polymer Entity Mutation Count = 0

131 Structures.

CATH:
( Lineage Identifier = "4.10.410.10" AND Annotation Type = "CATH" ) AND Data Collection Resolution <= 3 AND Polymer Entity Sequence Length = [ 50 - 80 ] AND Polymer Entity Mutation Count = 0

137 Structures.

INTERPRO:
( Lineage Identifier = "IPR036880" AND Annotation Type = "InterPro" ) AND Data Collection Resolution <= 3 AND Polymer Entity Sequence Length = [ 50 - 80 ] AND Polymer Entity Mutation Count = 0

140 Structures.

#GROUP THEM: returning representative of polymers entities clusters based on 70% of Sequence identity:
PFAM: 21 Representatives
CATH: 23 Representatives
INTERPRO: 24 Representatives

#Then create a Customize table and download it as CSV:
Entity ID (by default), Resolution (Å), Sequence, Auth Asym ID, Entry ID
'''

#Clean the custom report: remove double quotes around elements and deleting the lines without the entity ID because they are the same chain of the protein above and print them as a FASTA file.
cat ./PFAM/rcsb_pdb_custom_report_20240509023402.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!="") {print ">"$5"_"$4,"\n",$3}}' > ./PFAM/pdb_seq.fasta
cat ./CATH/rcsb_pdb_custom_report_20240509023512.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!="") {print ">"$
5"_"$4,"\n",$3}}' > ./CATH/pdb_seq.fasta
cat ./INTERPRO/rcsb_pdb_custom_report_20240509023548.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!="") {print ">"$5"_"$4,"\n",$3}}' > ./INTERPRO/pdb_seq.fasta

#Check
grep '>' PFAM/pdb_seq.fasta |wc
    # 21      21     189
grep '>' CATH/pdb_seq.fasta |wc
    # 23      23     207
grep '>' INTERPRO/pdb_seq.fasta |wc
    # 24      24     216
	
#Create the list of ids
cat ./PFAM/rcsb_pdb_custom_report_20240509023402.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!=""){print $5":"$4}}' > ./PFAM/pdb_seq.ids
cat ./CATH/rcsb_pdb_custom_report_20240509023512.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!=""){print $5":"$4}}' > ./CATH/pdb_seq.ids
cat ./INTERPRO/rcsb_pdb_custom_report_20240509023548.csv | tr -d '"' |tail -n +3|awk -F "," '{if ($1!=""){print $5":"$4}}' > ./INTERPRO/pdb_seq.ids
	
#-------------------#	
	
#BENCHMARK SETS -> SWISSPROT [https://ftp.uniprot.org/pub/databases/uniprot/] -> collection of sequences#
'''
Release of the db used:
UniProt Knowledgebase Release 2024_02 consists of:
UniProtKB/Swiss-Prot Release 2024_02 of 27-Mar-2024
UniProtKB/TrEMBL Release 2024_02 of 27-Mar-2024

Query for the Positives:
(reviewed:true) AND ((xref:gene3d-4.10.410.10) OR (xref:pfam-PF00014) OR (xref:interpro-IPR036880))
401 results

Query for the Negatives:
(reviewed:true) NOT ((xref:gene3d-4.10.410.10) OR (xref:pfam-PF00014) OR (xref:interpro-IPR036880))
570,881 results
'''
#rename the positives file:
mv mv /mnt/c/Users/carbi/Downloads/uniprotkb_reviewed_true_AND_xref_pfam_P_2024_05_09.fasta.gz ./swiss_positives.fasta.gz

#Count the positives sequences
zcat swiss_positives.fasta.gz |grep ">" |wc
	# 401

gunzip swiss_positives.fasta.gz

#Create a file with the positive ids
grep ">" swiss_positives.fasta |cut -d " " -f 1 | tr -d ">" | sort > all.ids

#rename the negatives file:
mv mv /mnt/c/Users/carbi/Downloads/uniprotkb_reviewed_true_NOT_xref_pfam_P_2024_05_09.fasta.gz ./swiss_negatives.fasta.gz

#Count the negatives sequences
zcat swiss_negatives.fasta.gz |grep ">" |wc
	# 570881

zcat -f swiss_negatives.fasta.gz |grep "^>" |cut -d "|" -f 2 > negatives.ids

#Unify the structures used for the 3 HMMs
cat ./PFAM/pdb_seq.fasta ./C
ATH/pdb_seq.fasta ./INTERPRO/pdb_seq.fasta > all_pdb.fasta
#68 seq, 26 are unique

#Create the database for blast 
./../../Programs/ncbi-blast-2.15.0+/bin/makeblastdb -in all_pdb.fasta -dbtype prot

#Run blastp to find the matches
./../../Programs/ncbi-blast-2.15.0+/bin/blastp -query swiss_positives.fasta -db all_pdb.fasta -out swiss_positives.blast -outfmt 7 &

#Take only sequences with a similarity above 95% and at least 55 positions aligned
grep -v "^#" swiss_positives.blast | awk '{if ($3>95 && $4>=55) {print $0}}'|sort -nk 4 |cut -f 1 |sort -u |wc
    #35  -> sequences to be removed 

#Put the ids to remove in a file
grep -v "^#" swiss_positives.blast | awk '{if ($3>95 && $4>=55) {print $0}}'|sort -nk 4 |cut -f 1 |sort -u > tobe_removed.ids 

#Retrieve the ids of the sequences to keep
comm -23 all.ids tobe_removed.ids|cut -d "|" -f 2 > selected.ids

#check
wc selected.ids
 #366 (401-35)
wc negatives.ids
 #570881
 
#Retrieve the sequences through the python script
py getSeqs_new.py --seq swiss_positives.fasta --ids selected.ids --pos 2 > positive_selected.fasta

#Check
grep ">" positive_selected.fasta |wc
    #366 
	
#SUMMARY: 
wc selected.ids
 #366 (401-35)
wc negatives.ids
 #570881

"""2. STRUCTURE MULTIPLE ALIGNMENT"""
#https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver
#Upload the file in the format id:chain
#Download the multiple structural alignment in the kunitz_3d.aln file

#Clean the file: put the sequences in one line
awk '{if (substr($0,1,1)==">") {print "\n"$1} else {printf $1}}' ./PFAM/kunitz_3d.aln > ./PFAM/clean_kunitz_3d.aln
awk '{if (substr($0,1,1)==">") {print "\n"$1} else {printf $1}}' ./CATH/kunitz_3d.aln > ./CATH/clean_kunitz_3d.aln
awk '{if (substr($0,1,1)==">") {print "\n"$1} else {printf $1}}' ./INTERPRO/kunitz_3d.aln > ./INTERPRO/clean_kunitz_3d.aln

"""3. HMM"""
#hmmbuild output input

#TEST
hmmbuild ./PFAM/clean_kunitz_3d.hmm ./PFAM/clean_kunitz_3d.aln 
hmmbuild ./CATH/clean_kunitz_3d.hmm ./CATH/clean_kunitz_3d.aln 
hmmbuild ./INTERPRO/clean_kunitz_3d.hmm ./INTERPRO/clean_kunitz_3d.aln 

#Clean the aln 
#The hmm starts from residue 20 and continues for 59 positions
awk '{if (substr($0,1,1)==">"){print $0} else {print toupper(substr($0,20,59))}}' ./PFAM/clean_kunitz_3d.aln > ./PFAM/cut_kunitz_3d.aln
awk '{if (substr($0,1,1)==">"){print $0} else {print toupper(substr($0,19,19))toupper(substr($0,39,20))toupper(substr($0,61,9))toupper(substr($0,75,12))}}' ./CATH/clean_kunitz_3d.aln > ./CATH/cut_kunitz_3d.aln
awk '{if (substr($0,1,1)==">"){print $0} else {print toupper(substr($0,20,34))toupper(substr($0,55,18))toupper(substr($0,74,7))}}' ./INTERPRO/clean_kunitz_3d.aln > ./INTERPRO/cut_kunitz_3d.aln

#Definitive HMM
hmmbuild ./PFAM/pfam.hmm ./PFAM/cut_kunitz_3d.aln 
#final length 59 positions
hmmbuild ./CATH/cath.hmm ./CATH/cut_kunitz_3d.aln
#final length 60 positions
hmmbuild ./INTERPRO/interpro.hmm ./INTERPRO/cut_kunitz_3d.aln 
#final length 59 positions

"""4. VALIDATION AND TEST"""
#Division of positive and negatives in train and test sets
wc selected.ids
 #366 (401-35)
wc negatives.ids
 #570881
 
'''
SET			#POSITIVES (366)	#NEGATIVES (570881)
Train			256					399616
  Set1_train	64					99904
  Set2_train	64					99904
  Set3_train	64					99904
  Set4_train	64					99904
Test			110					171265
'''

#Random shuffle the selected ids
sort -R selected.ids > r1_selected.ids
sort -R negatives.ids > r1_negatives.ids

#First division of the positives and negatives between train (70%) and test (30%)
head -n 256 r1_selected.ids > pos_train.ids
tail -n +257 r1_selected.ids > pos_test.ids

head -n 399616 r1_negatives.ids > neg_train.ids
tail -n +399617 r1_negatives.ids > neg_test.ids

#Check if they don't overlap
comm -12 <(sort pos_train.ids) <(sort pos_test.ids)
	#empty
comm -23 <(sort pos_train.ids) <(sort pos_test.ids)|wc
	#256
comm -13 <(sort pos_train.ids) <(sort pos_test.ids)|wc
	#110

comm -12 <(sort neg_train.ids) <(sort neg_test.ids)
	#empty
comm -23 <(sort neg_train.ids) <(sort neg_test.ids)|wc
	#399616
comm -13 <(sort neg_train.ids) <(sort neg_test.ids)|wc
	#171265

#Divide the train in the 4 sets for the 4-fold cross-validation
 
head -n 128 pos_train.ids | head -n 64 > pos_train_1.ids
head -n 128 pos_train.ids | tail -n 64 > pos_train_2.ids
tail -n +129 pos_train.ids | head -n 64 > pos_train_3.ids
tail -n +129 pos_train.ids | tail -n 64 > pos_train_4.ids

head -n 199808 neg_train.ids | head -n 99904 > neg_train_1.ids
head -n 199808 neg_train.ids | tail -n 99904 > neg_train_2.ids
tail -n 199808 neg_train.ids | head -n 99904 > neg_train_3.ids
tail -n 199808 neg_train.ids | tail -n 99904 > neg_train_4.ids

#From the ids retrieve the sequences and create the sets of sequences
##Train sets
for i in `seq 1 4`; do py getSeqs_new.py --seq positive_selected.fasta --ids pos_train_$i.ids --pos 1 > train_pos_$i.fasta; py getSeqs_new.py --seq <(zcat -f swiss_negatives.fasta.gz) --ids neg_train_$i.ids --pos 2 > train_neg_$i.fasta; cat train_neg_$i.fasta train_pos_$i.fasta > train_$i.fasta; done

##Test sets
py getSeqs_new.py --seq positive_selected.fasta --ids pos_test.ids --pos 1 > test_pos.fasta

py getSeqs_new.py --seq <(zcat -f swiss_negatives.fasta.gz) --ids neg_test.ids --pos 2 > test_neg.fasta

cat test_pos.fasta test_neg.fasta > test.fasta

#Checks:
#All train sets:
for i in `seq 1 4`; do echo Train set $i; grep ">" train_$i.fasta | wc; done
	'''
	Train set 1
	  99968  (64+99904)
	Train set 2
	  99968   
	Train set 3
	  99968   
	Train set 4
	  99968   
	'''

#Test set
grep ">" test.fasta |wc
	#171375 (110+171265)

#hmmsearch
for j in PFAM INTERPRO CATH; do hmmsearch -Z 1 --domZ 1 --noali --max --tblout ./$j/test_pos.out ./$j/$j.hmm test_pos.fasta; hmmsearch -Z 1 --domZ 1 --noali --max --tblout ./$j/test_neg.out ./$j/$j.hmm test_neg.fasta; for i in `seq 1 4`; do hmmsearch -Z 1 --domZ 1 --noali --max --tblout ./$j/train_pos_$i.out ./$j/$j.hmm train_pos_$i.fasta; hmmsearch -Z 1 --domZ 1 --noali --max --tblout ./$j/train_neg_$i.out ./$j/$j.hmm train_neg_$i.fasta; done; done

for j in PFAM INTERPRO CATH; do echo $j: Test Neg; grep -v "^#" ./$j/test_neg.out |wc; for i in `seq 1 4`; do echo $j: Train Neg $i; grep -v "^#" ./$j/train_neg_$i.out |wc; done; done

#CREATE THE PREDS FILES
#For the positives is easy, all have an e-value from hmmsearch and we create a file with 3 columns: id, e-value and label
for j in PFAM INTERPRO CATH; do grep -v "^#" ./$j/test_pos.out | awk '{print $1"\t"$8"\t"1}' >  test_$j.preds; for i in `seq 1 4`; do grep -v "^#" ./$j/train_pos_$i.out | awk '{print $1"\t"$8"\t"1}' >  ./$j/train_$i.preds; done; done

#For the negative we put them in a temporary file
for j in PFAM INTERPRO CATH; do grep -v "^#" ./$j/test_neg.out | awk '{print $1"\t"$8"\t"0}' |sort -grk 2 > ./$j/tmp_test_neg.preds; for i in `seq 1 4`; do grep -v "^#" ./$j/train_neg_$i.out | awk '{print $1"\t"$8"\t"0}' > ./$j/tmp_train_neg_$i.preds; done; done

#Identify the missing sequences 
for j in PFAM INTERPRO CATH; do echo $j: Test Neg; comm -23 <(sort neg_test.ids) <(cut -f 1 ./$j/tmp_test_neg.preds | sort) |wc; for i in `seq 1 4`; do echo $j: Train Neg $i; comm -23 <(sort neg_train_$i.ids) <(cut -f 1 ./$j/tmp_train_neg_$i.preds | sort) |wc; done; done

#Add the missing sequences to the temporary files (e-value default threshold is 1, we put 10 as value to the one not retrieved by hmmsearch)
for j in PFAM INTERPRO CATH; do comm -23 <(sort neg_test.ids) <(cut -f 1 ./$j/tmp_test_neg.preds | sort) | awk '{print $1"\t10\t0"}' >> ./$j/tmp_test_neg.preds; for i in `seq 1 4`; do comm -23 <(sort neg_train_$i.ids) <(cut -f 1 ./$j/tmp_train_neg_$i.preds | sort) | awk '{print $1"\t10\t0"}' >> ./$j/tmp_train_neg_$i.preds; done; done

#Unify the temporary files with the preds files with the positives
for j in PFAM INTERPRO CATH; do cat ./$j/tmp_test_neg.preds >> test_$j.preds; for i in `seq 1 4`; do cat ./$j/tmp_train_neg_$i.preds >> ./$j/train_$i.preds; done; done

#Checks
wc test_*.preds
wc ./*/train_*.preds

#-----------------------------#

##CROSS-VALIDATION
#Validation
for j in PFAM INTERPRO CATH; do for i in `seq 1 15`; do py performance_new.py --th 1e-$i <(cat ./$j/train_1.preds ./$j/train_2.preds ./$j/train_3.preds) ; done > ./$j/train_123.res; for i in `seq 1 15`; do py performance_new.py --th 1e-$i <(cat ./$j/train_1.preds ./$j/train_2.preds ./$j/train_4.preds) ; done > ./$j/train_124.res; for i in `seq 1 15`; do py performance_new.py --th 1e-$i <(cat ./$j/train_1.preds ./$j/train_3.preds ./$j/train_4.preds) ; done > ./$j/train_134.res; for i in `seq 1 15`; do py performance_new.py --th 1e-$i <(cat ./$j/train_2.preds ./$j/train_3.preds ./$j/train_4.preds) ; done > ./$j/train_234.res; done

#Train Test
for j in PFAM INTERPRO CATH; do echo Train 4 $j:; cat ./$j/train_123.res | sort -nrk 6 |cut -d ' ' -f 2 |head -n 1; echo Train 3 $j:; cat ./$j/train_124.res | sort -nrk 6 |cut -d ' ' -f 2 |head -n 1; echo Train 2 $j:; cat ./$j/train_134.res | sort -nrk 6 |cut -d ' ' -f 2 |head -n 1; echo Train 1 $j:; cat ./$j/train_234.res | sort -nrk 6 |cut -d ' ' -f 2 |head -n 1; done

'''
Train 4 PFAM:
1e-06
Train 3 PFAM:
1e-06
Train 2 PFAM:
1e-06
Train 1 PFAM:
1e-06
Train 4 INTERPRO:
1e-06
Train 3 INTERPRO:
1e-06
Train 2 INTERPRO:
1e-06
Train 1 INTERPRO:
1e-06
Train 4 CATH:
1e-06
Train 3 CATH:
1e-06
Train 2 CATH:
1e-05
Train 1 CATH:
1e-05
'''
for j in PFAM INTERPRO; do for i in 1 2 3 4; do py performance_new.py --th 1e-06 ./$j/train_$i.preds >> Train_Tests_$j.res; done; done
for i in 1 2; do py performance_new.py --th 1e-06 ./CATH/train_$i.preds >> Train_Tests_CATH.res; done
for i in 3 4; do py performance_new.py --th 1e-05 ./CATH/train_$i.preds >> Train_Tests_CATH.res; done 

#GLOBAL PERFORMANCES
for j in PFAM INTERPRO CATH; do echo $j:; cat Train_Tests_$j.res | sort -nrk 6 |cut -d ' ' -f 2 |head -n 1; done 

'''
PFAM:
1e-06
INTERPRO:
1e-06
CATH:
1e-06
'''

for j in PFAM INTERPRO CATH; do echo Test $j >> Global_performances.res; py performance_new.py --cm --th 1e-06 test_$j.preds >>  Global_performances.res; echo -------------------- >>  Global_performances.res; done 

##Analysis of False positive and False negatives: see if the annotation is wrong or the HMM doesn't recognise the domain

#Positive sorted from the highest to the lowest -> check the positives with an e-value above the threshold:
for j in PFAM INTERPRO CATH; do awk '{if ($2>1e-06 && $3==1) print $0}' test_$j.preds; done 
#The False Negative (predicted negative but from the label is positive) is the same in all the models: P36235 

#Negatives sorted from the lowest to the highest -> check the negatives with an e-value below the threshold:
for j in PFAM INTERPRO CATH; do awk '{if ($3==0 && $2<1e-06) print $0}' test_$j.preds; done 
#The False Positive is the same in the 2 models that have a false postive: P0DJ63 -> problem of annotation
#In the model without the false poositive (PFAM) this sequences has an e-value of 1.2e-06.
