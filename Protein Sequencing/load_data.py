################################################
#Soundarya Rajendran
#HIV Project
#This is the file that creates a table, and uploads data from
#a results.csv file with data downloaded online.
#BIOFINFORMATICS : SPRING 2014 ; HIV MUTATION SEQUENCING PROJECT
################################################5
import sys
import csv
import string
import MySQLdb as mdb
#####################################################################
#####################################################################
##This is the section of loading in data from cvs to the db
with open('result.csv','rb') as f:
    reader = csv.reader(f)
    with open('protein_data.sql','a') as writer:
        writer.seek(0)
        writer.truncate()
        for row in reader:
            test_string = ('INSERT INTO Proteins(HXB2_K03455,HXB2_Position,frame1,frame2,frame3,aminoAcid,aminoAcidDes,proteinNumber,gene,protein ,RNAFeature,proteinFeature ,secondFrame ,numbering ,proteinFeature2 ,thirdFrame ,numbering3 ,proteinFeature3 ,reference)')
            test_string+= ('VALUES(' + row[0] + ',' + row[1] + ',' + row[2] +  ','+ row[3]  + ','+ row[4] +  ','+ row[5] +  ','+ row[6] +  ','+ row[7] +  ','+ row[8] +  ','+ row[9] +  ','+ row[10] +  ','+ row[11] + ','+ row[12] +  ','+ row[13] +  ','+ row[14] +  ','+ row[15] +  ','+ row[16] +  ','+ row[17] +  ','+ row[18]+ ');\n')
            writer.write(test_string)
f.close()
writer.close()
