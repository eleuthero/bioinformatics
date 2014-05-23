import sys
import csv

with open('Finaldata_Sheet1.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
                print row



##file = open("coord.txt", 'r')
##
##row = file.readlines()
##
##for line in row:
##	print line
