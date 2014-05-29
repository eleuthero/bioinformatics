import sys
import csv
import string
import MySQLdb



db = MySQLdb.connect(host= "db.cs.wwu.edu", user = "rajends_writer" , passwd = "aqFsF7bQ", db = "rajends")
if db:
    print "connecte"
    cur = db.cursor()
else:
        print "The db could not be found"

##with open('Finaldata_Sheet1.csv','rb') as f:
##        
##        reader = csv.reader(f)
##        for row in reader:
##                print row

##file = open("coord.txt", 'r')
##
##row = file.readlines()
##
##for line in row:
##	print line
