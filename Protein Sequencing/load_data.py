import sys
import csv
import string
import MySQLdb as mdb


try:
    #Open a database connection
    con = mdb.connect('db.cs.wwu.edu', 'rajends_writer', 'aqFsF7bQ', 'rajends');
    #Prepare a cursor object using cursor() method
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS Proteins")
    #Create a table
    sql = """CREATE TABLE Proteins(HXB2_K03455 VARCHAR(3),
                 HXB2_Position VARCHAR(10),
                 frame1 VARCHAR(5),
                 frame2 VARCHAR(5),
                 frame3 VARCHAR(5),
                 aminoAcid VARCHAR(5),
                 aminoAcidDes VARCHAR(50),
                 proteinNumber VARCHAR(20),
                 gene VARCHAR(30),
                 protein VARCHAR(70),
                 RNAFeature VARCHAR(100),
                 proteinFeature VARCHAR(70),
                 secondFrame VARCHAR(5) DEFAULT NULL,
                 numbering VARCHAR(10),
                 proteinFeature2 VARCHAR(50),
                 thirdFrame VARCHAR(10),
                 numbering3 VARCHAR(6),
                 proteinFeature3 VARCHAR(70),
                 reference VARCHAR(50))"""
    cur.execute(sql)
#####################################################################
#####################################################################
##This is the section of loading in data from cvs to the db
    with open('result.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            cur.execute('INSERT INTO Proteins(HXB2_K03455,\
                                                 HXB2_Position,\
                                                 frame1,\
                                                 frame2,\
                                                 frame3,\
                                                 aminoAcid,\
                                                 aminoAcidDes,\
                                                 proteinNumber,\
                                                 gene,\
                                                 protein ,\
                                                 RNAFeature, \
                                                 proteinFeature ,\
                                                 secondFrame ,\
                                                 numbering ,\
                                                 proteinFeature2 ,\
                                                 thirdFrame ,\
                                                 numbering3 ,\
                                                 proteinFeature3 ,\
                                                 reference)'
                            'VALUES("%s","%s","%s","%s","%s","%s","%s",\
                                    "%s","%s","%s","%s","%s","%s","%s",\
                                    "%s","%s","%s","%s","%s")',row)                      
except mdb.Error, e:
    print "Error %d: %s" % (e.args[0],e.args[1])
    sys.exit(1)
if con:
    con.commit()
    con.close()
