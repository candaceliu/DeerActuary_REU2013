#!/usr/bin/python


import os
import glob
import csv
import sys

rowWriter = csv.writer(sys.stdout) #,delimter=",")
rowWriter.writerow(["time","P","alpha","x","m","sumx", \
                    "sumx2","summ","summ2","N"])

for theFile in glob.glob("*.dat"):
    fp=open(theFile,"r")
    rowReader = csv.reader(fp,delimiter=',')
    header = rowReader.next()
    for row in rowReader:
        rowWriter.writerow(row)

    fp.close()


