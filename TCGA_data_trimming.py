#!/usr/bin/python
#-*-coding:utf8-*-
#import os
import sys
#os.chdir("C://User/Ning/Desktop/melanoma_data") 
def trimming():
#
    file = sys.argv[1]
    source_con = dict()
	
    source_text = open(file,'r')
    header = source_text.readline()
    output_file = open('trimmed_data.txt','w')
	
    for line in source_text:
        line_li = line[:-2].split('\t')
        source_con[line_li[0]] = line_li[1:]
    for key in source_con.keys():
        count = 0
        for ele in source_con[key]:
            if float(ele) > 1:
                count += 1
        if count <= 6:
            del source_con[key]
    output_file.write(header)
    for key in source_con.keys():
        output_file.write(key+'\t'+'\t'.join(source_con[key])+'\n')
    source_text.close()
    output_file.close()
trimming()