#!/bin/bash
import os,re,sys
import xml.etree.ElementTree as ET
reload(sys)
sys.setdefaultencoding('utf-8')
   
rootDir = os.getcwd() + '/raw_clinical/'
outputDir = os.getcwd() + '/clinical_trimmed/'
def xml2csv(rootDir,outputDir):
    file_list = file_iter(rootDir)
    for fi in file_list:
        tree = ET.parse(rootDir+fi)
        root = tree.getroot()
        output = open(outputDir+fi[-16:-4]+'.csv','w')
        output.write('SampleID,'+fi[-16:-4]+'\n')
        for child1 in root:
            for child2 in child1:
                output.write(re.findall('(?<=})\w*',child2.tag)[0]+ ',' + str(child2.text) + '\n')
                for child3 in child2:
                    output.write(re.findall('(?<=})\w*',child3.tag)[0]+ ',' + str(child3.text) + '\n')
                    for child4 in child3:
                        output.write(re.findall('(?<=})\w*',child4.tag)[0]+ ',' + str(child4.text) + '\n')
                        for child5 in child4:
                            output.write(re.findall('(?<=})\w*',child5.tag)[0]+ ',' + str(child5.text) + '\n')
    output.close()
    return

def file_iter(rootDir):
    file_list = []
    for i in os.listdir(rootDir):
        if os.path.isfile(os.path.join(rootDir,i)):
            file_list += [i]
    return(file_list)
	
xml2csv(rootDir,outputDir)