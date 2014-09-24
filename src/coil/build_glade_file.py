#!/usr/bin/env python2.7
import sys
import os

inputfile='clwingtk.gladexml'
file_name=os.path.basename(inputfile)
output = open("coil/gladexmlfile.cpp", 'w+')
symbol_name=os.path.splitext(file_name)[0]
output.write("#include <string>\n")
output.write("extern const std::string "+symbol_name+";\n")
output.write("const std::string "+symbol_name+" = ")
for line in open(inputfile, 'r'):
    #Escape any backslashes
    line=line.replace('\\', '\\\\')
    #Escape any quotes
    line=line.replace('"', '\\"')
    #Write out the line wrapped in quotes with newlines in the data
    output.write('"'+line.replace('\n', '\\n"\n'))
output.write(";")
