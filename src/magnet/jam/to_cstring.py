#!/usr/bin/env python2
import sys
import os

if len(sys.argv) != 3:
    print "Usage: ./to_cstring.py input output"
    quit()

output = open(sys.argv[2], 'w+')
file_name=os.path.basename(sys.argv[2])
symbol_name=os.path.splitext(file_name)[0]
output.write("#include <string>\n")
output.write("extern const std::string "+symbol_name+";\n")
output.write("const std::string "+symbol_name+" = ")
for line in open(sys.argv[1], 'r'):
    #Escape any backslashes
    line=line.replace('\\', '\\\\')
    #Escape any quotes
    line=line.replace('"', '\\"')
    #Write out the line wrapped in quotes with newlines in the data
    output.write('"'+line.replace('\n', '\\n"\n'))
output.write(";")
