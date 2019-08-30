#! /usr/bin/python
import os
import sys

print 'Usage: replaceStringInPDF.py FileToBeReplaced.pdf oldstring newstring'
print 'Argument List:', str(sys.argv)
assert(len(sys.argv)==4)
#print sys.argv[1]

# a la https://stackoverflow.com/questions/9871585/how-to-find-and-replace-text-in-a-existing-pdf-file-with-pdftk-or-other-command
os.system("pdftk {} output uncompressed.pdf uncompress".format(sys.argv[1]))
os.system('sed -e "s/{}/{}/g" <uncompressed.pdf >modified.pdf'.format(sys.argv[2],sys.argv[3]))
os.system("pdftk modified.pdf output recompressed.pdf compress")
os.system("mv recompressed.pdf {}".format(sys.argv[1]))
os.system("rm uncompressed.pdf modified.pdf")
