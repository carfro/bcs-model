# Script that strips whitespaces from file, arg1 is file to read from, arg2 is file to write to
import sys
#WMAT_ORG=open('W_mat.dat','r')
#WMAT_FIN=open('W_F.dat','w')

WMAT_ORG=open(sys.argv[1],'r')
WMAT_FIN=open(sys.argv[2],'w')
for line in WMAT_ORG:
    if len(line)>2 :
        line=line.replace(" ","")
        line=line.replace(','," ")
    #print line
    WMAT_FIN.writelines(line)
