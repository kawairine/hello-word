# print('''please input a fasta file of the Beta-barrel protein
#         that you want to predict its topology:
#         i: inner-loop, o: outer-loop, m: transmembrane Beta strand
#         path:      ''')

import sys
path=sys.argv[1]

openfile=open(path,'r')
openfile=openfile.readlines()
ID=openfile[0].rstrip('\n').lstrip('>')
ID=str(ID)
# seqstring=openfile[1]
# seq=[]
#
# seq.append(openfile[1])

path=str(path)+'.pssm'

PSIfreq=[]
with open(path, 'r') as path:
    path = path.readlines()
    for line in path[3:-6]:
        line = line.split()
        freq=line[22:42]
        # print(freq)
        # print(freq[-1])
        # freq[1][:0]=zeroes
        # freq[-1].extend(zeroes)
        PSIfreq.append(freq)
    # PSSMlist.append(PSIfreq)

PSIfreqfloat=[]
for each_position in PSIfreq:
    eachres=[]
    for eachfreq in each_position:
        eachfreq=(float(eachfreq)/100)
        # eachfreq=round(eachfreq,2)
        eachres.append(eachfreq)
    PSIfreqfloat.append(eachres)


zeropad=([float(0)]*20)

ws=31
sp=int(ws/2)
temptrainingset=[]
# eachprotein=[]

proteintest=[]

eachproteintest=[]
for i in range(0,len(PSIfreqfloat)):
    if i<sp:
        wis=([zeropad]*(sp-i)) + PSIfreqfloat[0:ws-(sp-i)]
        wis=[a for b in wis for a in b]
        eachproteintest.append(wis)
    elif i>= (len(PSIfreqfloat)-sp):
        middle=PSIfreqfloat[(i-sp):ws-(sp-i)+1]
        if len(middle) != ws:
            middle.extend([zeropad]*(ws-len(middle)))
        middle=[a for b in middle for a in b]
        eachproteintest.append(middle)
    else:
        w=PSIfreqfloat[(i-sp):(i+sp+1)]
        w=[a for b in w for a in b]
        eachproteintest.append(w)#

from sklearn.externals import joblib
from sklearn import svm
clf=joblib.load('svmpredictor.pkl')
pred=clf.predict(eachproteintest)
predictionresult=list(pred)
dicfeastates={0:"i",1:"o",2:"m"}
predictionstring=[]
for eachpred in predictionresult:
    eachpred=dicfeastates[eachpred]
    predictionstring.append(eachpred)
predstring=''.join(predictionstring)
#

print("The prediction of %s is:  " %(ID))
print(" ")
print(predstring)
print("'i'is inner-loop, 'o' is outer-loop, 'm' is transmembrane beta strand")
print(" ")

with open('prediction_of_topology.txt', 'w') as fwrite:
    for i in range(0,len(openfile)):
        fwrite.write(openfile[i])
    fwrite.write("\n >%s_prediction_topology \n" %(ID))
    fwrite.write(predstring)
