from sklearn import svm
import numpy as np
from sklearn import datasets
import sys

####extracting the data
openfile=open('membrane-beta_3state.3line.txt','r')
openfile=openfile.readlines()
idlist=[]
seqlist=[]
fea=[]

for i in range(0, len(openfile), 5):
    idlist.append(openfile[i].rstrip('\n'))
    seqlist.append(openfile[i+1].upper().rstrip('\n'))
    fea.append(openfile[i+3].lower().rstrip('\n'))


#### PSIBLAST
PSSMlist=[]

for ID in idlist:
    PSIfreq=[]
    ID = str(ID)
    path='./fasta/PSSM/'+ID+'.fasta.pssm'
    with open(path, 'r') as path:
        path = path.readlines()
        for line in path[3:-6]:
            line = line.split()
            freq=line[22:42]
            PSIfreq.append(freq)
        PSSMlist.append(PSIfreq)

PSIfreqfloat=[]



for eachseq_PSSM in PSSMlist:
    seqPSSM=[]
    for eachseq_position in eachseq_PSSM:
        eachseqPSSM=[]
        for eachfreq in eachseq_position:
            eachfreq=(float(eachfreq)/100)
            # eachfreq=round(eachfreq,2)
            eachseqPSSM.append(eachfreq)
        seqPSSM.append(eachseqPSSM)
    PSIfreqfloat.append(seqPSSM)

fealist=[]
for feaindividual in fea:
    fealist.append(feaindividual.split(','))

###leave one out
numberofseq=len(PSIfreqfloat)
major_seq_train=[]
one_seq_test=[]
major_fea_train=[]
one_fea_test=[]
last=0.0

while last <numberofseq:
    one_seq_test.append(PSIfreqfloat[int(last):int(last+1)])
    major_seq_train.append(PSIfreqfloat[0:int(last)]+ PSIfreqfloat[int(last+1):])
    one_fea_test.append(fealist[int(last):int(last+1)])
    major_fea_train.append(fealist[0:int(last)]+ fealist[int(last+1):])
    last += 1


zeropad=([float(0)]*20)

ws=31
sp=int(ws/2)
temptrainingset=[]

for sett in major_seq_train:
    eachset=[]
    for seq in sett:
        for i in range(0,len(seq)):
            if i<sp:
                wis=([zeropad]*(sp-i)) + seq[0:ws-(sp-i)]
                wis=[a for b in wis for a in b]
                eachset.append(wis)
            elif i>= (len(seq)-sp):
                middle=seq[(i-sp):ws-(sp-i)+1]
                if len(middle) != ws:
                    middle.extend([zeropad]*(ws-len(middle)))
                middle=[a for b in middle for a in b]
                eachset.append(middle)
            else:
                w=seq[(i-sp):(i+sp+1)]
                w=[a for b in w for a in b]
                eachset.append(w)
    temptrainingset.append(eachset)


for j in one_seq_test:
    pears=[j for i in one_seq_test for j in i]

proteintest=[]
for seqtest in pears:
    eachproteintest=[]
    for i in range(0,len(seqtest)):
        if i<sp:
            wis=([zeropad]*(sp-i)) + seqtest[0:ws-(sp-i)]
            wis=[a for b in wis for a in b]
            eachproteintest.append(wis)
        elif i>= (len(seqtest)-sp):
            middle=seqtest[(i-sp):ws-(sp-i)+1]
            if len(middle) != ws:
                middle.extend([zeropad]*(ws-len(middle)))
            middle=[a for b in middle for a in b]
            eachproteintest.append(middle)
        else:
            w=seqtest[(i-sp):(i+sp+1)]
            w=[a for b in w for a in b]
            eachproteintest.append(w)
    proteintest.append(eachproteintest)

flattenseqinset = []
for k in major_fea_train: #for every set --> k, flatten
    for everyseq in k:
        flatseqinset = [everyseq for i in k for everyseq in i]
    flattenseqinset.append(flatseqinset)


dicfeastates={"i":0,"o":1,"l":2,"p":2,"m":2}
fea_train_settemp=[]
for everyset in flattenseqinset:
    setlist=[]
    for eachseq in everyset:
        for eachfea in eachseq:
            tempoo=dicfeastates[eachfea]
            setlist.append(tempoo)
    fea_train_settemp.append(setlist)


fea_train_set=[]
for everyset in fea_train_settemp:
    temp=np.array(everyset)
    fea_train_set.append(temp)


fea_test_set=[]
for j in one_fea_test:
    banana=[j for i in one_fea_test for j in i] #eachseq=a test set(presented as a list)

for eachfealist in banana: #eachseq as set
    fea_test = []
    for fea in eachfealist:
        for eachfea in fea:
            tempo=dicfeastates[eachfea]
            fea_test.append(tempo)
    fea_test_set.append(fea_test)

fea_test_setarray=[]
for x in fea_test_set:
    arrayfea=np.array(x)
    fea_test_setarray.append(arrayfea)


predictionlist=[]
truelist = []
clf=svm.LinearSVC(class_weight="balanced",C=5.218540993)

codelabel=[0,1,2]

for i in range(0,len(temptrainingset)):
    clf.fit(temptrainingset[i],fea_train_settemp[i])

from sklearn.externals import joblib
joblib.dump(clf,'linear_svmpredictor.pkl')
