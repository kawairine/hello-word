from sklearn import svm
import numpy as np
from sklearn import datasets

####extracting the data
openfile=open('membrane-beta_3state.3line.txt','r')
seqlist=[]
fea=[]
for line in openfile:
    low=line.lower()
    if "sequence" in low:
        seqlist.append((next(openfile).rstrip('\n')))
    if "observed" in low:
        fea.append((next(openfile).lower().rstrip('\n')))

### residues with zero
seqindividuallist=[]
for seq in seqlist:
    sequpp=seq.upper()
    sequpp='0'+seq+'0'
    seqindividuallist.append(sequpp.split(','))

#convert each as a list
fealist=[]
for feaindividual in fea:
    fealower=feaindividual.lower()
    fealist.append(fealower.split(','))


###leave one out
numberofseq=len(seqindividuallist)
major_seq_train=[]
one_seq_test=[]
major_fea_train=[]
one_fea_test=[]
last=0.0

while last <numberofseq:
    one_seq_test.append(seqindividuallist[int(last):int(last+1)])
    major_seq_train.append(seqindividuallist[0:int(last)]+ seqindividuallist[int(last+1):])
    one_fea_test.append(fealist[int(last):int(last+1)])
    major_fea_train.append(fealist[0:int(last)]+ fealist[int(last+1):])
    last += 1

# print (len(major_seq_train[0]))
# print (len(major_fea_train[0]))
# print (len(one_seq_test[0]))
# print (len(one_fea_test[0]))

##slidingggggggggggggggggg
set_train = []
for j in range(0,len(major_seq_train)): #eachset
    apples = [j for i in major_seq_train[j] for j in i]

    temp_seq = []
    for eachseq in apples: #each sequence in the set

        a =[eachseq[i:i+3] for i in range(len(eachseq)-2)]

        temp_seq.append(a)

    set_train.append(temp_seq)



# print (len(set_train), len(set_train[0][0]))


# print(one_seq_test)

for j in one_seq_test:
    pears=[j for i in one_seq_test for j in i] #eachseq=a test set(presented as a list)
    seq_test = []
    for eachseqlist in pears: #eachseq as set
        for seq in eachseqlist:
            p =[seq[i:i+3] for i in range(len(seq)-2)]
            seq_test.append(p)
    # print(seq_test)
        # tritest.append(seq_test)

# print (seq_test)
# print (len(seq_test))

# ##### numerical value

mapa={'0':0,'R':1,'H':2,'K':3,'D':4,'E':5,'S':6,'T':7,
                    'N':8,'Q':9,'C':10,'W':11,'G':12,'P':13,
                    'A':14,'V':15,'I':16,'L':17,'M':18,'F':19,
                    'Y':20}

#mapping trilisttrain to numbers and sparsing encoding

from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder(n_values=21)
#
# print (set_train)
#
settrainmap=[]
for everyset in set_train:
    # print (everyset)
    mappedtri=[]
    for map_everyseq in everyset: #seq as list
        for everywindow in map_everyseq:
            tmp=[]
            for everyres in everywindow:
                tmp1=mapa[everyres]
                tmp.append(tmp1)
            mappedtri.append(tmp)
    settrainmap.append(mappedtri)
#
# print(len(settrainmap), len(settrainmap[2]))

arraytrain_list=[]
for everyset in settrainmap: #each set
    mappedarray=enc.fit_transform(everyset).toarray()
    arraytrain_list.append(mappedarray)


# #mapping tritest to numbers
settestmap=[]
for everyseq in seq_test:  #every seq as a list
    mappedtritest=[]
    for everywindow in everyseq:
        tmptest=[]
        for everyres in everywindow:
            tmp1test=mapa[everyres]
            tmptest.append(tmp1test)
        mappedtritest.append(tmptest)
    settestmap.append(mappedtritest)


# print (len(settestmap[2]))

# print (settestmap)
arraytest_list=[]
for everyset in settestmap:
    mappedarraytest=enc.fit_transform(everyset).toarray()
    arraytest_list.append(mappedarraytest)

# print (len(arraytrain_list[0]))

# ###fea########################################


# print (major_fea_train) #1st level set, 2nd level is sequence in the set
flattenseqinset = []
for k in major_fea_train: #for every set --> k, flatten
    for everyseq in k:
        flatseqinset = [everyseq for i in k for everyseq in i]
    flattenseqinset.append(flatseqinset)



# print(flattenseqinset[0])
dicfeastates={"i":0,"o":1,"l":2,"p":2,"m":2}
fea_train_set=[]
for everyset in flattenseqinset:
    setlist=[]
    for eachseq in everyset:
        for eachfea in eachseq:
            tempoo=dicfeastates[eachfea]
            setlist.append(tempoo)
    fea_train_set.append(setlist)
# print(fea_train_set[0])
# print (len(fea_train_set[0]))

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


# ##################SVM

predictionlist=[]
clf=svm.SVC()
for i in range(0,len(arraytrain_list)):
    clf.fit(arraytrain_list[i],fea_train_set[i])
    predictresult=clf.predict(arraytest_list[i])
    predictionlist.append(predictresult)

print(predictionlist)


# clf=svm.SVC()
# clf.fit(mappedtrifortrainarray,featrain)
# predictresult=clf.predict(mappedtrifortestarray)
# probability=clf.predict_proba(predictresult)
# print(probability)
