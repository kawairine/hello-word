from sklearn import svm
import numpy as np

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


### residues with window
seqindividuallist=[]
for seq in seqlist:
    sequpp=seq.upper()
    sequpp='0'+seq+'0'
    seqindividuallist.append(sequpp.split(','))

##2 train vs 1 test sets :
train_threefirst=[]
train_threesecond=[]
test_threefirst=[]

for e_seq1 in seqindividuallist[0:15]:
    train_threefirst.append(e_seq1)
for e_seq2 in seqindividuallist[15:30]:
    train_threesecond.append(e_seq2)
for e_seq3 in seqindividuallist[30:43]:
    test_threefirst.append(e_seq3)


trifortrain = []
trifortest=[]
for eachseq_1 in train_threefirst:
    for eachaa_1 in eachseq_1:
        aatri_1=[eachaa_1[i:i+3] for i in range(len(eachaa_1)-2)]
        trifortrain.append(aatri_1)


for eachseq_2 in train_threesecond:
    for eachaa_2 in eachseq_2:
        aatri_2=[eachaa_2[i:i+3] for i in range(len(eachaa_2)-2)]
        trifortrain.append(aatri_2)


trifortrain = [j for i in trifortrain for j in i]


for eachseq_3 in test_threefirst:
    for eachaa_3 in eachseq_3:
        aatri_3=[eachaa_3[i:i+3] for i in range(len(eachaa_3)-2)]
        trifortest.append(aatri_3)
trifortest = [j for i in trifortest for j in i]

###############convert to numerical values

mapa={'0':0,'R':1,'H':2,'K':3,'D':4,'E':5,'S':6,'T':7,
                    'N':8,'Q':9,'C':10,'U':11,'G':12,'P':13,
                    'A':14,'V':15,'I':16,'L':17,'M':18,'F':19,
                    'Y':20,'W':21}

mappedtrifortrain=[]
for everywindow in trifortrain:
    tmp=[]
    for everyres in everywindow:
        tmp1=mapa[everyres]
        tmp.append(tmp1)
    mappedtrifortrain.append(tmp)

mappedtrifortest=[]
for everywindowtest in trifortest:
    tmptest=[]
    for everyrestest in everywindowtest:
        tmp2=mapa[everyrestest]
        tmptest.append(tmp2)
    mappedtrifortest.append(tmptest)

#### convert to 01 format
from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder()
mappedtrifortrainarray=enc.fit_transform(mappedtrifortrain).toarray()

from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder()
mappedtrifortestarray=enc.fit_transform(mappedtrifortrain).toarray()


### 2nd structure
dicfeastates={"i":0,"o":1,"l":2,"p":2,"m":2}

#convert each as a list and convert to numerical term
fea1=[]
fea2=[]
fea3=[]
for e_fea1 in fea[0:15]:
    fea1.append(e_fea1)
for e_fea2 in fea[15:30]:
    fea2.append(e_fea2)
for e_fea3 in fea[30:43]:
    fea3.append(e_fea3)

featrain=[]
for feaindividual1 in fea1:
    for eachfeature1 in feaindividual1:
        fealower1=eachfeature1.lower()
        tempo1=[]
        tempoo1=dicfeastates[fealower1]
        tempo1.append(tempoo1)
        featrain.append(tempo1)
for feaindividual2 in fea2:
    for eachfeature2 in feaindividual2:
        fealower2=eachfeature2.lower()
        tempo2=[]
        tempoo2=dicfeastates[fealower2]
        tempo2.append(tempoo2)
        featrain.append(tempo2)

featest=[]
for feaindividual3 in fea3:
    for eachfeature3 in feaindividual3:
        fealower3=eachfeature3.lower()
        tempo3=[]
        tempoo3=dicfeastates[fealower3]
        tempo3.append(tempoo3)
        featest.append(tempo3)


#########################################################
### if doesn't take the last and last residue...
# featurewindowtotal=[]
# for seqfea in feaindividuallist:
#     for aafea in seqfea:
#         featurewindow=[aafea[1:-1]]
#         featurewindowtotal.append(featurewindow)
#########################################################


#### CV sets
#
# from sklearn.model_selection import StratifiedKFold
# X = np.array(mappedbigtri)
# y = np.array(feasingle)
# skf = StratifiedKFold(n_splits=3)
# skf.get_n_splits(X, y)
# for train_index, test_index in skf.split(X, y):
#     print("TRAIN:", train_index, "TEST:", test_index)
#     X_train, X_test = X[train_index], X[test_index]
#     y_train, y_test = y[train_index], y[test_index]
# # TRAIN: [1 3] TEST: [0 2]
# TRAIN: [0 2] TEST: [1 3]




#

# aaarray=np.array(aalist)
#

#
#
#
# feaindividualarray=np.array(feaindividual)
#
# array2= np.array([aaarray,feaindividualarray])
# array2= array2.T

# print (aalist)
# with open('featurefile.txt', 'w+') as fwrite:
#      for state in feastatesmap:
#          for A in aalist:
#              fwrite.write(state+A+'\n')

# np.column_stack(aaarray,feaindividualarray)


# xsplit=[]
# for seq in aalist:
#     xsplit.append(seq.split(','))
#
# xsplitsplit=[]
# for AA in xsplit:
#     xsplitsplit.append(AA.split())
#     # num=ord(AA)
#     # lisAA.append(num)
#
# print xsplitsplit
# # x='MSL'
# y = 'iow'
# xsplit=list(x)
# ysplit=list(y)
# lisAA=[]
# lisfeatures=[]
# for AA in xsplit:
    #  num=ord(AA)
    #  lisAA.append(num)
#
# for features in ysplit:
#     lisfeatures.append(features)
# #print (lis)
# arrAA = np.array(lisAA).reshape(-1,1)
# arrfea= np.array(lisfeatures)
# #print (arr)
# X = [[1,2], [2]]
# y = ['i', 'o']
# clf = svm.SVC()
# clf.fit(X, y)
# #
# #
# #
# #
# #
# print (clf.predict([1]))

#print ord(xsplit)
