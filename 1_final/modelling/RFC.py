from sklearn import svm
from sklearn.ensemble import RandomForestClassifier as RFC
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
            # print(freq)
            # print(freq[-1])
            # freq[1][:0]=zeroes
            # freq[-1].extend(zeroes)
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
# eachprotein=[]

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
            # eachprotein=[j for i in eachprotein for j in i]
        # eachset.append(eachprotein)
        # eachset=[j for i in eachset for j in i]
        # eachprotein=[]
    temptrainingset.append(eachset)


for j in one_seq_test:
    pears=[j for i in one_seq_test for j in i]

proteintest=[]
for seqtest in pears:
    # sys.exit()
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
#

# print (major_fea_train) #1st level set, 2nd level is sequence in the set
flattenseqinset = []
for k in major_fea_train: #for every set --> k, flatten
    for everyseq in k:
        flatseqinset = [everyseq for i in k for everyseq in i]
    flattenseqinset.append(flatseqinset)


# print(major_fea_train[0])
# print(len(major_fea_train))
# print(len(major_fea_train[0]))


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
clf=RFC(n_estimators=10,class_weight="balanced")
scorelist=[]
p=[]
r=[]
f=[]
sup=[]
# clf=svm.SVC(kernel='rbf')
codelabel=[0,1,2]
from sklearn.metrics import precision_recall_fscore_support as scoreprfs
# from sklearn.metrics import classification_report
target_names=[]

for i in range(0,len(temptrainingset)):
    # print(len(temptrainingset[i]),len(fea_train_settemp[i]))
    clf.fit(temptrainingset[i],fea_train_settemp[i])
    pred=clf.predict(proteintest[i])
    predictresult=list(pred)
    # predictionlist.append(predictresult)
    # truelist.extend(fea_test_set[i])
    s=float(clf.score(proteintest[i],fea_test_set[i]))
    scorelist.append(s)
    precision,recall,fscore,support=scoreprfs(fea_test_set[i],pred,labels=codelabel)
    # classification_report(fea_test_set[i],pred)
    # p.append(list(precision))
    r.append(list(recall))
    f.append(list(fscore))
    # # sup.append(list(support))
# TP=0
# FP=0
# TN=0
# FN=0
# for tf in range(len(predictionlist)):
#     if truelist[i]==predictionlist[i]


print ("RFC balanced, window=%s" %(ws))
import statistics
trainingsamplestd=round(statistics.stdev(scorelist),3)
print("standard")
print(trainingsamplestd)
averagescore=round(sum(scorelist)/len(scorelist),3)
print("averagescore")
print(averagescore)
#
# # p_class0=[]
# # p_class1=[]
# # p_class2=[]
r_class0=[]
r_class1=[]
r_class2=[]
f_class0=[]
f_class1=[]
f_class2=[]
# # sup_class0=[]
# # sup_class1=[]
# # sup_class2=[]
# #
for i in range(0,len(f)):
# #     p_class0.append(p[i][0])
# #     p_class1.append(p[i][1])
# #     p_class2.append(p[i][2])
    r_class0.append(r[i][0])
    r_class1.append(r[i][1])
    r_class2.append(r[i][2])
    f_class0.append(f[i][0])
    f_class1.append(f[i][1])
    f_class2.append(f[i][2])
# #     sup_class0.append(sup[i][0])
# #     sup_class1.append(sup[i][1])
# #     sup_class2.append(sup[i][2])
# #
# # p_class0avg=statistics.mean(p_class0)
# # p_class1avg=statistics.mean(p_class1)
# # p_class2avg=statistics.mean(p_class2)
r_class0avg=statistics.mean(r_class0)
r_class1avg=statistics.mean(r_class1)
r_class2avg=statistics.mean(r_class2)
f_class0avg=statistics.mean(f_class0)
f_class1avg=statistics.mean(f_class1)
f_class2avg=statistics.mean(f_class2)
# # sup_class0avg=statistics.mean(sup_class0)
# # sup_class1avg=statistics.mean(sup_class1)
# # sup_class2avg=statistics.mean(sup_class2)
# #
# # print (p_class0avg)
# # print (p_class1avg)
# # print (p_class2avg)
print ('r_class0avg')
print (r_class0avg)
print ('r_class1avg')
print (r_class1avg)
print ('r_class2avg')
print (r_class2avg)
print ('f_class0avg')
print (f_class0avg)
print ('f_class1avg')
print (f_class1avg)
print ('f_class2avg')
print (f_class2avg)
# print (sup_class0avg)
# print (sup_class1avg)
# print (sup_class2avg)
#
#
#
# # #     # accuracy average
# #     # -> whole set -> confusion matrix
#
#
# # from sklearn.metrics import precision_recall_fscore_support as score
# #
#
# # print('precision:{}'.format(precision))
# # print('recall:{}'.format(recall))
# # print('fscore:{}'.format(fscore))
# # print('support:{}'.format(support))
#
# #
# # f1_score(truelist,predictionlist,average='weighted')
# # #
# # import itertools
# # import matplotlib.pyplot as plt
# # from sklearn.metrics import confusion_matrix
# # def plot_confusion_matrix(cm, classes,
# #                           normalize=True,
# #                           title='Confusion matrix',
# #                           cmap=plt.cm.Blues):
# #     plt.imshow(cm, interpolation='nearest', cmap=cmap)
# #     plt.title(title)
# #     plt.colorbar()
# #     tick_marks = np.arange(len(classes))
# #     plt.xticks(tick_marks, classes, rotation=45)
# #     plt.yticks(tick_marks, classes)
# #
# #     if normalize:
# #         cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
# #         print("Normalized confusion matrix")
# #     else:
# #         print('Confusion matrix, without normalization')
# #
# #     print(cm)
# #
# #     thresh = cm.max() / 2.
# #     for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
# #         plt.text(j, i, cm[i, j],
# #                  horizontalalignment="center",
# #                  color="white" if cm[i, j] > thresh else "black")
# #
# #     plt.tight_layout()
# #     plt.ylabel('True label')
# #     plt.xlabel('Predicted label')
# #
# # # Compute confusion matrix
# # #codelabel=[0,1,2]
# # cnf_matrix = confusion_matrix(truelist, predictionlist,labels=codelabel)
# # np.set_printoptions(precision=2)
# #
# # # Plot normalized confusion matrix
# # plt.figure()
# # plot_confusion_matrix(cnf_matrix, classes=codelabel, normalize=True,
# #                       title='Normalized confusion matrix,window size=3')
# #
# #
# # plt.show()
# # #
# # # from sklearn.utils import compute_class_weights
