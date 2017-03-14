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


### residues with zero
seqindividuallist=[]
for seq in seqlist:
    # sequpp='0'*23+seq+'0'*23
    seqindividuallist.append(seq.split(','))

# freq[1][:0]=zeroes
# freq[-1].extend(zeroes)

#convert each as a list
fealist=[]
for feaindividual in fea:
    fealist.append(feaindividual.split(','))

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


# print(len(major_seq_train[0][0][0])) #-> residues
# print(len(major_seq_train[0][0])) #-> seq
# print(len(major_seq_train[0])) #-> set
# print(len(major_seq_train))
# sys.exit()

mapa={'0':0,'R':1,'H':2,'K':3,'D':4,'E':5,'S':6,'T':7,
                    'N':8,'Q':9,'C':10,'W':11,'G':12,'P':13,
                    'A':14,'V':15,'I':16,'L':17,'M':18,'F':19,
                    'Y':20}

from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder(n_values=21)

zeropad=[0]

ws=77
sp=int(ws/2)
temptrainingset=[]
eachprotein=[]

for sett in major_seq_train:
    eachset=[]
    for seq in sett:
        for res in seq:
            res=[mapa[i] for i in res]
            for i in range(0,len(res)):
                if i<sp:
                    win=(zeropad*(sp-i))+res[0:ws-(sp-i)]
                    # print(win)
                    # sys.exit()
                    # win=[a for b in win for a in b]
                    eachprotein.append(win)
                elif i>= (len(res)-sp):
                    middle=res[(i-sp):ws-(sp-i)+1]
                    if len(middle) != ws:
                        middle.extend(zeropad*(ws-len(middle)))
                    # middle=[c for d in middle for c in d]
                    eachprotein.append(middle)
                else:
                    wind= res[(i-sp):(i+sp+1)]
                    # win=[a for b in win for a in b]
                    eachprotein.append(wind)
                    # eachprotein=[j for i in eachprotein for j in i]
                # print(eachprotein[-1])
            eachset.append(eachprotein)
            eachprotein=[]
    temptrainingset.append(eachset)
            # eachset=[j for i in eachset for j in i]
#             eachprotein=[]

# print(len(temptrainingset))
# print(len(temptrainingset[0]))
# print(len(temptrainingset[0][0]))
# print(temptrainingset[0][2][-2:])

# print(len(temptrainingset))
# print(len(temptrainingset[0]))
# print(len(temptrainingset[0][0]))
trainingset=[]
for sett in temptrainingset:
    setwhole=[]
    for seq in sett:
        for window in seq:
            # window=[j for i in window for j in i]
            # print(window)
            setwhole.append(window)
    trainingset.append(setwhole)

# print(trainingset)
# sys.exit()

#
# print(len(trainingset))
# print(len(trainingset[0])) #21287
# print(len(trainingset[0][0])) #60

for j in one_seq_test:
    pears=[j for i in one_seq_test for j in i]
# print(len(pears[0][0])) #each residues
# print(len(pears[0])) #each seq.

temptestset=[]


for seq_ in pears:
    # eachsett=[]
    for seqtest in seq_:
        seqtest=[mapa[it] for it in seqtest]
        eachproteintest=[]
        for i in range(0,len(seqtest)):
            if i<sp:
                wis=(zeropad*(sp-i)) + seqtest[0:ws-(sp-i)]
                # wis=[a for b in wis for a in b]
                eachproteintest.append(wis)
            elif i>= (len(seqtest)-sp):
                middle=seqtest[(i-sp):ws-(sp-i)+1]
                if len(middle) != ws:
                    middle.extend(zeropad*(ws-len(middle)))
                # middle=[a for b in middle for a in b]
                eachproteintest.append(middle)
            else:
                w=seqtest[(i-sp):(i+sp+1)]
                # w=[a for b in w for a in b]
                eachproteintest.append(w)
    temptestset.append(eachproteintest)

# print(len(temptestset))
# sys.exit()

arraytrain_list=[]
for everyset in trainingset: #each set
    mappedarray=enc.fit_transform(everyset).toarray()
    arraytrain_list.append(mappedarray)

arraytest_list=[]
for everyset in temptestset:
    mappedarraytest=enc.fit_transform(everyset).toarray()
    arraytest_list.append(mappedarraytest)

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


fea_test_setarray=[]
for x in fea_test_set:
    arrayfea=np.array(x)
    fea_test_setarray.append(arrayfea)

# print(type(fea_train_set))

# print(fea_test_setarray[0:3])
# print(len(fea_test_setarray))
# print(len(arraytest_list))
# print(fea_test_setarray[0])
# print(arraytest_list[0])
# # ##################SVM
#
# # testarrayseq=arraytrain_list[0:3]
# # testarrayfea=fea_train_set[0:3]
# # testlistt=arraytest_list[0:3]
# # fea_test_setarraytest=fea_test_setarray[0:3]
predictionlist=[]
truelist = []
# clf=svm.LinearSVC(class_weight="balanced")
scorelist=[]
p=[]
r=[]
f=[]
sup=[]
clf=svm.SVC(kernel='linear')
codelabel=[0,1,2]
from sklearn.metrics import precision_recall_fscore_support as scoreprfs
# from sklearn.metrics import classification_report
target_names=[]
for i in range(0,len(arraytrain_list)):
    clf.fit(arraytrain_list[i],fea_train_set[i])
    pred=clf.predict(arraytest_list[i])
    predictresult=list(pred)
    # predictionlist.append(predictresult)
    predictionlist.extend(predictresult)
    truelist.extend(fea_test_set[i])
    # standard=clf.score(arraytest_list[i],fea_test_setarray[i])
    # t=standard.std()
    s=float(clf.score(arraytest_list[i],fea_test_setarray[i]))
    scorelist.append(s)
    precision,recall,fscore,support=scoreprfs(fea_test_set[i],pred,labels=codelabel)
    # classification_report(fea_test_set[i],pred)
    # p.append(list(precision))
    # r.append(list(recall))
    f.append(list(fscore))
    # sup.append(list(support))


print ("linear kernel, window=%s" %(ws))
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
# # r_class0=[]
# # r_class1=[]
# # r_class2=[]
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
# #     r_class0.append(r[i][0])
# #     r_class1.append(r[i][1])
# #     r_class2.append(r[i][2])
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
# # r_class0avg=statistics.mean(r_class0)
# # r_class1avg=statistics.mean(r_class1)
# # r_class2avg=statistics.mean(r_class2)
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
# # print (r_class0avg)
# # print (r_class1avg)
# # print (r_class2avg)
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
