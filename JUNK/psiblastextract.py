from sklearn import svm
import numpy as np
from sklearn import datasets

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
    path='/Users/KWL/Documents/project2017sem1/fasta/fasta/PSSM/'+ID+'.fasta.pssm'
    with open(path, 'r') as path:
        path = path.readlines()
        for line in path[4:-6]:
            line = line.split()
            freq=line[22:42]
            # print(freq)
            # print(freq[-1])
            # freq[1][:0]=zeroes
            # freq[-1].extend(zeroes)
            PSIfreq.append(freq)
        PSSMlist.append(PSIfreq)

# print(PSSMlist)
# zeroes='0'*20
# zerolist=[]
# zerolist.extend(zeroes)

zerolist=list(np.zeros(20))

PSIfreqfloat=[]



for eachseq_PSSM in PSSMlist:
    eachseq_PSSM.insert(0,zerolist)
    eachseq_PSSM.append(zerolist)
    seqPSSM=[]
    for eachseq_position in eachseq_PSSM:
        eachseqPSSM=[]
        for eachfreq in eachseq_position:
            eachfreq=(float(eachfreq)/100)
            eachseqPSSM.append(eachfreq)
        seqPSSM.append(eachseqPSSM)
    PSIfreqfloat.append(seqPSSM)



# print(PSIfreqfloat)
#
fealist=[]
for feaindividual in fea:
    fealist.append(feaindividual.split(','))
# #
#
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

# print (len(major_seq_train[0]))
# print (len(major_fea_train[0]))
# print (len(one_seq_test[0]))
# print (len(one_fea_test[0]))


#slidingggggggggggggggggg

# major_seq_train has 4 layers
# 1: residue represented frequency
# 2: seq
# 3: set
test=major_seq_train[0:2]
# print(len(test)) # number of set
# print(len(test[0])) #number of seq in a set
# print(test[0][0]) # residues of 1st seq
# set_train = []

# for eachset in test:
#     settri=[]
#     for eachseq_res in eachset:
#         # seq=[]
#         for j in range(0, len(eachseq_res)):
#             a=[eachseq_res[j:j+3] for j in range(len(eachseq_res)-2)]
#             for eachwindow in a:
#                 b=[item for res in eachwindow for item in res]
#                 #flatten the window (len each: 60)
#                 settri.append(b)

# print(settri)


for j in one_seq_test:
    pears=[j for i in one_seq_test for j in i] #eachseq=a test set(presented as a list)
    seq_test = []
    for eachseqlist in pears:
        for j in range(0, len(eachseqlist)):
            flat=[eachseqlist[j:j+3] for j in range(len(eachseqlist)-2)]
            for each in flat:
                res_in_window=[item for res in each for item in res]
                seq_test.append(res_in_window)

                # print(seq_test)
            # p =[seqres[i:i+3] for i in range(len(eachseqlist)-2)]
            # print(p)
            # seq_test.append(p)
    # print(seq_test)
        # tritest.append(seq_test)

# print (len(seq_test))

# #

# # print (settestmap)
# arraytest_list=[]
# for everyset in settestmap:
#     mappedarraytest=enc.fit_transform(everyset).toarray()
#     arraytest_list.append(mappedarraytest)
#
# # print (len(arraytrain_list[0]))
#
# # ###fea########################################
#
#
# # print (major_fea_train) #1st level set, 2nd level is sequence in the set
# flattenseqinset = []
# for k in major_fea_train: #for every set --> k, flatten
#     for everyseq in k:
#         flatseqinset = [everyseq for i in k for everyseq in i]
#     flattenseqinset.append(flatseqinset)
#
#
#
# # print(flattenseqinset[0])
# dicfeastates={"i":0,"o":1,"l":2,"p":2,"m":2}
# fea_train_set=[]
# for everyset in flattenseqinset:
#     setlist=[]
#     for eachseq in everyset:
#         for eachfea in eachseq:
#             tempoo=dicfeastates[eachfea]
#             setlist.append(tempoo)
#     fea_train_set.append(setlist)
# # print(fea_train_set[0])
# # print (len(fea_train_set[0]))
#
# fea_test_set=[]
# for j in one_fea_test:
#     banana=[j for i in one_fea_test for j in i] #eachseq=a test set(presented as a list)
#     for eachfealist in banana: #eachseq as set
#         fea_test = []
#         for fea in eachfealist:
#             for eachfea in fea:
#                 tempo=dicfeastates[eachfea]
#                 fea_test.append(tempo)
#         fea_test_set.append(fea_test)
#
#
# fea_test_setarray=[]
# for x in fea_test_set:
#     arrayfea=np.array(x)
#     fea_test_setarray.append(arrayfea)
#
#
#
# # print(fea_test_setarray[0:3])
# # print(len(fea_test_setarray))
# # print(len(arraytest_list))
# # print(fea_test_setarray[0])
# # print(arraytest_list[0])
# # ##################SVM
#
# testarrayseq=arraytrain_list[0:3]
# testarrayfea=fea_train_set[0:3]
# testlistt=arraytest_list[0:3]
# fea_test_setarraytest=fea_test_setarray[0:3]
# predictionlist=[]
# truelist = []
# clf=svm.LinearSVC(class_weight="balanced")
# for i in range(0,len(testarrayseq)):
#     clf.fit(arraytrain_list[i],fea_train_set[i])
#     predictresult=list(clf.predict(arraytest_list[i]))
#     # predictionlist.append(predictresult)
#     predictionlist.extend(predictresult)
#     truelist.extend(fea_test_set[i])
#     # accuracy average
#     # -> whole set -> confusion matrix
#
# # print (truelist)
# # print (predictionlist)
# # print(predictionlist[0])
# # predictionlistarray=np.array(predictionlist)
# # print(predictionlistflatten)
# # print(predictionlistflatten[0])
# # testfeatrue=fea_test_setarray[0:3]
# # print(testfeatrue)
# # print(testfeatrue[0])
# # print(predictionlistarray)
# # print(len(predictionlistarray))
# # print(type(predictionlistarray))
#
# # print(fea_test_setarray[0:3])
# # print(len(fea_test_setarray[0:3]))
# # print(type(fea_test_setarray))
#
# import itertools
# import matplotlib.pyplot as plt
# from sklearn.metrics import confusion_matrix
# def plot_confusion_matrix(cm, classes,
#                           normalize=True,
#                           title='Confusion matrix',
#                           cmap=plt.cm.Blues):
#     plt.imshow(cm, interpolation='nearest', cmap=cmap)
#     plt.title(title)
#     plt.colorbar()
#     tick_marks = np.arange(len(classes))
#     plt.xticks(tick_marks, classes, rotation=45)
#     plt.yticks(tick_marks, classes)
#
#     if normalize:
#         cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
#         print("Normalized confusion matrix")
#     else:
#         print('Confusion matrix, without normalization')
#
#     print(cm)
#
#     thresh = cm.max() / 2.
#     for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
#         plt.text(j, i, cm[i, j],
#                  horizontalalignment="center",
#                  color="white" if cm[i, j] > thresh else "black")
#
#     plt.tight_layout()
#     plt.ylabel('True label')
#     plt.xlabel('Predicted label')
#
# # Compute confusion matrix
# codelabel=[0,1,2]
# cnf_matrix = confusion_matrix(truelist, predictionlist,labels=codelabel)
# np.set_printoptions(precision=2)
#
# # Plot normalized confusion matrix
# plt.figure()
# plot_confusion_matrix(cnf_matrix, classes=codelabel, normalize=True,
#                       title='Normalized confusion matrix,window size=3')
#
#
# plt.show()
#
# from sklearn.utils import compute_class_weights
