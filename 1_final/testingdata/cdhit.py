openfile=open('listopm.txt.full_topology','r')
openfile=openfile.readlines()


with open('CDwithall.txt', 'w') as fwrite:
    for i in range(0, len(openfile), 3):
            pos=openfile[i].find("|")
            tmp=openfile[i][:pos]
            fwrite.write(tmp+'\n')
            fwrite.write(openfile[i+1])

with open('CDwithallplustop.txt', 'w') as fwrite:
    for i in range(0, len(openfile), 3):
            pos=openfile[i].find("|")
            tmp=openfile[i][:pos]
            fwrite.write(tmp+'\n')
            fwrite.write(openfile[i+1])
            fwrite.write(openfile[i+2])


xfile=open('CDwithallplustop.txt','r')
yfile=open('30idcdhitresult.txt','r')

seqandtop=xfile.readlines()
thirty=yfile.readlines()


count=0
with open('50proteinstest.txt', 'w') as pwrite:
    for i in range(0,len(thirty)):
        if ">" in thirty[i]:
            tmpp=thirty[i]
            for i in range(0,len(seqandtop)):
                if count <= 49:
                    if tmpp in seqandtop[i]:
                            pwrite.write(seqandtop[i])
                            pwrite.write(seqandtop[i+1])
                            pwrite.write(seqandtop[i+2])
                            count+=1

#
# counn=0
# with open('50proteinstest.txt', 'r') as test:
#     x=test.readlines()
#     for i in range(0,len(x)):
#         if ">" in x[i]:
#             counn+=1
#             print(counn)
