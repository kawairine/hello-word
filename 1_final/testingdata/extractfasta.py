# openfile=open('50proteinstest.txt','r')
# openfile=openfile.readlines()
#
# idlist = []
# seqlist=[]
# fea=[]
# for i in range(0, len(openfile), 2):
#     idlist.append(openfile[i].rstrip('\n'))
#     seqlist.append(openfile[i+1].rstrip('\n'))
#     a = openfile[i].rstrip('\n')
#     filename = 'fasta/'+str(a)+'.fasta'
#     with open(filename, 'w+') as apples:
#         apples.write(str(openfile[i]))
#         apples.write(str(openfile[i+1].rstrip('\n')))


openfile=open('50proteinstest.txt','r')
openfile=openfile.readlines()


for i in range(0, len(openfile),3):
    a = openfile[i].rstrip('\n').lstrip('>')
    filename = 'fasta/'+str(a)+'.fasta'
    with open(filename, 'w+') as apples:
        apples.write(str(openfile[i]))
        apples.write(str(openfile[i+1].rstrip('\n')))
