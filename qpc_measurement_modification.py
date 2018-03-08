from bsm_lib import *

#find the maximum logical Bell measurement efficiency for the QPC(2,2)
#in the no loss case
print('QPC(2,2)')
maximum=0
for info in itertools.product('xyz',repeat=4):
    data=general_safe_bm([[1,1,0,0],[1,1,1,1]],[[1,1,0,0],[0,0,1,1]],[1,0,1,0],#z-type input
                         [[1,0,1,0],[1,1,0,0],[0,0,1,1]],[[1,1,1,1]],[1,1,0,0],#x-type input
                         info)
    if data>maximum:
        maximum=data
print('maximum number of decodable codeword combinations'+str(maximum))

count=0
info2=[]
for info in itertools.product('xyz',repeat=4):
    data=general_safe_bm([[1,1,0,0],[1,1,1,1]],[[1,1,0,0],[0,0,1,1]],[1,0,1,0],#z-type input
                         [[1,0,1,0],[1,1,0,0],[0,0,1,1]],[[1,1,1,1]],[1,1,0,0],#x-type input
                         info)
    if data==maximum:
        info2.append(info)
        count+=1
print('number of information formations which give the maximum no los efficiency:'+str(count))

for info in info2:
    print('efficiency of'+str(general_safe_bm_erasure([[1,1,0,0],[1,1,1,1]],[[1,1,0,0],[0,0,1,1]],[1,0,1,0],#z-type input
                         [[1,0,1,0],[1,1,0,0],[0,0,1,1]],[[1,1,1,1]],[1,1,0,0],#x-type input
                         info))+' with an information formation of '+str(info))
print('__________________')
print(general_safe_bm_erasure([[1,1,0,0],[1,1,1,1]],[[1,1,0,0],[0,0,1,1]],[1,0,1,0],#z-type input
                         [[1,0,1,0],[1,1,0,0],[0,0,1,1]],[[1,1,1,1]],[1,1,0,0],#x-type input
                         ['y','x','z','z']))
print(general_safe_bm_erasure([[1,1,0,0],[1,1,1,1]],[[1,1,0,0],[0,0,1,1]],[1,0,1,0],#z-type input
                         [[1,0,1,0],[1,1,0,0],[0,0,1,1]],[[1,1,1,1]],[1,1,0,0],#x-type input
                         ['z','z','z','z']))
print('__________________')
#testing this safe information formation for QPC(3,2)
print('QPC(3,2)')
print('new formation')
print(general_safe_bm_erasure([[1,1,1,1,0,0],[0,0,1,1,1,1],[1,1,0,0,0,0]],[[1,1,0,0,0,0],[0,0,1,1,0,0],[0,0,0,0,1,1]],[1,0,1,0,1,0],
                [[1,1,0,0,0,0],[0,0,1,1,0,0],[0,0,0,0,1,1],[1,0,1,0,1,0]],[[1,1,1,1,0,0],[0,0,1,1,1,1]],[1,1,0,0,0,0],
                ['y','x','z','z','z','z']))
print('old only ZZ=1 formation')
print(general_safe_bm_erasure([[1,1,1,1,0,0],[0,0,1,1,1,1],[1,1,0,0,0,0]],[[1,1,0,0,0,0],[0,0,1,1,0,0],[0,0,0,0,1,1]],[1,0,1,0,1,0],
                [[1,1,0,0,0,0],[0,0,1,1,0,0],[0,0,0,0,1,1],[1,0,1,0,1,0]],[[1,1,1,1,0,0],[0,0,1,1,1,1]],[1,1,0,0,0,0],
                ['z','z','z','z','z','z'])) 
