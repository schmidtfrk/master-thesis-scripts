#checks all possible combinations of products of golay-code stabilizers and looks at their weigth
data=[[0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
[1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
[0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
[1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
[1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
others=[]
zwoelf=0
for q1 in range(2):
    for q2 in range(2):
        for q3 in range(2):
            for q4 in range(2):
                for q5 in range(2):
                    for q6 in range(2):
                        for q7 in range(2):
                            for q8 in range(2):
                                for q9 in range(2):
                                    for q10 in range(2):
                                        for q11 in range(2):
                                            test=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                                            for i in range(len(data[0])):
                                                test[i]=(test[i]+q1*data[0][i]+q2*data[1][i]+q3*data[2][i]+q4*data[3][i]+q5*data[4][i]+q6*data[5][i]+q7*data[6][i]+q8*data[7][i]+q9*data[8][i]+q10*data[9][i]+q11*data[10][i])%2
                                            check=0
                                            for k in range(len(test)):
                                                check+=test[k]
                                            if check not in others:
                                                others.append(check)
                                            if check ==12:
                                                zwoelf+=1
print(others)
print(zwoelf)
#others is a list of all possible weigths of stabilizer combinations                                                                                                        
    
