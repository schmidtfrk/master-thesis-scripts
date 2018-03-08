import itertools
import copy
from math import factorial as fact

#warning: carefully read the documentation wheter a function counts succesfull or failed logical Bell measurements
#auxillary functions used in the scripts calculating efficiencies
def binomial(n, k):
    """ Calculates the binomial coefficient n over k
    """
    try:
        binom = fact(n) // fact(k) // fact(n - k)
    except ValueError:
        binom = 0
    return binom

def optical_condition(codeword,conditions,informationplaces=0,addinfo=0):
    ''' Determines if the given codition is meet.

        Input:
        codeword: list of 0 or 1s which correspond to possible
                   measurement outcomes
        conditions: list of 0 or 1s which show which information is needed
                    for logical information
        informationplaces:list of 0 or 1s, determines the X/Z type of
                          of the needed measurement for additional info
                          (0 means safe ZZ)
        addinfo: list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        Boolean
    '''
    #Function was modifiedin a ad-hoc manner in order to obtain compatibility
    #with the old function bsmprobability, but it wasn't well tested.
    if addinfo == 0:
        addinfo=[1]*len(codeword)
    if informationplaces == 0:#needed for compatibility with old function
        informationplaces=[0]*len(codeword)#bsmprobability, not well tested
    for i in range(len(conditions)):
        localboolean=1
        for j in range(len(codeword)):
            if conditions[i][j]==1:
                if informationplaces[j]!=1:
                    if codeword[j]!=addinfo[j] :
                        localboolean=0
                        break
        if localboolean==1:
            return True
    if localboolean==0:
        return False

def create_logicals(stabilizer,log):
    """ Creates all logical operators if one represantive is given. """
    liste=list(span_generator_mod2(stabilizer))
    if liste!=[]:
        for i in liste:
            for j in range(len(i)):
                i[j]+=log[j]
                i[j]%=2
    else:
        liste.append(log)
    return liste

def create_erasure_patterns(numberofqubits):
    """ Creates a list of all possible erasure patterns with numberofqubits
        qubits, which is more efficient than
        list(span_generator_mod2(ones(numberofqubits)))
    """
    y=[[int(x) for x in list('{0:0b}'.format(i))]for i in range(2**numberofqubits)]
    length=len(y[-1])
    for i in y:
        while len(i)<length:
            i.insert(0,0)
    return y

def span_generator_mod2(liste):
    """ Takes a list of vectors and returns all linearcombinations(modulo 2)

        Input:
        liste: list of lists of numbers

        Output:
        generator of Lists of numbers
    """
    if liste==[]:
        return liste
    n = len(liste[0])
    transpose = list(zip(*liste))
    coefficients = itertools.product(range(2), repeat=len(liste))
    for coeff in coefficients:
        yield [sum((a * c) for a, c in zip(transpose[i], coeff)) % 2 for i in range(n)]

def ones(n):
    """ Creates a list of canonical basis vectors of R^n

        Input:
        n: integer number

        Output:
        results: List of lists

    """
    results=[]
    for i in range(n):
        zeroes=[0]*n
        zeroes[i]=1
        results.append(zeroes)
    return results

def ones2(basis):
    """ Input:
        basis: List of Integers (0 or 1)
        Output:
        List of lists of Integers (each list has one 1 on a place like in basis)
    """
    results=[]
    for i in range(len(basis)):
        zeroes=[0]*len(basis)
        if basis[i]==1:
            zeroes[i]=1
            results.append(zeroes)
    return results

#scripts for calculating the logical Bell efficiency with ideal Bell measurements
def erasure(xstabilizer,xlog,zstabilizer,zlog):
    """ Counts uncorrectable losses.
        It assumes a [n,1,d]-CSS(C_X,C_Z)-code.

        Input:
        xstablizer: list of lists of the support of the X-stabilizer
                    generators (0 ist qubit not in support ,1 otherwise)
        xlog: list of the support of one representative of logical X-operators
        zstablizer: list of lists of the support of the Z-stabilizer
                    generators (0 ist qubit not in support ,1 otherwise)
        zlog: list of the support of one representative of logical Z-operators

        Output:
        dummy: List of integers counting uncorrectable errors
    """
    xlogical=create_logicals(xstabilizer,xlog)
    zlogical=create_logicals(zstabilizer,zlog)
    n=len(xstabilizer[0])
    dummy=[0]*(n+1)
    xerrors=[]
    errors=create_erasure_patterns(n)
    for i in range(len(errors)):
        breakvar=0
        for j in range(len(xlogical)):
            breakvarx=0
            for k in range(n):
                if 1==xlogical[j][k]:
                    if 1!=errors[i][k]:
                        breakvarx=1
            if 0==breakvarx:
                xerrors.append(errors[i])
                breakvar=1
                break
        if breakvar==0:
            for j in range(len(zlogical)):
                breakvarz=0
                for k in range(n):
                    if 1==zlogical[j][k]:
                        if 1!=errors[i][k]:
                            breakvarz=1
                if 0==breakvarz:
                    xerrors.append(errors[i])
                    break
    for i in range(len(xerrors)):
        dummy[sum(xerrors[i])]+=1
    return dummy

def erasure_color(xstabilizer,xlog):
    """ Counts uncorrectable losses.
        It assumes a [n,1,d]-CSS(C_X,C_X)-code.

        Input:
        xstablizer: list of lists of the support of the stabilizer
                    generators (0 ist qubit not in support ,1 otherwise)
        xlog: list of the support of one representative of logical operators

        Output:
        dummy: List of integers counting uncorrectable errors
    """
    xlogical=create_logicals(xstabilizer,xlog)
    n=len(xstabilizer[0])
    dummy=[0]*(n+1)
    xerrors=[]
    errors=create_erasure_patterns(n)
    for i in range(len(errors)):
        breakvar=0
        for j in range(len(xlogical)):
            breakvarx=0
            for k in range(n):
                if 1==xlogical[j][k]:
                    if 1!=errors[i][k]:
                        breakvarx=1
            if 0==breakvarx:
                xerrors.append(errors[i])
                breakvar=1
                break
    for i in range(len(xerrors)):
        dummy[sum(xerrors[i])]+=1
    return dummy

def erasure_color_fast(xstabilizer,xlog):
    """ Counts correctable losses, such that all losses can be threated
        with the usage of algebra. It assumes a [n,1,d]-CSS(C_X,C_X)-code.

        Input:
        xstablizer: list of lists of the support of the stabilizer
                    generators (0 ist qubit not in support ,1 otherwise)
        xlog: list of the support of one representative of logical operators

        Output:
        dummy: List of integers counting correctable errors
    """
    xlogical=create_logicals(xstabilizer,xlog)
    n=len(xstabilizer[0])
    dummy=[0]*(n+1)
    errors=create_erasure_patterns(n)
    for error in errors:
        if sum(error)>(len(xlog)-1)/2:
            continue
        breakvar=0
        for xlog in xlogical:
            breakvarx=0
            for k in range(n):
                if 1==xlog[k]:
                    if 1==error[k]:
                        breakvarx=1
                        break
            if 0==breakvarx:
                dummy[sum(error)]+=1
                break
    for i in range((n+1)//2,n):#using theorem 2
        dummy[i]=binomial(n,i)-dummy[n-i]        
    return dummy

#logical Bell measurement efficiency with linear optical Bell measurements
def general_safe_bm(zcode,zstabgen,zloggen,xcode,xstabgen,xloggen,infos,addinfo=0):
    """ Counts the number of possible succesfull codeword combinations of a
        logical Bell measurement without photon loss and arbitrary safe information
        formations are can be used

        Input:
        zcode: list of lists of a basis of codewords of the code C_Z
        zstabgen: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zloggen: list of of support of one representative of logical Z
        xcode: list of lists of a basis of codewords of the code C_X
        xstabgen: list of lists of support of X-stabilizers generators
        xloggen: list of of support of one representative of logical X
        infos: list of 'x','y','z', each individual string denotes the
                safe information of this safe information Bell measurement
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        count: number of possible succesfull logical BSM measurements
    """
    n=len(zcode[0])
    if addinfo==0:
        addinfo=[1]*n
    count=0
    xcodewords=list(span_generator_mod2(xcode))
    zcodewords=list(span_generator_mod2(zcode))
    xlog=create_logicals(xstabgen,xloggen)
    zlog=create_logicals(zstabgen,zloggen)
    ylog=[]
    for i in xlog:
        for j in zlog:
            ylog.append([i,j])
    for xword in xcodewords:
        for zword in zcodewords:
            xinfo=[0]*n
            zinfo=[0]*n
            yinfo=[0]*n
            for index,k in enumerate(infos):
                if k == 'x':
                    xinfo[index]=1
                    if xword[index]==addinfo[index]:
                        zinfo[index]=1
                        yinfo[index]=1
                if k== 'z':
                    zinfo[index]=1
                    if zword[index]==addinfo[index]:
                        xinfo[index]=1
                        yinfo[index]=1
                if k =='y':
                    yinfo[index]=1
                    if (xword[index]+zword[index])%2==addinfo[index]:
                        xinfo[index]=1
                        zinfo[index]=1
            for i in xlog:
                xboolean=True
                for index,k in enumerate(i):
                    if k==1 and xinfo[index]==0:
                        xboolean=False
                        break
                if xboolean:
                    break
            for i in zlog:
                zboolean=True
                for index,k in enumerate(i):
                    if k==1 and zinfo[index]==0:
                        zboolean=False
                        break
                if zboolean:
                    break
            if xboolean and zboolean:
                count+=1
                continue
            for i in ylog:
                yboolean=True
                for index in range(n):
                    if ((i[0][index]==1 and i[1][index]==0 and xinfo[index]==0) or
                        (i[1][index]==1 and i[0][index]==0 and zinfo[index]==0) or
                        (i[0][index]==1 and i[1][index]==1 and yinfo[index]==0)):
                        yboolean=False
                        break
                if yboolean:
                    break
            if (xboolean and yboolean) or (zboolean and yboolean) :
                count+=1
                continue
    return count                   


def general_safe_bm_erasure(zcode,zstabgen,zloggen,xcode,xstabgen,xloggen,infos,addinfo=0):
    """ Counts the number of possible succesfull codeword combinations of a
        logical Bell measurement with photon loss and arbitrary safe information
        formations are can be used

        Input:
        zcode: list of lists of a basis of codewords of the code C_Z
        zstabgen: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zloggen: list of of support of one representative of logical Z
        xcode: list of lists of a basis of codewords of the code C_X
        xstabgen: list of lists of support of X-stabilizers generators
        xloggen: list of of support of one representative of logical X
        infos: list of 'x','y','z', each individual string denotes the
                safe information of this safe information Bell measurement
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        count: list of numbers of possible succesfull logical BSM
               measurements (0'th element gives no loss cases,
               1'st elements gives the number of all valid logical BSM
               measurements with 1 lost photon...)

    """
    n=len(xloggen)
    if addinfo==0:
        addinfo=[1]*n
    losses=create_erasure_patterns(n)
    count=[0]*(n+1)
    xcodewords=list(span_generator_mod2(xcode))
    zcodewords=list(span_generator_mod2(zcode))
    xlog=create_logicals(xstabgen,xloggen)
    zlog=create_logicals(zstabgen,zloggen)
    ylog=[]
    for i in xlog:
        for j in zlog:
            ylog.append([i,j])
    for errors in losses:
        for xword in xcodewords:
            for zword in zcodewords:
                xinfo=[0]*n
                zinfo=[0]*n
                yinfo=[0]*n
                for index,k in enumerate(infos):
                    if k == 'x' and errors[index]==0:
                        xinfo[index]=1
                        if xword[index]==addinfo[index]:
                            zinfo[index]=1
                            yinfo[index]=1
                    if k== 'z' and errors[index]==0:
                        zinfo[index]=1
                        if zword[index]==addinfo[index]:
                            xinfo[index]=1
                            yinfo[index]=1
                    if k =='y' and errors[index]==0:
                        yinfo[index]=1
                        if (xword[index]+zword[index])%2==addinfo[index]:
                            xinfo[index]=1
                            zinfo[index]=1
                for i in xlog:
                    xboolean=True
                    for index,k in enumerate(i):
                        if k==1 and xinfo[index]==0:
                            xboolean=False
                            break
                    if xboolean:
                        break
                for i in zlog:
                    zboolean=True
                    for index,k in enumerate(i):
                        if k==1 and zinfo[index]==0:
                            zboolean=False
                            break
                    if zboolean:
                        break
                if xboolean and zboolean:
                    count[sum(errors)]+=1
                    continue
                for i in ylog:
                    yboolean=True
                    for index in range(n):
                        if ((i[0][index]==1 and i[1][index]==0 and xinfo[index]==0) or
                            (i[1][index]==1 and i[0][index]==0 and zinfo[index]==0) or
                            (i[0][index]==1 and i[1][index]==1 and yinfo[index]==0)):
                            yboolean=False
                            break
                    if yboolean:
                        break
                if (xboolean and yboolean) or (zboolean and yboolean) :
                    count[sum(errors)]+=1
                    continue
    return count                    


def general_safe_bm_erasure_advanced(zcode,zstabgen,zloggen,xcode,xstabgen,xloggen,infos,addinfo=0):
    """ Counts the number of possible succesfull codeword combinations of a
        logical Bell measurement with photon loss and arbitrary safe information
        formations are can be used and advanced Bell emasurements are allowed.

        Input:
        zcode: list of lists of a basis of codewords of the code C_Z
        zstabgen: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zloggen: list of of support of one representative of logical Z
        xcode: list of lists of a basis of codewords of the code C_X
        xstabgen: list of lists of support of X-stabilizers generators
        xloggen: list of of support of one representative of logical X
        infos: list of 'x','y','z', each individual string denotes the
                safe information of this safe information Bell measurement
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        count: list of list of list of  integers.
                count[j] consists of of a list of list of three integers and
                describes events with j erasures
                count[j][k][0] gives the number of measurements where
                additional information was given.(k is only an index counting
                these triple integers)
                count[j][k][1] gives the number of measurements where additional
                information could have been given.
                count[j][k][2] gives the number of succesful identifications
                of the logical Bell state.

    """
    n=len(xloggen)
    if addinfo==0:
        addinfo=[1]*n
    losses=create_erasure_patterns(n)
    count=[]
    [count.append([]) for i in range(n+1)]
    xcodewords=list(span_generator_mod2(xcode))
    zcodewords=list(span_generator_mod2(zcode))
    xlog=create_logicals(xstabgen,xloggen)
    zlog=create_logicals(zstabgen,zloggen)
    ylog=[]
    for i in xlog:
        for j in zlog:
            ylog.append([i,j])
    for errors in losses:
        for xword in xcodewords:
            for zword in zcodewords:
                foo=[0]*n
                xinfo=[0]*n
                zinfo=[0]*n
                yinfo=[0]*n
                for index,k in enumerate(infos):
                    if k == 'x' and errors[index]==0:
                        xinfo[index]=1
                        if xword[index]==addinfo[index]:
                            zinfo[index]=1
                            yinfo[index]=1
                    if k== 'z' and errors[index]==0:
                        zinfo[index]=1
                        if zword[index]==addinfo[index]:
                            xinfo[index]=1
                            yinfo[index]=1
                    if k =='y' and errors[index]==0:
                        yinfo[index]=1
                        if (xword[index]+zword[index])%2==addinfo[index]:
                            xinfo[index]=1
                            zinfo[index]=1
                for i in range(n):
                    if errors[i]==0 and (xinfo[i]==0 or yinfo[i]==0 or zinfo[i]==0):
                        foo[i]=1
                neededaddinfo=list(span_generator_mod2(ones2(foo)))
                if len(neededaddinfo)== 0:
                    neededaddinfo=[foo]
                for info in neededaddinfo:
                    xinfo2=xinfo[:]
                    yinfo2=yinfo[:]
                    zinfo2=zinfo[:]
                    for i in range(n):
                        if info[i]==1:
                            xinfo2[i]=1
                            yinfo2[i]=1
                            zinfo2[i]=1
                    for i in xlog:
                        xboolean=True
                        for index,k in enumerate(i):
                            if k==1 and xinfo2[index]==0:
                                xboolean=False
                                break
                        if xboolean:
                            break
                    for i in zlog:
                        zboolean=True
                        for index,k in enumerate(i):
                            if k==1 and zinfo2[index]==0:
                                zboolean=False
                                break
                        if zboolean:
                            break
                    if xboolean and zboolean:
                        toappend=True
                        for g in count[sum(errors)]:
                            if g[0]==sum(info) and g[1]==sum(foo):
                                g[2]+=1
                                toappend=False
                                break
                        if toappend:
                            count[sum(errors)].append([sum(info), sum(foo),1])
                        continue
                    for i in ylog:
                        yboolean=True
                        for index in range(n):
                            if ((i[0][index]==1 and i[1][index]==0 and xinfo2[index]==0) or
                                (i[1][index]==1 and i[0][index]==0 and zinfo2[index]==0) or
                                (i[0][index]==1 and i[1][index]==1 and yinfo2[index]==0)):
                                yboolean=False
                                break
                        if yboolean:
                            break
                    if (xboolean and yboolean) or (zboolean and yboolean) :
                        toappend=True
                        for g in count[sum(errors)]:
                            if g[0]==sum(info) and g[1]==sum(foo):
                                g[2]+=1
                                toappend=False
                                break
                        if toappend:
                            count[sum(errors)].append([sum(info), sum(foo),1])
    return count                  

#old functions for caclulating logical Bell measurement efficiencies
#which lack some features


def bsm_prob_mix(zcodewords,zstabilizer,zlog,xcodewords,xstabilizer,xlog,xinformation,addinfo=0):
    """ Counts the number of possible succesfull logical Bell measurement
        without photon loss.

        Input:
        zcodewords: list of lists of a basis of codewords of the code C_Z
        zstabilizers: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zlog: list of of support of one representative of logical Z
        xcodewords: list of lists of a basis of codewords of the code C_X
        xstabilizers: list of lists of support of X-stabilizers generators
        xlog: list of of support of one representative of logical X
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        count: number of possible succesfull logical BSM measurements

    """
    n=len(zcodewords[0])
    zinformation=[]
    count=0
    for i in range(n):
        zinformation.append((xinformation[i]+1)%2)
    zliste=create_logicals(zstabilizer,zlog)
    xliste=create_logicals(xstabilizer,xlog)
    zcodewordslist=list(span_generator_mod2(zcodewords))
    xcodewordslist=list(span_generator_mod2(xcodewords))
    for xword in xcodewordslist:
        for zword in zcodewordslist:
            if optical_condition(zword,xliste,xinformation,addinfo):
                if optical_condition(xword,zliste,zinformation,addinfo):
                    count+=1
    return count

def bsm_mix_erasure(zcodewords,zstabilizer,zlog,xcodewords,xstabilizer,xlog,xinformation,addinfo=0):
    """ Counts the number of possible succesfull logical Bell measurement
        with photon loss.

        Input:
        zcodewords: list of lists of a basis of codewords of the code C_Z
        zstabilizers: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zlog: list of of support of one representative of logical Z
        xcodewords: list of lists of a basis of codewords of the code C_X
        xstabilizers: list of lists of support of X-stabilizers generators
        xlog: list of of support of one representative of logical X
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        count: list of numbers of possible succesfull logical BSM
               measurements (0'th element gives no loss cases,
               1'st elements gives the number of all valid logical BSM
               measurements with 1 lost photon...)
    """
    n=len(zcodewords[0])
    zinformation=[]
    for i in range(n):
        zinformation.append((xinformation[i]+1)%2)
    zliste=create_logicals(zstabilizer,zlog)
    xliste=create_logicals(xstabilizer,xlog)
    count=[0]*(n+1)
    zcodewordslist=list(span_generator_mod2(zcodewords))
    xcodewordslist=list(span_generator_mod2(xcodewords))
    errors=create_erasure_patterns(n)
    for error in errors:
        xwordlist_copy=copy.deepcopy(xcodewordslist)
        zwordlist_copy=copy.deepcopy(zcodewordslist)
        zinformation_copy=zinformation[:]#get a workcopy of z and xinformation 
        xinformation_copy=xinformation[:]                                        
        for i in range(n):                                                      
            if error[i]==1:                                                     
                zinformation_copy[i]=0                                          
                xinformation_copy[i]=0                                          
        for xword in xwordlist_copy:
            for zword in zwordlist_copy:
                for i in range(n):
                    if error[i]==1:
                        xword[i]=0
                        zword[i]=0
                if optical_condition(zword,xliste,xinformation_copy,addinfo):
                    if optical_condition(xword,zliste,zinformation_copy,addinfo):
                        count[sum(error)]+=1
    return count


def bsm_mix_erasure_advanced(xcodeword, xstabilizer, xlog, zcodeword, zstabilizer, zlog, xinformation, addinfo=0):#todo implement addinfo
    """ Counts the number of possible succesfull logical Bell measurement
        with photon loss and allows advanced bsms.

        Input:
        zcodewords: list of lists of a basis of codewords of the code C_Z
        zstabilizers: list of lists of support of Z-stabilizers generators
                      (0 if qubit not in support, 1 if qubit in support)
        zlog: list of of support of one representative of logical Z
        xcodewords: list of lists of a basis of codewords of the code C_X
        xstabilizers: list of lists of support of X-stabilizers generators
        xlog: list of of support of one representative of logical X
        addinfo: (opt) list of 0 or 1s, determines the needed measurement
                outcome for add. info

        Output:
        results[a][b] gives the number of logical measurement successes
        with a lost photons and b photons where additional imformation
        was gained by the use of advanced bsms.
        
    """
    n=len(xcodeword[0])
    if addinfo==0:
        addinfo=[1]*n
    zinformation=[(xinformation[i]+1)%2 for i in range(n)]
    xlogical=create_logicals(xstabilizer,xlog)
    zlogical=create_logicals(zstabilizer,zlog)
    errors=create_erasure_patterns(n)
    results=[]
    [results.append([]) for i in range(n+1)]
    xwords=list(span_generator_mod2(xcodeword))
    zwords=list(span_generator_mod2(zcodeword))
    for error in errors:
        xwords_copy=copy.deepcopy(xwords)
        zwords_copy=copy.deepcopy(zwords)
        for xword in xwords_copy:
            for zword in zwords_copy:
                foo=[0]*n
                for i in range(n):
                    if error[i]==0==xinformation[i] and zword[i]!=addinfo[i]:
                        foo[i]=1
                    if error[i]==0==zinformation[i] and xword[i]!=addinfo[i]:
                        foo[i]=1
                for i in range(n):
                    if error[i]==1:
                        xword[i]=(addinfo[i]+1)%2
                        zword[i]=(addinfo[i]+1)%2
                neededaddinfo=list(span_generator_mod2(ones2(foo)))
                if  not len(neededaddinfo)>=1:
                    neededaddinfo=[foo]
                for info in neededaddinfo:
                    xinformation_copy=xinformation[:]
                    zinformation_copy=zinformation[:]
                    number_of_possible_info=sum(foo)
                    succ=sum(info)
                    for j in range(n):
                        if info[j]==1:
                            xinformation_copy[j]=1
                            zinformation_copy[j]=1
                        if error[j]==1:
                            xinformation_copy[j]=0
                            zinformation_copy[j]=0
                    if optical_condition(zword, xlogical, xinformation_copy,addinfo):
                        if optical_condition(xword, zlogical, zinformation_copy,addinfo):
                            toappend=True
                            for g in results[sum(error)]:
                                if g[0]==succ and g[1]==number_of_possible_info:
                                    g[2]+=1
                                    toappend=False
                                    break
                            if toappend:
                                results[sum(error)].append([succ, number_of_possible_info,1])
    return results

def bsmprobability(kernel,stabilizer,log):
    """ Counts the number of Z codewords which allow a identification of
    logical Bell states. Since it ignores X-codewords, it's faster than using
    general_safe_bm for this task.
    
    kernel: list of lists which are a basis of all Z-codewords
    stabilizer: list of list whit the support of X-stabilizer generators
    log: list with the support of logical X 
    """
    allowedbsm=list(span_generator_mod2(kernel))
    logicalcond=list(span_generator_mod2(stabilizer))
    for i in range(len(logicalcond)):
        for j in range(len(logicalcond[0])):
            logicalcond[i][j]=(logicalcond[i][j]+log[j])%2
    count=0
    for i in range(len(allowedbsm)):
        if optical_condition(allowedbsm[i],logicalcond):
            count+=1
    print(str(count)+' of '+str(2**len(kernel))+' possible BSM results(*2^(n-number of x stabilizer)) can be detected with linear optics')
