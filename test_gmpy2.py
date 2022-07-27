import scipy.signal
import gmpy2
import numpy as np

def smooth(R):
    S = [];
    S.append(R[0])
    S.append(R[1])
    for i in range(2,len(R)-2):
        S.append((R[i-2]+R[i-1]+R[i]+R[i+1]+R[i+2])/5)
    S.append(R[len(R)-2])
    S.append(R[len(R)-1])
    return S



def F_to_C(RT,df2):
    # 数据初始化
    deltaz = []
    time_z =[]
    R = []
    C = []
    den_end = []
    num_end = []
    RZ =[]
    Cauer_C = []
    Cauer_R = []
    dz = []
    gmpy2.get_context().precision = 200
    # 导入数据
    # pd.set_option('display.unicode.east_asian_width',True)
    # df = pd.read_excel(r'C:\Users\SmartZZH\Desktop\T3ster\RT.xlsx',dtype=str)     # 读取excel文件
    # df2 = pd.read_excel(r'C:\Users\SmartZZH\Desktop\T3ster\time_z.xlsx',dtype=str)     # 读取excel文件python coeffs
    # RT = pd.DataFrame(df)   #贝叶斯计算完毕后的数据
    # 高精度转换
    for i in range(len(df2)-1):
        RZ.append(RT[i])
        RZ[i] = gmpy2.mpfr(RZ[i])
    # 将t转换为z

    for i in range(len(df2)-1):
        dz.append(gmpy2.mpfr(df2[i]))

    # 数据离散化获取Δz
    length = len(dz)

    deltaz = (dz[length -1] - dz[0])/length

    for i in range(length):
        R.append(RZ[i] * deltaz) # R = R(zi) * Δz
        C.append(gmpy2.exp(dz[i])/R[i])
    # print(R)
    #分母系数
    den_end = [gmpy2.mpfr('1')]
    for i in range(length):
        den_end = scipy.signal.convolve(den_end,[gmpy2.exp(dz[i]),gmpy2.mpfr('1')])
    num_end = gmpy2.mpfr('0')
    for i in range(length):
        w = [R[i]]
        for j in range(length):
            if i != j:
                w = scipy.signal.convolve(w,[gmpy2.exp(dz[j]),gmpy2.mpfr('1')])
        num_end = num_end + w
    D = np.copy(den_end)
    N = np.copy(num_end)
    # 计算2
    for i in range(length-5):
        n = np.hstack([D, N]).max()
        D = D / n
        N = N / n
        Cauer_C = np.append(Cauer_C,D[0] / N[0])        # Consecutive C's in the chain
        Nnew = N * N[0]
        Dnew = D * N[0] - np.hstack([N, [gmpy2.mpfr("0")]]) * D[0]
        Dnew = Dnew[1:]     # in MATLAB: Dnew(1) = []
        Cauer_R = np.append(Cauer_R,Nnew[0] / Dnew[0])  # Consecutive R's in the chain
        N = Nnew * Dnew[0] - Dnew * Nnew[0]
        N = N[1:]           # in MATLAB: N(1) = []
        D = Dnew * Dnew[0]

    for i in range(1,len(Cauer_R)):
        Cauer_C[i] = Cauer_C[i-1] + Cauer_C[i]
        Cauer_R[i] = Cauer_R[i-1] + Cauer_R[i]
    Cauer_C = np.hstack((np.array([0]),Cauer_C))
    Cauer_R = np.hstack((np.array([0]),Cauer_R))
    Cauer_C = smooth(Cauer_C)
    Cauer_C = smooth(Cauer_C)
    Cauer_R = smooth(Cauer_R)
    Cauer_R = smooth(Cauer_R)
    return np.vstack((Cauer_R[0:length - 50],Cauer_C[0:length - 50])).astype(float).tolist()
# print(Cauer_C)
# print(Cauer_R)
# data = {'R':Cauer_R,"C":Cauer_C}
# dataexcel = pd.DataFrame.from_dict(data,orient='index')
# writer = pd.ExcelWriter(r'C:\Users\SmartZZH\Desktop\T3ster\data.xlsx')
# dataexcel.to_excel(writer)
# writer.save()
