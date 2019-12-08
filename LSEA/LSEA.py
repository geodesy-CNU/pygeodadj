# 最小二乘法相關程序
import numpy as np
from numpy.linalg import inv
import math

# General purpose matrix printing function
def Matrix_print(X,matrix_title="Matrix"):
    print("\n%s"%matrix_title)
    
    if len(X.shape)==1:
        for x in X:
            print('%10.4f'%(x))
    else:
        print(X)
        
# General purpose matrix printing function for two vectors
def Matrix_print2(X,sX,matrix_title="Matrix",unit=""):
    print("\n%s"%matrix_title)
    
    for x,sx in zip(X,sX):
        if unit is None:
            print('%.4f +/- %.4f'%(x,sx))
        else:
            print('%.4f%s +/- %.4f (%s)'%(x, unit, sx, unit))
            

def LSEA(A,L, p=None, px=None, unit=None):
    m,n=A.shape
    #print(m,n)

    #% 若未賦予權值，以等權計算
    if p is None:
        P=np.diag(np.ones(m),0)
    else:
        P=np.diag(p)  

    #% 未知數的答解 
    #% 法方程式之反矩陣
    # N=A'*P*A;
    N=np.dot(A.T,P)
    N=np.dot(N,A)

    #% 如果有提供未知數加權參數
    if px is not None:
        px=diag(px)
        N=N+px

    # %N_inv=inv(N);
    # U=A'*P*L;
    U=np.dot(A.T,P)
    U=np.dot(U,L)

    N_inv=inv(N)
    X=np.dot(N_inv,U)

    #% 誤差分析
    #  V=A*X-L;
    V=np.dot(A,X)-L

    #%單位權標準誤差   
    #sigma0=sqrt((V'*P*V)./(m-n));
    q=np.dot(V.T,P)
    q=np.dot(q,V)
    sigma0=math.sqrt(q/(m-n))

    #%未知數精度矩陣
    DX=N_inv*(sigma0**2)
    sX=np.diag(DX)**0.5

    #%觀測量精度分析   
    DL=np.dot(A,N_inv)
    DL=sigma0**2*np.dot(DL,A.T)
    sL=(np.diag(DL))**0.5

    return X, V, sigma0, DX, sL, P, N, U

def LSEA_Print(A,L, X, V, sigma0, DX, sL, P, N, U):
    print('***** 最小二乘法計算模組報表 *****')
    
    Matrix_print(A,'設計矩陣A')
    Matrix_print(L,'觀測矩陣L')
    
    Matrix_print(P,'加權矩陣P')
    Matrix_print(N,'法方程式 N矩陣')
    Matrix_print(U,'法方程式 U矩陣')
    Matrix_print(X,'未知數X')
    Matrix_print(V,'誤差矩陣V')

    print('\n單位權標準誤差：+/- %.4f'%(sigma0) )

    sX=np.diag(DX)**0.5
    Matrix_print2(X,sX,'未知數及精度分析',unit='m')
    Matrix_print(sL,'觀測量精度矩陣')