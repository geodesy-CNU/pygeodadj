
import numpy as np
from numpy.linalg import inv
import math

# Compare sting
def compare_names(name,obs_pts):
    id=0
    for i in range(len(obs_pts)):
        if (name==obs_pts[i]):
            id=i
            break
    
    return id

# 由輸入的sigma，組成上三角矩陣
def tri_back(sigma):
    n=len(sigma)
    m=n
    for i in range(n):
        m=m-i
        if (m==0):
            break

    tri=np.zeros([i,i])
    id=0
    id2=0
    for k in range(i,0,-1):
        #print(tri[id,id:id+k])
        tri[id,id:id+k]=sigma[id2:id2+k]
        id=id+1
        id2=id2+k

    return tri
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
            
def GNSS_LSEA(A,L,P):
    m,n=A.shape
    #print(m,n)
    
    #% 未知數的答解 
    #% 法方程式之反矩陣
    # N=A'*P*A;
    N=np.dot(A.T,P)
    N=np.dot(N,A)

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

def GNSS_Report(wolf_gnss_net, A,L, X, V, sigma0, DX, sL, P, N, U):
    print('***** GNSS 最小二乘法 網形平差報表 *****')

    print('\nGNSS Net專案名稱: %s'%(wolf_gnss_net.ProjectName))
    
    print('\n控制點數: %d'%(wolf_gnss_net.num_control_pts))
    print('\n控制點座標')
	
	Xp=np.zeros([wolf_gnss_net.num_control_pts+wolf_gnss_net.num_obs_pts,3])
    for i in range(wolf_gnss_net.num_control_pts):
        ctrl=wolf_gnss_net.ctrl_pts[i]
        print('%d %s\t%14.5f\t%14.5f\t%14.5f'%(ctrl[0],ctrl[1],ctrl[2],ctrl[3],ctrl[4]))

    print('\n觀測點數: %d'%(wolf_gnss_net.num_obs_pts))
    print('\n觀測點概略座標')
    for i in range(wolf_gnss_net.num_obs_pts):
        obs=wolf_gnss_net.obs_pts[i]
        print('%d %s\t%14.5f\t%14.5f\t%14.5f'%(obs[0],obs[1],obs[2],obs[3],obs[4]))

    print('\n觀測基線數: %d'%(wolf_gnss_net.num_obs_baselines))
    print('\n基線觀測資料')
    for i in range(wolf_gnss_net.num_obs_baselines):
        BSL=wolf_gnss_net.baselines[i]
        print('%s(%d)  ->  %s(%d)\t%12.4f\t%12.4f\t%12.4f'%
              (BSL[2],BSL[0],BSL[3],BSL[1],BSL[4],BSL[5],BSL[6]))
        
    Matrix_print(A,'設計矩陣A')
    Matrix_print(L,'觀測矩陣L')

    Matrix_print(P,'加權矩陣P')
    Matrix_print(N,'法方程式 N矩陣')

    print('\n單位權標準誤差：+/- %.5f'%(sigma0) )

    sX=np.diag(DX)**0.5
    print('\n控制點及觀測點座標')
    for i in range(wolf_gnss_net.num_control_pts):
        ctrl=wolf_gnss_net.ctrl_pts[i]
        print('%d %s %14.5f %14.5f %14.5f'%(ctrl[0],ctrl[1],ctrl[2],ctrl[3],ctrl[4]))
		Xp[i,0]=ctrl[2]
		Xp[i,1]=ctrl[3]
		Xp[i,2]=ctrl[4]
    
    id=0    
    for i in range(0,wolf_gnss_net.num_obs_pts*3,3):
        obs=wolf_gnss_net.obs_pts[id]
        print('%d %s %14.5f %14.5f %14.5f %6.5f %6.5f %6.5f'
              %(obs[0],obs[1],obs[2]+X[i],obs[3]+X[i+1],obs[4]+X[i+2],
                sX[i],sX[i+1],sX[i+2]))
		Xp[id+wolf_gnss_net.num_control_pts,0]=obs[2]+X[i]
		Xp[id+wolf_gnss_net.num_control_pts,1]=obs[3]+X[i+1]
		Xp[id+wolf_gnss_net.num_control_pts,2]=obs[4]+X[i+2]
        id+=1
    print('\n基線平差結果與精度分析')
    id=0    
    for i in range(0,wolf_gnss_net.num_obs_baselines*3,3):
        BSL=wolf_gnss_net.baselines[id]
        dx=BSL[4]+V[i]
        dy=BSL[5]+V[i+1]
        dz=BSL[6]+V[i+2]
        print('%s(%d) -> %s (%d) %12.5f %12.5f %12.5f %6.5f %6.5f %6.5f'
              %(BSL[2],BSL[0],BSL[3],BSL[1],dx,dy,dz,sL[i],sL[i+1],sL[i+2]))
        id+=1
    #Matrix_print2(X,sX,'未知數及精度分析',unit='m')
    #Matrix_print(sL,'觀測量精度矩陣')
    Matrix_print(U,'法方程式 U矩陣')
    Matrix_print(X,'未知數X')
    Matrix_print(V,'誤差矩陣V')
	
	