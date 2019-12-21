import numpy as np
from numpy.linalg import inv
import math
from .GNSSNet_utils import tri_back
from .GNSSNet_utils import compare_names

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
            
# GNSS_net клеє
class GNSS_net:
    
    # private variables
    
    # private methods
    #self.__GNSS_design_matrix()
    
    def __init__(self):
		self._ProjectName='[Empty GNSS Net Project]'
		self._num_control_pts = 0
		self._num_obs_pts = 0
		self._num_obs_baselines = 0
		self._ctrl_pts=[]
		self._obs_pts=[]
		self._pts=[]
		self._Xp=[]
		self._baselines=[]
		self._pts_name=[]
		self._A=np.zeros([1,1])
		self._L=np.zeros([1,1])
		self._P=np.zeros([1,1])
		self._sigma0=0
    
    def GNSS_read_Adat(self, path, sigma0, ComputeObs=False):
        self.__read_Adat(path, sigma0, ComputeObs)
        
    def GNSS_design_matrix(self):
        self.__GNSS_design_matrix()

    def GNSS_Report(self,X, V, sigma0, DX, sL, N, U):
        self.__GNSS_Report(X, V, sigma0, DX, sL, N, U)
        
    # Setup GNSS_design_matrix according to input data.
    def __GNSS_design_matrix(self):
        m = self._num_obs_pts*3
        n = self._num_obs_baselines*3
        k = self._num_control_pts + self._num_obs_pts
		
        #print(m,n,k)

        self._A=np.zeros([n,m])
        self._L=np.zeros([n,1])
        self._P=np.zeros([n,n])
		self._Xp=np.zeros([k,3])
		
        id=0
        for i in range(0,n,3):
            sigma=self._baselines[id][7] 
            sigma=inv(np.dot(sigma,1.0/self._sigma0**2))
            self._P[i:i+3,i:i+3]=sigma
            
            id1=self._baselines[id][0]
            id2=self._baselines[id][1]
            Fid=(id1 - self._num_control_pts)*3
            Tid=(id2 - self._num_control_pts)*3

            #print(id1,id2,Fid,Tid)
            if (id1 >= self._num_control_pts):
                self._A[i][Fid]=-1
                self._A[i+1][Fid+1]=-1
                self._A[i+2][Fid+2]=-1

            if (id2 >= self._num_control_pts):
                self._A[i][Tid]=1
                self._A[i+1][Tid+1]=1
                self._A[i+2][Tid+2]=1
            
            FromID = self._baselines[id][0]
            ToID   = self._baselines[id][1]
            
            X1=self._pts[FromID][2]
            Y1=self._pts[FromID][3]
            Z1=self._pts[FromID][4]
            X2=self._pts[ToID][2]
            Y2=self._pts[ToID][3]
            Z2=self._pts[ToID][4]
            
            DX=X2-X1
            DY=Y2-Y1
            DZ=Z2-Z1
            self._L[i,0]  = self._baselines[id][4]-DX        
            self._L[i+1,0]= self._baselines[id][5]-DY
            self._L[i+2,0]= self._baselines[id][6]-DZ
            id=id+1                    
            

    # initial a GNSS Network from Aadj file.
    def __read_Adat(self, path,sigma0,ComputeObs):
        self._sigma0=sigma0
        self._pts_name=[]
        f=open(path, "r")
        self._ProjectName = f.readline()
        items=f.readline().split(' ')
        
        # Load Project infomation: ProjectName
        self._num_control_pts = int(items[0])
        self._num_obs_pts = int(items[1])
        self._num_obs_baselines = int(items[2])
        
        self._ctrl_pts=[]
        # Load control points.
        # Point ID start from 0
        for i in range(self._num_control_pts):
            items=f.readline().split(' ')
            ID   = i
            name = items[0]
            X    = float(items[1])
            Y    = float(items[2])
            Z    = float(items[3])
            ctrl = [ID, name, X, Y, Z]
            self._ctrl_pts.append(ctrl)
            self._pts.append(ctrl)
            self._pts_name.append(name)
        
        # load baselines
        self._baselines=[]
        for i in range(self._num_obs_baselines):
            items   = f.readline().split(' ')
            FromID  = 0
            ToID    = 0
            fromName= items[0]
            toName  = items[1]
            dX      = float(items[2])
            dY      = float(items[3])
            dZ      = float(items[4])
            sigma   = np.zeros(6)
            for j in range(5,11):
                sigma[j-5]=float(items[j])
            
            # put sigma data as upper triangle matrix
            sigma_tri = tri_back(sigma)
            sigma_full= sigma_tri+sigma_tri.T
            for i in range(3):
                sigma_full[i,i]=sigma_full[i,i]/2.0
                
            baseline=[FromID,ToID,fromName,toName,dX,dY,dZ,sigma_full]
            self._baselines.append(baseline)
        
        #print(self._baselines)
        # load initial observation points.
        # Observation ID start after number of control points.
        self._obs_pts=[]
        if (ComputeObs==True):
            for i in range(self._num_obs_pts):
                items=f.readline().split(' ')
                ID   = self._num_control_pts + i
                name = items[0]
                X    = float(items[1])
                Y    = float(items[2])
                Z    = float(items[3])
                obs  = [ID, name, X, Y, Z]
                self._obs_pts.append(obs)
                self._pts.append(obs)
                self._pts_name.append(name)
        else:
            for i in range(self._num_obs_pts):
                items=f.readline().split(' ')
                ID   = self._num_control_pts + i
                name = items[0]
                X    = float(items[1])
                Y    = float(items[2])
                Z    = float(items[3])
                obs  = [ID, name, X, Y, Z]
                self._obs_pts.append(obs)
                self._pts.append(obs)
                self._pts_name.append(name)
        #print(self._obs_pts)
        
        f.close()
        # 把base_lines裡面的測站名稱& ID放進obs_pts內
        pts_name=self._pts_name

        for i in range(self._num_obs_baselines):
            FromName = self._baselines[i][2]
            ToName   = self._baselines[i][3]
            self._baselines[i][0]=compare_names(FromName, pts_name)
            self._baselines[i][1]=compare_names(ToName, pts_name)
        #print(self._baselines)
        
        if (ComputeObs==True):
            # 用已知值估算未知測站概略座標
            for i in range(self._num_obs_baselines):
                id1=self._baselines[i][0]
                id2=self._baselines[i][1]
                #print(i,id1,id2,self._pts[id1][2],self._pts[id2][2])
                if (id1 >= self._num_control_pts):
                    if (self._pts[id1][2]==0):
                        self._pts[id1][1]=self._baselines[i][2]
                        self._pts[id1][2]=self._pts[id2][2]-self._baselines[i][4]
                        self._pts[id1][3]=self._pts[id2][3]-self._baselines[i][5]
                        self._pts[id1][4]=self._pts[id2][4]-self._baselines[i][6]
                        #print(self._pts[id1])
                        
                if (id2 >= self._num_control_pts):    
                    if (self._pts[id2][2]==0):
                        self._pts[id2][1]=self._baselines[i][3]
                        self._pts[id2][2]=self._pts[id1][2]+self._baselines[i][4]
                        self._pts[id2][3]=self._pts[id1][3]+self._baselines[i][5]
                        self._pts[id2][4]=self._pts[id1][4]+self._baselines[i][6]
                        #print(self._pts[id2])
            
            for i in range(self._num_obs_pts):
                k=i+self._num_control_pts
                self._obs_pts[i][2]=self._pts[k][2]
                self._obs_pts[i][3]=self._pts[k][3]
                self._obs_pts[i][4]=self._pts[k][4]
            #print(self._obs_pts)
         
    def __GNSS_Report(self,X, V, sigma0, DX, sL, N, U):
        print('***** GNSS 最小self._baselines二乘法 網形平差報表 *****')

        print('\nGNSS Net專案名稱: %s'%(self._ProjectName))

        print('\n控制點數: %d'%(self._num_control_pts))
        print('\n控制點座標')
		
        for i in range(self._num_control_pts):
            ctrl=self._ctrl_pts[i]
            print('%d %s\t%14.5f\t%14.5f\t%14.5f'%(ctrl[0],ctrl[1],ctrl[2],ctrl[3],ctrl[4]))

        print('\n觀測點數: %d'%(self._num_obs_pts))
        print('\n觀測點概略座標')
        for i in range(self._num_obs_pts):
            obs=self._obs_pts[i]
            print('%d %s\t%14.5f\t%14.5f\t%14.5f'%(obs[0],obs[1],obs[2],obs[3],obs[4]))

        print('\n觀測基線數: %d'%(self._num_obs_baselines))
        print('\n基線觀測資料')
        for i in range(self._num_obs_baselines):
            BSL=self._baselines[i]
            print('%s(%d)  ->  %s(%d)\t%12.4f\t%12.4f\t%12.4f'%
                  (BSL[2],BSL[0],BSL[3],BSL[1],BSL[4],BSL[5],BSL[6]))

        Matrix_print(self._A,'設計矩陣A')
        Matrix_print(self._L,'觀測矩陣L')

        Matrix_print(self._P,'加權矩陣P')
        Matrix_print(N,'法方程式 N矩陣')

        print('\n單位權標準誤差：+/- %.5f'%(sigma0) )

        sX=np.diag(DX)**0.5
        print('\n控制點及觀測點座標')
        for i in range(self._num_control_pts):
            ctrl=self._ctrl_pts[i]
            print('%d %s %14.5f %14.5f %14.5f'%(ctrl[0],ctrl[1],ctrl[2],ctrl[3],ctrl[4]))
			self._Xp[i,0]=ctrl[2]
			self._Xp[i,1]=ctrl[3]
			self._Xp[i,2]=ctrl[4]

        id=0    
        for i in range(0,self._num_obs_pts*3,3):
            obs=self._obs_pts[id]
            print('%d %s %14.5f %14.5f %14.5f %6.5f %6.5f %6.5f'
                  %(obs[0],obs[1],obs[2]+X[i],obs[3]+X[i+1],obs[4]+X[i+2],
                    sX[i],sX[i+1],sX[i+2]))
			self._Xp[id+self._num_control_pts,0]=obs[2]+X[i]
			self._Xp[id+self._num_control_pts,1]=obs[3]+X[i+1]
			self._Xp[id+self._num_control_pts,2]=obs[4]+X[i+2]					
            id+=1
        print('\n基線平差結果與精度分析')
        id=0    
        for i in range(0,self._num_obs_baselines*3,3):
            BSL=self._baselines[id]
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
        
    # ProjectName
    @property
    def ProjectName(self):
        return self._ProjectName
    @ProjectName.setter
    def ProjectName(self, ProjectName):
        self._ProjectName = ProjectName

    # num_control_pts
    @property
    def num_control_pts(self):
        return self._num_control_pts
    @num_control_pts.setter
    def num_control_pts(self,num_control_pts):
        self._num_control_pts = num_control_pts

    # num_obs_pts
    @property
    def num_obs_pts(self):
        return self._num_obs_pts
    @num_obs_pts.setter
    def num_obs_pts(self,num_obs_pts):
        self._num_obs_pts = num_obs_pts

    # num_obs_baselines
    @property
    def num_obs_baselines(self):
        return self._num_obs_baselines
    @num_obs_baselines.setter
    def num_obs_baselines(self,num_obs_baselines):
        self._num_obs_baselines = num_obs_baselines
        
    @property
    def ctrl_pts(self):
        return self._ctrl_pts
    
    @property
    def baselines(self):
        return self._baselines
    
    @property
    def obs_pts(self):
        return self._obs_pts
		
	@property
    def Xp(self):
        return self._Xp
 
    @property
    def sigma0(self):
        return self._sigma0
    
    @property
    def A(self):
        return self._A
    
    @property
    def L(self):
        return self._L
    
    @property
    def P(self):
        return self._P
    
