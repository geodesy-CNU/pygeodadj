import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

pi=np.pi

# 度分秒轉換函數
def dms(radian):
    ss=radian*206264.8063
    deg=math.floor(ss/3600.0)
    tt=ss-deg*3600.0
    mm=math.floor(tt/60.0)
    sec=tt-mm*60.0
    return deg, mm, sec

# 度分秒轉度度
def dms2d(deg,min,sec):
    dd=deg+min/60.0+sec/3600.0
    return dd

# 方位角距離計算函數
def azdis(xa,ya,xb,yb):
    dx=float(xb-xa)
    dy=float(yb-ya)
    dis=math.sqrt(float(dx**2+dy**2))  

    # 東西南北四個象限   
    if ((abs(dx) < 0.000001) & (dy > 0.0) ):
        alpha=0.0
        dd,mm,ss=dms(alpha)
        return dd, mm, ss, dis
    elif ((abs(dx) < 0.000001) & (dy < 0.0) ):
        alpha=pi
        dd,mm,ss=dms(alpha)
        return dd, mm, ss, dis  
    elif ((abs(dy) < 0.000001) & (dx > 0.0)):
        alpha=pi/2
        dd,mm,ss=dms(alpha)
        return dd, mm, ss, dis
    elif ((abs(dy) < 0.000001) & (dx < 0.0)):
        alpha=3*pi/2
        dd,mm,ss=dms(alpha)
        return dd, mm, ss, dis

    #一般情形
    ang=math.atan(abs(dx)/abs(dy));    
    if((dx > 0.0) & (dy > 0.0) ):
        alpha=ang  
    elif ((dx > 0.0) & (dy < 0.0) ):
        alpha=pi-ang  
    elif ((dx < 0.0) & (dy < 0.0) ):
        alpha=pi+ang  
    elif ((dx < 0.0) & (dy > 0.0) ):
        alpha=2*pi-ang

    dd,mm,ss=dms(alpha)
    return dd, mm, ss, dis

# 極座標法定位計算核心
def azdis_compute(x0,y0,az0,dd,dist):
	az=az0+dd
	if (az>360.0):
		z=az-360.0

	az=az*pi/180.0
	dx=dist*math.sin(az)
	dy=dist*math.cos(az)
	xp=x0+dx
	yp=y0+dy
	
	return xp,yp


# 讀取觀測檔案
def azdis_pos_df(path):
	# 讀取控制點
	#data_dir='drive/My Drive/Colab Notebooks/Courses/adjustment/'
	#df=pd.read_csv(data_dir+'azdis_pos_cnu.txt',header=0)
	df=pd.read_csv(path,header=0)

	# 觀測點 (CNU DG:02) TWD97
	xa=df.loc[0,'xp']
	ya=df.loc[0,'yp']

	# 參考原方向 (CNU: GD52)
	xb=df.loc[1,'xp']
	yb=df.loc[1,'yp']

	# 呼叫計算方位角距離程式
	deg,mm,sec,dis=azdis(xa,ya,xb,yb)
	az0=dms2d(deg,mm,sec)
	df.loc[1,'deg']=deg
	df.loc[1,'min']=mm
	df.loc[1,'sec']=sec
	df.loc[1,'dist']=dis

	# 點位計算
	#print(df.shape[0])
	for index in range(2,df.shape[0]):
		pt_name=df.loc[index,'PtName']
		dd=df.loc[index,'deg']
		mm=df.loc[index,'min']
		ss=df.loc[index,'sec']
		dist=df.loc[index,'dist']
		dd=dms2d(dd,mm,ss)

		xp,yp=azdis_compute(xa,ya,az0,dd,dist)
		#print(xp,yp)
		df.loc[index,'xp']=xp
		df.loc[index,'yp']=yp

	return df

# 方位角距離計算繪圖
def plot_azdis(xa,ya,xb,yb,dd,mm,ss,dis):
    # plot
    #xa=float(xa)
    #ya=float(ya)
    #xb=float(xb)
    #yb=float(yb)

    dy=float((ya+yb)/2.0);
    dx=float((xb-xa)/50.0);

    X=np.array([xa,xb])
    Y=np.array([ya,yb])
    Y2=np.array([ya,ya+dy])
    plt.figure(figsize=(10,10))
    plt.plot(X,Y,marker='o',markersize=12,color='blue')
    plt.plot(X,Y,'r',linewidth=2.0);

    plt.title('Azimuth & distance ')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')

    plt.text(xa+dx,ya,'A',size=16)
    plt.text(xb,yb+dx,'B',size=16)
    txt_str='Length= %.2f (m)'%(dis)
    plt.text((xa+xb)/2+dx,(ya+yb)/2,txt_str,va="baseline", size=16)
    txt_str='Azimuth=%d$^{o}$%d\'%.1f\"'%(dd,mm,ss)
    plt.text(xa+dx*2,ya+dy/5,txt_str,va="baseline", size=16)
    plt.arrow(xa,ya,0,dy,overhang=3)
    plt.show()
    
# 方位角距離計算繪圖
def plot_azdis_pos(xa,ya,xb,yb,dd,mm,ss,dis):
    # plot
    #xa=float(xa)
    #ya=float(ya)
    #xb=float(xb)
    #yb=float(yb)

    dy=float((ya+yb)/2.0);
    dx=float((xb-xa)/50.0);

    X=np.array([xa,xb])
    Y=np.array([ya,yb])
    Y2=np.array([ya,ya+dy])
    plt.figure(figsize=(10,10))
    plt.plot(X,Y,marker='o',markersize=12,color='blue')
    plt.plot(X,Y,'r',linewidth=2.0);

    plt.title('Azimuth & distance ')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')

    plt.text(xa+dx,ya,'A',size=16)
    plt.text(xb,yb+dx,'B',size=16)
    txt_str='Length= %.2f (m)'%(dis)
    plt.text((xa+xb)/2+dx,(ya+yb)/2,txt_str,va="baseline", size=16)
    txt_str='Azimuth=%d$^{o}$%d\'%.1f\"'%(dd,mm,ss)
    plt.text(xa+dx*2,ya+dy/5,txt_str,va="baseline", size=16)
    plt.arrow(xa,ya,0,dy,overhang=3)
    plt.show()