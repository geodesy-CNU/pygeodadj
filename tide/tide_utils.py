#
#從ttide預報結果取出天文潮參數(潮名、振福、相位角)
#==============================================
# ttide的輸出變數: out，可用下列指令檢查其dict結構
#
# for key in out: 
#    print(key)
from pandas.plotting import register_matplotlib_converters
import matplotlib.pyplot as plt
import numpy as np

def get_const(out,const_name):
  amp=0.0
  phase=0.0
  freq=0.0

  const_name=const_name.upper()
  if (const_name=='Z0'):
    ast='Z0'
    amp=out['z0']

  else:
    ast_name=out['nameu']
    amp_phase=out['tidecon']
    FQ=out['fu']
    n=len(ast_name)
    for i in range(n):
      ast=str(ast_name[i],'utf-8')
      if (ast.strip()==const_name):
        amp=amp_phase[i][0]
        phase=amp_phase[i][2]
        freq=FQ[i]
        break

  return [ast, amp, phase, freq]
  
def plot_obs(water_level,Station_Name):
	#water_level = df[' Water Level']['2015-01-01':]
	#water_level = df[' Water Level'][start_date:]
	ax = water_level.plot(figsize=(13, 3.5),title=Station_Name)
	ret=ax.set(ylabel='Height (m)')

def plot_analysis(out,water_level,Station_Name,datum):	
	register_matplotlib_converters()

	Z0=out['z0']
	plt.figure(figsize=[18,10])

	plt.subplot(3,1,1)
	#water_level = df[' Water Level'][start_date:]
	plt.plot(water_level, label='Observations')
	plt.legend(numpoints=1, loc='lower right')
	plt.ylim([-1,3])
	plt.title(Station_Name)

	plt.subplot(3,1,2)
	plt.title('Prediction')
	xout=out['xout']+Z0
	n=len(water_level)
	plt.plot(xout[:n])
	#plt.plot(xout[:n], label='Prediction')
	#plt.legend(numpoints=1, loc='lower right')
		
	DL=datum['DL']
	xline=[0,n]
	yline=[DL,DL]
	plt.plot(xline,yline,'r',alpha=0.5)
	
	SR=datum['SR']
	xline=[0,n]
	yline=[SR,SR]
	plt.plot(xline,yline,'g',alpha=0.5)
	
	NR=datum['NR']
	xline=[0,n]
	yline=[NR,NR]
	plt.plot(xline,yline,'b',alpha=0.5)
	plt.legend(['DL','SR','NR'],loc='upper left')
	plt.ylim([-1,3])

	plt.subplot(3,1,3)
	res=out['xres']-Z0
	plt.plot(res[:n], alpha=0.5, label='Residue')
	plt.legend(numpoints=1, loc='lower right')
	ret=plt.ylim([-0.5,0.5])
	
	return
	

def get_datum(out):
#  計算潮汐基準面
#=================================
#       MN: 平均潮差
#     MHWI: 平均高潮間隙
#       SR: 平均大潮升
#       NR: 平均小潮升
#     MLAD: 基準面起之平均海面  
#       DL: 平均基準面
#       MR: 平均潮升
#     HHWL: 平均高高潮面
#=================================
# 程式設計： Dr. Ming-Jer Huang (geodesy.cnu@gmail.com)
# 日期： 2019/10/24
#  	
	Z0=get_const(out,'Z0')
	O1=get_const(out,'O1')
	K1=get_const(out,'K1')
	M2=get_const(out,'M2')
	S2=get_const(out,'S2')
	M4=get_const(out,'M4')
	M6=get_const(out,'M6')

	DL = Z0[1] - (M2[1] + S2[1] + K1[1] + O1[1])
	SR = M2[1]*2.0 + S2[1]*2.0 + K1[1] + O1[1]
	NR = M2[1]*2.0 + K1[1] + O1[1]
	MR = M2[1]*2.0 + S2[1] + K1[1] + O1[1]

	HHWL = Z0[1] + (M2[1] + S2[1] + K1[1] + O1[1])
	MLAD = Z0[1] - DL 

	B4 = -M4[1]*np.sin((M2[2]*2.0-M4[2])*np.pi/180.0)
	B6 = -M6[1]*np.sin((M2[2]*3.0 - M6[2])*np.pi/180.0)
	MHWI = (M2[2]/28.98)+(B4*4.0/M2[1])+(B6*6.0/M2[1])
	MN = M2[1]*2.0

	dat={'Z0':Z0,'DL':DL, 'SR':SR,'NR':NR,'MR':MR,'HHWL':HHWL,'MLAD':MLAD,'MHWI':MHWI,'MN':MN}
	return dat