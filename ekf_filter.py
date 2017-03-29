# -*- coding: utf-8 -*-
from math  import sqrt,atan2,asin,sqrt
from numpy  import  pi ,array
import numpy as  np


def FastSqrtI(n=1):
	if n ==0:
		#print "norm ==0 !!"
		return 1
	else:
		return  1.0/sqrt(n)

def Sqr(n):
    return n**2
def qa_to_erula_angle(q):
	q0,q1,q2,q3=q
	RollAngle = atan2(2.0*q2*q3 - 2.0*q0*q1, 2.0*Sqr(q0) + 2.0*Sqr(q3) - 1.0);
	num = (2.0*q1*q2 - 2.0*q0*q2);
	if (num<-1):
		#print num
		num=-1.0
	elif (num>1) :
		#print num
		num=1.0
	PitchAngle = asin(num)
	YawAngle = atan2(2.0*q1*q2 - 2.0*q0*q3,  2.0*Sqr(q0) + 2.0*Sqr(q1) - 1.0);
	return np.array([PitchAngle,RollAngle,YawAngle])

def EKF_IMUGetAngle(X):
	CBn=np.zeros(9)
	rpy =np.zeros(3)
	CBn[0] = 2.0 * (X[0] * X[0] + X[1] * X[1]) - 1.0
	CBn[1] = 2.0 * (X[1] * X[2] + X[0] * X[3])
	CBn[2] = 2.0 * (X[1] * X[3] - X[0] * X[2])
	# //CBn[3] = 2.0f * (X[1] * X[2] - X[0] * X[3]);
	# //CBn[4] = 2.0f * (q0q0 + X[2] * X[2]) - 1.0f;
	#CBn[5] = 2.0 * (X[2] * X[3] + X[0] * X[1])
	CBn[5] = 2.0 * (X[2] * X[3] + X[0] * X[1])
	# //CBn[6] = 2.0f * (X[1] * X[3] + X[0] * X[2]);
	# //CBn[7] = 2.0f * (X[2] * X[3] - X[0] * X[1]);
	CBn[8] = 2.0 * (X[0] * X[0] + X[3] * X[3]) - 1.0

	#//roll

	RollAngle = atan2(CBn[5], CBn[8])

	RollAngle -=np.pi
	if(RollAngle<-np.pi):
		RollAngle+=2*np.pi
	RollAngle=-RollAngle
	#//pitch
	if (CBn[2] >= 1.0):
		PitchAngle = -np.pi/2.0
	elif (CBn[2] <= -1.0):
		PitchAngle = np.pi/2.0
	else:
		PitchAngle = asin(-CBn[2])
	#//yaw
	YawAngle = atan2(CBn[1], CBn[0])
	if (YawAngle < 0.0):
		YawAngle += np.pi*2

	if (YawAngle >= np.pi*2):
		YawAngle = 0.0
	return np.array([PitchAngle,RollAngle,YawAngle])




class ekf_filter(object):
	def __init__(self):
		self.status_dm =7
		self.measure_dm =6
		self.init_flag = 0

		self.qq=0.045 ##
		self.qw=0.025


		self.ra= 0.07 
		self.rm=0.105

		self.pq = 0.000001
		self.pw =0.000010
        ### 
        ### 对应陀螺漂移误差		
		self.Q =    np.array([
								[self.qq, 0, 0, 0, 0, 0, 0],
								[0, self.qq, 0, 0, 0, 0, 0],
								[0, 0, self.qq, 0, 0, 0, 0],
								[0, 0, 0, self.qq, 0, 0, 0],
								[0, 0, 0, 0, self.qw, 0, 0],
								[0, 0, 0, 0, 0, self.qw, 0],
								[0, 0, 0, 0, 0, 0, self.qw],
			                 ])
		### 重力和磁场误差
		self.R =  np.array([	
							[self.ra, 0, 0, 0, 0, 0],
							[0, self.ra, 0, 0, 0, 0],
							[0, 0, self.ra, 0, 0, 0],
							[0, 0, 0, self.rm, 0, 0],
							[0, 0, 0, 0, self.rm, 0],
							[0, 0, 0, 0, 0, self.rm],
			              ]) ##6,6
###
		self.I =np.eye(self.status_dm)
		self.F =np.eye(self.status_dm)
		self.H =np.zeros((self.measure_dm,self.status_dm))
		#self.H =np.zeros((self.status_dm,self.status_dm)) ##7,6
		self.P =np.array([ 
							[self.pq, 0, 0, 0, 0, 0, 0],
							[0, self.pq, 0, 0, 0, 0, 0],
							[0, 0, self.pq, 0, 0, 0, 0],
							[0, 0, 0, self.pq, 0, 0, 0],
							[0, 0, 0, 0, self.pw, 0, 0],
							[0, 0, 0, 0, 0, self.pw, 0],
							[0, 0, 0, 0, 0, 0, self.pw],
        	            ])

		self.X=np.zeros(self.status_dm)
		self.Y=np.zeros(self.measure_dm)
		self.K=np.zeros((self.status_dm,self.measure_dm))
	def filter_init(self,gyro,accel):
		nedVector = np.array([0, 0 , -1.0])
		accelVector=np.zeros(3)
		crossVector=np.zeros(3)
		#//unit accel
		norm =   FastSqrtI(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]) 

		## or use  norm = np.linalg.norm(accel)
		accelVector[0] = accel[0] * norm
		accelVector[1] = accel[1] * norm
		accelVector[2] = accel[2] * norm
     
		#cross product between accel and reference
		crossVector[0] = accelVector[1] * nedVector[2] - accelVector[2] * nedVector[1]
		crossVector[1] = accelVector[2] * nedVector[0] - accelVector[0] * nedVector[2]
		crossVector[2] = accelVector[0] * nedVector[1] - accelVector[1] * nedVector[0]

		sinwi = FastSqrtI(crossVector[0] * crossVector[0] + crossVector[1] * crossVector[1] + crossVector[2] * crossVector[2])
		crossVector[0] *= sinwi
		crossVector[1] *= sinwi
		crossVector[2] *= sinwi
		 #the angle between accel and reference is the dot product of the two vectors
		cosw = accelVector[0] * nedVector[0] + accelVector[1] * nedVector[1] + accelVector[2] * nedVector[2]
		# coshalfw = sqrt((1.0- cosw) * 0.5)
		# sinhalfw = sqrt((1.0 +cosw) * 0.5)

		coshalfw = sqrt((1.0 + cosw) * 0.5)
		sinhalfw = sqrt((1.0 -cosw) * 0.5)

		self.X[0] = coshalfw
		self.X[1] = crossVector[0] * sinhalfw
		self.X[2] = crossVector[1] * sinhalfw
		self.X[3] = crossVector[2] * sinhalfw

		norm = FastSqrtI(self.X[0] * self.X[0] + self.X[1] * self.X[1] + self.X[2] * self.X[2] + self.X[3] * self.X[3])

		self.X[0] *= norm
		self.X[1] *= norm
		self.X[2] *= norm
		self.X[3] *= norm

	def filter_update(self,gyro,accel,dt=0.002):
		halfdt =dt*0.5
		### 0.5*t(w - wbias)
		##self.X[4:6]存储的是w的
		halfdx = halfdt * (gyro[0] - self.X[4])
		halfdy = halfdt * (gyro[1] - self.X[5])
		halfdz = halfdt * (gyro[2] - self.X[6])

		neghalfdx = -halfdx
		neghalfdy = -halfdy
		neghalfdz = -halfdz

		q0 = self.X[0]
		q1 = self.X[1]
		q2 = self.X[2]
		q3 = self.X[3]
		# //////////////////////////////////////////////////////////////////////////
		# //Extended Kalman Filter: Prediction Step
		# //state time propagation
		# //Update Quaternion with the new gyroscope measurements
		#x = x +∫f(x)dt 四元素方程微分形式
		#这里相当于执行 xk=A*xk-1
		self.X[0] = q0 - halfdx * q1 - halfdy * q2 - halfdz * q3
		self.X[1] = q1 + halfdx * q0 + halfdz * q2 - halfdy * q3 
		self.X[2] = q2 + halfdy * q0 - halfdz * q1 + halfdx * q3 
		self.X[3] = q3 + halfdz * q0 + halfdy * q1 - halfdx * q2 

		halfdtq0 = halfdt * q0
		neghalfdtq0 = -halfdtq0
		halfdtq1 = halfdt * q1
		neghalfdtq1 = -halfdtq1
		halfdtq2 = halfdt * q2
		neghalfdtq2 = -halfdtq2
		halfdtq3 = halfdt * q3
		neghalfdtq3 = -halfdtq3
        #生成预测矩阵 F
        #有关系x = F*x
        #对陀螺输出值进行积分有次计算姿态角
        # q0 = q0 - wx*q1/2 -  wy*q2/2 -（wz-wbias)*q3/2 
		#self.F[0][0] =1
		self.F[0][1] = neghalfdx
		self.F[0][2] = neghalfdy
		self.F[0][3] = neghalfdz
		self.F[0][4] = halfdtq1
		self.F[0][5] = halfdtq2
		self.F[0][6] = halfdtq3

		self.F[1][0] = halfdx
		#self.F[1][1] =1.0
		self.F[1][2] = halfdz	
		self.F[1][3] = neghalfdy
		self.F[1][4] = neghalfdtq0 
		self.F[1][5] = halfdtq3
		self.F[1][6] = neghalfdtq2


		self.F[2][0] = halfdy
		self.F[2][1] = neghalfdz
	    #self.F[2][2]=1
		self.F[2][3] = halfdx
		self.F[2][4] = neghalfdtq3
		self.F[2][5] = neghalfdtq0 
		self.F[2][6] = halfdtq1

		self.F[3][0] = halfdz
		self.F[3][1] = halfdy
		self.F[3][2] = neghalfdx
	    # self.F[3][3] = 1
		self.F[3][4] = halfdtq2
		self.F[3][5] = neghalfdtq1
		self.F[3][6] = neghalfdtq0

		# //P = F*P*F' + Q
		self.P= np.dot(np.dot(self.F,self.P),self.F.T) + self.Q


		#  //////////////////////////////////////////////////////////////////////////
		# //measurement update
		# //kalman gain calculation
		# //K = P * H' / (R + H * P * H')
		# print self.H
		_2q0,_2q1,_2q2,_2q3 = 2.0 * self.X[:4]
		## 更新观察矩阵 
		### H=J(y,x)   Jacobian
	    # y = h(x)=Cbn*x
	    #计算Cbn 对X的Jacbian方程即可。H阵为6,7

	    #即计算 四元素矩阵的Jacobin 公式

		self.H[0][0] = _2q2
		self.H[0][1] = -_2q3
		self.H[0][2] = _2q0 
		self.H[0][3] = -_2q1

		self.H[1][0] = -_2q1
		self.H[1][1] = -_2q0 
		self.H[1][2] = -_2q3 
		self.H[1][3] = -_2q2

		self.H[2][0] = -_2q0 
		self.H[2][1] = _2q1 
		self.H[2][2] = _2q2 
		self.H[2][3] = -_2q3

		#  //K=P*H'/(H*P*H'+R)
		
		S=np.dot(np.dot(self.H,self.P),self.H.T) + self.R
		SI=np.linalg.inv(S)
		self.K=np.dot(self.P,np.dot(self.H.T,SI))

		# //state measurement update
		# 	//Y = Z-h(x)  //h(x) is just   
		# 	//X = X + K * Y;
		# 	//acccel on  earth  [0 0  1]  is standard g vector 
		# 	//计算重力在机体上的投影
		#   为何为负值  对应H中负值
		self.Y[0] = -2.0 * (self.X[1] * self.X[3] - self.X[0] * self.X[2])
		self.Y[1] = -2.0 * (self.X[2] * self.X[3] + self.X[0] * self.X[1])
		self.Y[2] =  1.0 -2.0 * (self.X[0] * self.X[0] + self.X[3] * self.X[3])
		#//normalize accel
		#计算实际观测到的重力坐标
		norm = FastSqrtI(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2])
		acc=np.array(accel)
		acc = acc*norm
		# accel[0] *= norm;
		# accel[1] *= norm;
		# accel[2] *= norm;
		# y-h(x) 算作误观测误差，这里本质上为用重力来校准陀螺误差
		# 对于磁力有类似的算法
		#Y阵当前为[accx,accy,accy,0,0,0]
		#待添加磁校准形成观察阵[accx,accy,accy,mx,my,mz]
		self.Y[:3]  = acc - self.Y[0:3]
		# // Update State Vector
		# X=X + K(Y-h(X))
		self.X = self.X  + np.dot(self.K,self.Y)
		
		#//normalize quaternion
		norm=FastSqrtI(self.X[0]**2 + self.X[1]**2+ self.X[2]**2 + self.X[3]**2 )
		self.X[0] *= norm
		self.X[1] *= norm
		self.X[2] *= norm
		self.X[3] *= norm
		# # 
		# 	//covariance measurement update
		# 	//P = (I - K * H) * P
		# 	//P = P - K * H * P
		# 	//or
		# 	//P=(I - K*H)*P*(I - K*H)' + K*R*K'
		self.P = np.dot(self.I - np.dot(self.K,self.H),self.P)

		q=[self.X[0],self.X[1] ,self.X[2] ,self.X[3]]
		#print q
		return  EKF_IMUGetAngle(q)
		#return qa_to_erula_angle(q)


if __name__ == '__main__': 
	filter=ekf_filter()
	gyro=np.array([0,0,0])
	acc=np.array([0,0,1])
	a=filter.filter_update(gyro,acc,0.002)
	erual= np.array([a])
	print erual*180/np.pi

