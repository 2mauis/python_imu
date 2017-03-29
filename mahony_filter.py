# -*- coding: utf-8 -*-
from math  import sqrt,atan2,asin
from numpy import pi
def Sqr(x):
	return x**2
def qa_to_erula_angle(q):
	q0,q1,q2,q3=q
	RollAngle = atan2(2.0*q2*q3 - 2.0*q0*q1, 2.0*Sqr(q0) + 2.0*Sqr(q3) - 1.0);
	num = (2.0*q1*q2 - 2.0*q0*q2);
	if (num<-1):
		print num
		num=-1.0
	elif (num>1) :
		print num
		num=1.0
	PitchAngle = asin(num)
	YawAngle = atan2(2.0*q1*q2 - 2.0*q0*q3,  2.0*Sqr(q0) + 2.0*Sqr(q1) - 1.0);
	return [PitchAngle,RollAngle,YawAngle]
class mahony_filter(object):
	"""mahony filer """
	def __init__(self, kp=15.0,ki=0.005,kd=0.1):
		self.MKp=  kp
		self.MKi=  ki
		self.MKd = kd
		self.q0=1.0
		self.q1=0.0
		self.q2=0.0
		self.q3=0.0
		self.g=999
		self.exInt=0.0
		self.eyInt=0.0
		self.ezInt=0.0
		self.ex_pre =0;
		self.ey_pre =0;
		self.ez_pre =0;
		self.dT=1.0/1000 ##sample rate 1000HZ 2000°/s
		self.dTOn2 = self.dT/2.0  ##控制陀螺采样速度
		self.factor = 16.4*180.0/pi
	def filter_update(self,gyro,acc):
		if(len(gyro)<3):
			q=[self.q0,self.q1 ,self.q2 ,self.q3 ]
			return   qa_to_erula_angle(q)
		gx,gy,gz=gyro
		ax,ay,az=acc
		gx = gx/self.factor 
		gy = gy/self.factor 
		gz = gz/self.factor 
		norm = sqrt(Sqr(ax)+Sqr(ay)+Sqr(az))
		if(norm==0):
		  norm =1.0
		ax=ax/norm
		ay=ay/norm
		az=az/norm
		# estimated direction of gravity
		# vx = 2.0*(self.q1*self.q3 - self.q0*self.q2);
		# vy = 2.0*(self.q0*self.q1 + self.q2*self.q3);
		#vz = Sqr(self.q0) - Sqr(self.q1) - Sqr(self.q2) + Sqr(self.q3);

		vx = (self.q1*self.q3 - self.q0*self.q2);  
		vy = (self.q0*self.q1 + self.q2*self.q3);  
		vz = (self.q0*self.q0 - 0.5 + self.q3*self.q3);  

		''' error is sum of cross product between reference direction of field 
		and direction measured by sensor'''
		ex = (ay*vz - az*vy);
		ey = (az*vx - ax*vz);
		ez = (ax*vy - ay*vx);

		# integral error scaled integral gain
		self.exInt += ex*self.MKi#*self.dT;
		self.eyInt += ey*self.MKi#*self.dT;
		self.ezInt += ez*self.MKi#*self.dT;

		# adjusted gyroscope measurements
		gx += self.MKp*ex + self.exInt + self.MKd*(ex-self.ex_pre);
		gy += self.MKp*ey + self.eyInt + self.MKd*(ey-self.ey_pre);
		gz += self.MKp*ez + self.ezInt + self.MKd*(ez-self.ez_pre);
     
		self.ex_pre =ex 
		self.ey_pre =ey 
		self.ey_pre =ez


		# integrate quaternion rate and normalise
		self.q0 += (-self.q1*gx - self.q2*gy - self.q3*gz)*self.dTOn2;
		self.q1 += (self.q0*gx + self.q2*gz - self.q3*gy)*self.dTOn2;
		self.q2 += (self.q0*gy - self.q1*gz + self.q3*gx)*self.dTOn2;
		self.q3 += (self.q0*gz + self.q1*gy - self.q2*gx)*self.dTOn2;
		q=[self.q0,self.q1 ,self.q2 ,self.q3 ];
		#print q
		# normalise quaternion
		ReNorm =sqrt(Sqr(self.q0) + Sqr(self.q1) + Sqr(self.q2) + Sqr(self.q3));

		self.q0 = self.q0 / ReNorm;
		self.q1 = self.q1 / ReNorm;
		self.q2 = self.q2 / ReNorm;
		self.q3 = self.q3 / ReNorm;
		q=[self.q0,self.q1 ,self.q2 ,self.q3 ];
		print q
		return  qa_to_erula_angle(q)
if __name__ == '__main__': 
	filer=mahony_filter()
	print filer.filter_update([0,0,0],[0,0,1])