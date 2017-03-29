
# -*- coding: utf-8 -*-
from math  import sqrt,atan2,asin
from numpy  import  pi,array

def Sqr(x):
	return x**2
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
	return array([PitchAngle,RollAngle,YawAngle])
###def madgwick_filter
class madgwick_filter(object):
	"""docstring for ClassName"""
	def __init__(self, gyroMeasError = 0.1):  ##0.1 for mpu 9250
		self.sampleFreq = 100 ## sample rate type 10Hz
		self.dotime =1.0/self.sampleFreq
		
		self.q0  =1.0
		self.q1  =0
		self.q2  =0
		self.q3  =0
		self.factor = 16.4*180.0/pi
		#Gyroscope measurement error (in degrees per second). 
		#//Compute beta.  
		self.beta=  gyroMeasError#sqrt(3.0 / 4.0) * (pi * (gyroMeasError / 180.0));  
	def filter_update(self,gyro,acc):
		if (len(gyro)!=3) :
			q = [self.q0,self.q1 ,self.q2 ,self.q3 ]
			return   qa_to_erula_angle(q)
		gx,gy,gz=gyro
		gx = gx/self.factor
		gy = gy/self.factor
		gz = gz/self.factor
		ax,ay,az=acc
		# recipNorm;
		# s0, s1, s2, s3
		# qDot1, qDot2, qDot3, qDot4
		# _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3
        
		if ((gx==0.0) and (gy==0.0) and (gz==0.0)) or (len(gyro)==0):
			q=[self.q0,self.q1 ,self.q2 ,self.q3 ]
			return   qa_to_erula_angle(q)
	   # Rate of change of quaternion from gyroscope
		qDot1 = 0.5 * (-self.q1 * gx - self.q2 * gy - self.q3 * gz)
		qDot2 = 0.5 * (self.q0 * gx + self.q2 * gz - self.q3 * gy)
		qDot3 = 0.5 * (self.q0 * gy - self.q1 * gz + self.q3 * gx)
		qDot4 = 0.5 * (self.q0 * gz + self.q1 * gy - self.q2 * gx)

	# Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
		if(((ax != 0.0) and (ay != 0.0) and (az != 0.0))):

			# Normalise accelerometer measurement
			recipNorm = 1.0/sqrt(ax * ax + ay * ay + az * az)
			ax *= recipNorm
			ay *= recipNorm
			az *= recipNorm  

			# Auxiliary variables to avoid repeated arithmetic
			_2q0 = 2.0 * self.q0
			_2q1 = 2.0 * self.q1
			_2q2 = 2.0 * self.q2
			_2q3 = 2.0 * self.q3
			_4q0 = 4.0 * self.q0
			_4q1 = 4.0 * self.q1
			_4q2 = 4.0 * self.q2
			_8q1 = 8.0 * self.q1
			_8q2 = 8.0 * self.q2
			q0q0 = self.q0 * self.q0
			q1q1 = self.q1 * self.q1
			q2q2 = self.q2 * self.q2
			q3q3 = self.q3 * self.q3

			# Gradient decent algorithm corrective step
			s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
			s1 = _4q1 * q3q3 - _2q3 * ax + 4.0 * q0q0 * self.q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
			s2 = 4.0 * q0q0 * self.q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az
			s3 = 4.0 * q1q1 * self.q3 - _2q1 * ax + 4.0 * q2q2 * self.q3 - _2q2 * ay
			recipNorm = 1.0/sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) # normalise step magnitude
			s0 *= recipNorm
			s1 *= recipNorm
			s2 *= recipNorm
			s3 *= recipNorm

			# Apply feedback step
			qDot1 -= self.beta * s0
			qDot2 -= self.beta * s1
			qDot3 -= self.beta * s2
			qDot4 -= self.beta * s3

		# Integrate rate of change of quaternion to yield quaternion
		self.q0 += qDot1 * self.dotime
		self.q1 += qDot2 * self.dotime
		self.q2 += qDot3 * self.dotime
		self.q3 += qDot4 * self.dotime
		#Normalise quaternion
		recipNorm = 1.0/sqrt(self.q0 * self.q0 + self.q1 * self.q1 + self.q2 * self.q2 + self.q3 * self.q3)
		self.q0 *= recipNorm
		self.q1 *= recipNorm
		self.q2 *= recipNorm
		self.q3 *= recipNorm
		q=[self.q0,self.q1 ,self.q2 ,self.q3 ]
		
		return  qa_to_erula_angle(q)
if __name__ == '__main__': 
	filter=madgwick_filter()
	print [filter.q0,filter.q1 ,filter.q2 ,filter.q3 ]
	filter.filter_update([0,0,0],[0,0,1])



			