# -*- 2oding: utf-8 -*-
import com_d as ser

import time as tm
import numpy  as np
com_port =ser.com_init()
'''
import mahony_filter as mhf
filter =mhf.mahony_filter(10,0.01,1)

import madgwick_filter as mdf
filter =mdf.madgwick_filter(0.4) 
'''
import ekf_filter as ekf
filter =ekf.ekf_filter()


while 1:
	data = ser.get_data(com_port)
	acc=data[:3]
	gyro=data[3:6]
	erula_cal=data[6:]
	#print acc,gyro
	#tm.sleep (1)
	if (filter.init_flag ==0):
		filter.filter_init(gyro,acc)
		filter.init_flag=1	
	else:
		a=filter.filter_update(gyro,acc,0.002)
		erual= np.array([a])
		print erual*180/np.pi,erula_cal
