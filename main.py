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
filter1 = ekf.ekf_filter()

import madgwick_filter as mdf
filter2 =mdf.madgwick_filter(0.1) 


import time as tm
time_str = tm.strftime('%Y-%m-%d-%H-%M-%S',tm.localtime(tm.time()))
file_name='arhs-'+time_str+'.txt'

fo = open(file_name, "wb")
while 1:
	data = ser.get_data(com_port)
	acc=data[:3]
	gyro=data[3:6]
	erula_cal=data[6:]
	#print acc,gyro
	#tm.sleep (1)
	if (filter1.init_flag ==0):
		filter1.filter_init(gyro,acc)
		filter1.init_flag=1	
	else:
		data=filter1.filter_update(gyro,acc,0.001)
	   
		erual =data*180/np.pi
		stra=' '
		stra+=str(erual[0]) +' '
		stra+=str(erual[1]) +' '
		stra+=str(erual[2]) +' '



		data=filter2.filter_update(gyro,acc)
		erual =data*180/np.pi
		stra+=str(erual[0]) +' '
		stra+=str(erual[1]) +' '
		stra+=str(erual[2]) +' '
         
		stra+=str(erula_cal[0]) +' '
		stra+=str(erula_cal[1]) +' '
		stra+=str(erula_cal[2]) +' '
		stra+='\n'
		print stra
		fo.write(stra);



