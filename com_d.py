#! /usr/bin/env python
# -*- coding: utf-8 -*-

import serial
import time
import numpy as np
def hexShow(argv):  
    result = ''  
    hLen = len(argv)  
    for i in xrange(hLen):  
        hvol = ord(argv[i])  
        hhex = '%02x'%hvol  
        result += hhex+' '  
    print ':',result  


def com_init():  
  com= serial.Serial('COM3',115200,timeout=0)
  return com  


def get_num(input):
  val=0
  if(input> 0x7fff):
    val = -((~input)&0x7fff -1)  ## 获取负数补码
  else:
    val=input
  return val

def decode_erula(data_src):
  off_set=29
  raw_data= data_src[off_set:]
  data=[]
  for x in xrange(0,7):
    data.append(ord(raw_data[x]))
  if(data[0]==0xd0):
    erual_x=(data[1]+data[2]*256)
    erual_y=(data[3]+data[4]*256)
    erual_z=(data[5]+data[6]*256)
    return  np.array([get_num(erual_x)/100.0,get_num(erual_y)/100.0,get_num(erual_z)/10.0])



def decode_acc(data_src):
  off_set = 8
  raw_data= data_src[off_set:]
  data=[]
  for x in xrange(0,7):
    data.append(ord(raw_data[x]))
  if(data[0]==0xa0):
    data_x=(data[1]+data[2]*256)
    data_y=(data[3]+data[4]*256)
    data_z=(data[5]+data[6]*256)
  else:
    data_x=0
    data_y=0
    data_z=0
  return  np.array([get_num(data_x),get_num(data_y),get_num(data_z)])

def decode_gyro(data_src):
  off_set =15
  raw_data= data_src[off_set:]
  data=[]
  for x in xrange(0,7):
    data.append(ord(raw_data[x]))
  if(data[0]==0xb0):
    data_x=(data[1]+data[2]*256)
    data_y=(data[3]+data[4]*256)
    data_z=(data[5]+data[6]*256)
  else:
    data_x=0
    data_y=0
    data_z=0
  return  np.array([get_num(data_x),get_num(data_y),get_num(data_z)])


def get_data(com):
  while(com.inWaiting()<36):
      time.sleep(0.0005)
    # time.sleep(0.001)
  data = com.readall()
  #return 
  #return hexShow(data)
  return np.append(decode_acc(data),[decode_gyro(data),decode_erula(data)])
  #  print str  
if __name__ == '__main__':  
  comp=com_init()
  import time as tm
  time_str = tm.strftime('%Y-%m-%d-%H-%M-%S',tm.localtime(tm.time()))
  file_name='arhs-'+time_str+'.txt'
 # fo = open(file_name, "wb")
  np.set_printoptions(precision=6,suppress=True)
  while 1:
      data= get_data(comp)
      stra=str(data)+'\n'
      print stra
   # get_data(comp)
    #fo.write(stra);