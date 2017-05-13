#! /usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.preprocessing import normalize

# x = np.array([
# 	            [1,2],
# 	            [3,4]
# 	         ])
# print x.T*x





# norm1 = x / np.linalg.norm(x)
# print np.linalg.norm(x)
# print norm1
# norm2 = normalize(x[:,np.newaxis], axis=0).ravel()
# print np.all(norm1 == norm2)

import numpy as np
# b=np.ones((7,6))
# c =np.ones((6,7))
# a = np.dot(b,c)
# print a.shape
stra ="-0.094727,0.008301,1.003906,0.011706,-0.006385,0.007449"
a=np.fromstring(stra, dtype=float, sep=', ')
print a
print str(a)+'a'
