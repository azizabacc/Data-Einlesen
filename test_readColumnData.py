# script test_readColumnData.py

''' test data input from text file with module PhyPraKit.readColumnData
..  author:: Guenter Quast <g.quast@kit.edu>
'''

from PhyPraKit import readColumnData
import numpy as np, matplotlib.pyplot as plt

fname='xyData.dat'
ncols=2

data_array, info_dict =\
  readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
x=data_array[:,0] # 1st column
ex=data_array[:,1] # 2nd column
y=data_array[:,2] # 3rd column
ey=data_array[:,3] # 4th column

print "Title=", info_dict['*TITLE']
print "x= ", x
print "y= ", y
print "ex= ", ex
print "ey= ", ey
    
print "keywords found:"
for key in info_dict:
  if (info_dict[key]!=None): print key,':', info_dict[key]

# make a plot
fig=plt.figure(1, figsize=(5.,5.))
plt.errorbar(x , y, xerr=ex, yerr=ey, fmt='b.')
# show systematic errors
exs = x * np.float32(info_dict['*xRelCor'])
eys = y * np.float32(info_dict['*yRelCor'])
plt.errorbar(x , y, xerr=exs, yerr=eys, fmt='r.')
# and finally, the axis labels
plt.xlabel(info_dict['*xLabel'] +' / ' + info_dict['*xUnit'], size='x-large')
plt.ylabel(info_dict['*yLabel'], size='x-large')

plt.show()

