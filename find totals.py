# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 10:18:24 2016

@author: mohithmanohara
After getting constants, use this to get the 
initial population sizes
"""
n = 102
x1 = [.42, .42, .56, .73, .58, .78, .66, .48, .41]#values of x1
x2 = [.52, .47, .60, .76, .71, .79, .79, .57, .46]#values of x2
y1 = [.22, .30, .30, .41, .46, .61, .63, .59, .52]#values of y1
y2 = [.28, .32, .33, .47, .51, .66, .76, .68, .57]#values of y2
xcenter = []#define the stuff needed for the initial condition
xradius = []
xdistance = []
ycenter = []
yradius = []
ydistance = []
xr = zip(x1, x2)#so that x1 and x2 are grouped per region
yr = zip(y1, y2)
for x, y in xr:
    center = 0.5 * (x + y) #find the center of the region
    radius = 0.5 * (y - x) #find the radius of the region wrt x-axis
    distance = y - x #total x-width of the region
    xcenter.append(center)#do the above for each region and make a list from it
    xradius.append(radius)
    xdistance.append(distance)
for x, y in yr:
    center = 0.5 * (x + y)#same thing as above but for y
    radius = 0.5 * (y - x)
    distance = y - x
    ycenter.append(center)
    yradius.append(radius)
    ydistance.append(distance)
xvalues, yvalues = meshgrid(arange(0, 1.02, .01), arange(0, 1.02, .01))#create an array of x and y values
t = zeros([n, n])#define variables as 102x102 arrays
c0 = zeros([n, n])
c1 = zeros([n, n])
c2 = zeros([n, n])
c3 = zeros([n, n])
cco = zeros([n, n])
ctotal = zeros([n, n])
#below is the initial condition
c1 = 0.212618271575 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.259049438359 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.240741950531 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.187806106978 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.264195022836 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.216115788849 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))\
     + 0.18831285661 * exp(-((xvalues-xcenter[6])**2)/(xradius[6]**2) - ((yvalues - ycenter[6])**2)/(yradius[6]**2))\
     + 0.249550585375 * exp(-((xvalues-xcenter[7])**2)/(xradius[7]**2) - ((yvalues - ycenter[7])**2)/(yradius[7]**2))\
     + 0.225401596604 * exp(-((xvalues-xcenter[8])**2)/(xradius[8]**2) - ((yvalues - ycenter[8])**2)/(yradius[8]**2))
c2 = 0.323411575148 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.257447504943 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.631633039032 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.566447175922 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.50658600036 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.359019785716 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))\
     + 0.41393745248 * exp(-((xvalues-xcenter[6])**2)/(xradius[6]**2) - ((yvalues - ycenter[6])**2)/(yradius[6]**2))\
     + 0.534814448632 * exp(-((xvalues-xcenter[7])**2)/(xradius[7]**2) - ((yvalues - ycenter[7])**2)/(yradius[7]**2))\
     + 0.559586960997 * exp(-((xvalues-xcenter[8])**2)/(xradius[8]**2) - ((yvalues - ycenter[8])**2)/(yradius[8]**2))
c3 = 0.23148258867 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.218227417138 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.308149434751 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.284753112437 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.29713664012 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.254135214096 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))\
     + 0.245345621326 * exp(-((xvalues-xcenter[6])**2)/(xradius[6]**2) - ((yvalues - ycenter[6])**2)/(yradius[6]**2))\
     + 0.28433185191 * exp(-((xvalues-xcenter[7])**2)/(xradius[7]**2) - ((yvalues - ycenter[7])**2)/(yradius[7]**2))\
     + 0.251935908717 * exp(-((xvalues-xcenter[8])**2)/(xradius[8]**2) - ((yvalues - ycenter[8])**2)/(yradius[8]**2))

ttotal = 0#define variables of total cell type
c0total = 0
c1total = 0
c2total = 0
c3total = 0
ccototal = 0
for x in range(1, n-1):#sum up the total each cell type.
        for y in range(1, n-1):
            ttotal = ttotal + t[x, y]
            c0total = c0total + c0[x, y]
            c1total = c1total + c1[x, y]
            c2total = c2total + c2[x, y]
            c3total = c3total + c3[x, y]
            ccototal = ccototal + cco[x, y]        
print(c0total)
print(ccototal)
print(c1total)
print(c2total)
print(c3total)
print(ttotal)