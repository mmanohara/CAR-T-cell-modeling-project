from pylab import *
import matplotlib
matplotlib.use('TkAgg')
'''
This program is used to find the total populations
at the end of CAR therapy for PBT025
'''
x1 = [.23, .41, .40, .42, .62, .59]#list of x1 values (for each region)
x2 = [.80, .42, .42, .56, .70, .81]#list of x2 values (for each region)
y1 = [.12, .39, .48, .42, .42, .62]#list of y1 values (for each region)
y2 = [.33, .41, .51, .57, .57, .94]#list of y2 values (for each region)
xcenter = []#define the lists needed for initial condition
xradius = []
xdistance = []
ycenter = []
yradius = []
ydistance = []
xr = zip(x1, x2)#group the x values from each region
yr = zip(y1, y2)#group the y values from each region
for x, y in xr:
    center = 0.5 * (x + y)#find the center of the region wrt x-axis
    radius = 0.5 * (y - x)#radius of the region on x-axis
    distance = y - x#width of the region on x-axis
    xcenter.append(center)#do for each region and make a list of them
    xradius.append(radius)
    xdistance.append(distance)
for x, y in yr:
    center = 0.5 * (x + y)#do the same as above but for y
    radius = 0.5 * (y - x)
    distance = y - x
    ycenter.append(center)
    yradius.append(radius)
    ydistance.append(distance)
xvalues, yvalues = meshgrid(arange(0, 1, 0.01), arange(0, 1, 0.01))#define x and y values in space
c1 = zeros([100, 100])#define the variables as 100x100 arrays
c2 = zeros([100, 100])
c3 = zeros([100, 100])
#below is the initial condition
c1 = 0.240285449163 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.290015076236 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.240285449163 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.256862254407 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.23792699968 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.272847730808 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))
c2 = 0.429355952143 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.436515491985 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.344325475769 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.41449222030 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.482494125477 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.483961689156 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))
c3 = 0.210334022555 * exp(-((xvalues-xcenter[0])**2)/(xradius[0]**2) - ((yvalues - ycenter[0])**2)/(yradius[0]**2))\
     + 0.284565645848 * exp(-((xvalues-xcenter[1])**2)/(xradius[1]**2) - ((yvalues - ycenter[1])**2)/(yradius[1]**2))\
     + 0.195839270132 * exp(-((xvalues-xcenter[2])**2)/(xradius[2]**2) - ((yvalues - ycenter[2])**2)/(yradius[2]**2))\
     + 0.206594393398 * exp(-((xvalues-xcenter[3])**2)/(xradius[3]**2) - ((yvalues - ycenter[3])**2)/(yradius[3]**2))\
     + 0.223476822196 * exp(-((xvalues-xcenter[4])**2)/(xradius[4]**2) - ((yvalues - ycenter[4])**2)/(yradius[4]**2))\
     + 0.211561909763 * exp(-((xvalues-xcenter[5])**2)/(xradius[5]**2) - ((yvalues - ycenter[5])**2)/(yradius[5]**2))
t = 25 * exp(-((xvalues - .5)**2 + (yvalues - .5)**2) / (.05**2))#first t-cell injection
c1max = 0#define max value variables
c2max = 0
c3max = 0
for x in range(100):
    for y in range(100):
        if c1[x, y] > c1max:#find max value by comparing max with each point
            c1max = c1[x, y]
        if c2[x, y] > c2max:
            c2max = c2[x, y]
        if c3[x, y] > c3max:
            c3max = c3[x, y]
c1total = 0#define total value variables
c2total = 0
c3total = 0
ttotal = 0
for x in range(100):#sum up each variable for every point
    for y in range(100):
        c1total = c1total + c1[x, y]
        c2total = c2total + c2[x, y]
        c3total = c3total + c3[x, y]
        ttotal = ttotal + t[x, y]
print(c1total)#print values of total c1, c2, c3, t and the total cancer population
print(c2total)
print(c3total)
print(ttotal)
ctotal = c1total + c2total + c3total
print(ctotal)
#plot the distributions
clf()#idk why but I need this so that the program doesn't do weird stuff
subplot(1, 3, 1)
imshow(c1, vmin = 0, vmax = c1max)#show c1
colorbar(fraction=0.046, pad=0.04)#make the colorbar the same size as the graph
title('HER2')
subplot(1, 3, 2)
imshow(c2, vmin = 0, vmax = c2max)
colorbar(fraction=0.046, pad=0.04)
title('EGFR')
subplot(1, 3, 3)
imshow(c3, vmin = 0, vmax = c3max)
colorbar(fraction=0.046, pad=0.04)
title('IL13ra2')
savefig('After-treatment.png')#save the figure as this filename
tight_layout()#make everything evenly spaced
show()