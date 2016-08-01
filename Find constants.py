# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 16:21:55 2016

@author: mohithmanohara

For PBT027

This program is meant to find the constant that is
to be applied to the gaussian distribution equation. 
For each region, change the region number and change the
region letter to obtain constants for each of the three
cancer species. Then, multiply the constant by the gaussian
equation for each region for the initial condition for 
each cancer population.
"""
k = 350#arbitrary carrying capacity
Ac1 =58.4465/k#population densities, but normalized
Ac2 =88.9024/k
Ac3 =63.6321/k
Bc1 =71.2173/k
Bc2 =70.7769/k
Bc3 =59.9946/k
Cc1 =66.1774/k
Cc2 =173.6292/k
Cc3 =84.7070/k
Dc1 =51.6259/k
Dc2 =155.7103/k
Dc3 =78.2756/k
Ec1 =72.6244/k
Ec2 =139.2551/k
Ec3 =81.6797/k
Fc1 =49.3379/k
Fc2 =81.9620/k
Fc3 =58.0175/k
Gc1 =51.7652/k
Gc2 =113.7870/k
Gc3 =67.4429/k
Hc1 =68.5988/k
Hc2 =147.0148/k
Hc3 =78.1598/k
Ic1 =61.9605/k
Ic2 =153.8245/k
Ic3 =69.2545/k
x1 = [.42, .42, .56, .73, .58, .78, .66, .48, .41]#values of x1
x2 = [.52, .47, .60, .76, .71, .79, .79, .57, .46]#values of x2
y1 = [.22, .30, .30, .41, .46, .61, .63, .59, .52]#values of y1
y2 = [.28, .32, .33, .47, .51, .66, .76, .68, .57]#values of y2
A = [Ac1, Ac2, Ac3]#list of each cancer species population within the region
B = [Bc1, Bc2, Bc3]
C = [Cc1, Cc2, Cc3]
D = [Dc1, Dc2, Dc3]
E = [Ec1, Ec2, Ec3]
F = [Fc1, Fc2, Fc3]
G = [Gc1, Gc2, Gc3]
H = [Hc1, Hc2, Hc3]
I = [Ic1, Ic2, Ic3]
xvalues, yvalues = meshgrid(arange(0, 1, 0.01), arange(0, 1, 0.01))
r = 8#Region number. Reminder that region numbers start from zero.
for a in I:#change the region if you want    
    J = (1 / ((x2[r] - x1[r])*(y2[r]-y1[r]))) * 0.01 * 0.01 * sum(exp(-(((xvalues-(0.5 * (x1[r] + x2[r])))**2)/((0.5 * (x1[r] - x2[r]))**2)) - (((yvalues - (0.5 * (y2[r] + y1[r])))**2)/((0.5 * (y2[r] - y1[r]))**2))))#integrate the population density from the equation.
    K = a / J #this is the value of the constant.
    print(str(K))
