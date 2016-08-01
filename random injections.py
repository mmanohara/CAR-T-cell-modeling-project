# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:22:18 2016

@author: mohithmanohara

This program is meant to simulate a model
of the interaction between CAR T-cells and 
3 different cancer species. HER2 will be
represented in the program as c1, EGFR 
as c2, and Il13ra2 as c3. Images are
to be saved and compiled into a movie file.

The t-cells are placed randomly in the space in this program.
"""
import matplotlib
matplotlib.use('TkAgg')
from pylab import *
mu1 = 0.005#rate of movement of T-cells
mu2 = 0
mu3 = 0.005
muco = 0.005
d = 0.01#death rate of t cells
n = 102#size of space
k = 350 #Carrying capacity of the cancer cells
K = 1
a = 0.35#rate of growth of each cancer species
b1 = 0.6#death rate of HER2
b2 = 0#death rate of EGFR
b3 = 0.6#death rate of IL13a2
bd = 0.6#death rate of cancer showing coexpression.
Dt = 0.001#diffusion rates of the T-cells and the cancer species
Dc = 0.002
dh = .01 #spacial resolution
dt = 0.01#temporal resolution
q =  25#T-cell distribution coefficient
xcenter = random()
ycenter = random()
xvalues, yvalues = meshgrid(arange(0, 1.02, .01), arange(0, 1.02, .01))
def initialize():
    global t, c0, c1, c2, c3, cco, ctotal, nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal, ttotal, c0total, c1total, c2total, c3total, ccototal, tplot, c0plot, c1plot, c2plot, c3plot, ccoplot, step, check, tmax, c0max, c1max, c2max, c3max, ccomax, tpend, c0pend, c1pend, c2pend, c3pend, ccopend, steppend
    t = zeros([n, n])#set zero values for each variable in the space
    c0 = zeros([n, n])
    c1 = zeros([n, n])
    c2 = zeros([n, n])
    c3 = zeros([n, n])
    cco = zeros([n, n])
    ctotal = zeros([n, n])#below we set up the initial condition.
    c1 = .396606 * exp(-((xvalues-0.405)**2)/(.135**2) - ((yvalues - .33)**2)/(.24**2))\
         + .373651 * exp(-((xvalues-0.575)**2)/(.035**2) - ((yvalues - .39)**2)/(.05**2))\
         + .403440 * exp(-((xvalues-0.65)**2)/(.04**2) - ((yvalues - .385)**2)/(.085**2))\
         + .739448 * exp(-((xvalues-0.57)**2)/(.02**2) - ((yvalues - .53)**2)/(.02**2))\
         + .553523 * exp(-((xvalues-0.6)**2)/(.01**2) - ((yvalues - .515)**2)/(.015**2))\
         + .500782 * exp(-((xvalues-0.31)**2)/(.03**2) - ((yvalues - .595)**2)/(.025**2))\
         + .478544 * exp(-((xvalues-0.45)**2)/(.03**2) - ((yvalues - .625)**2)/(.035**2))\
         + .490583 * exp(-((xvalues-0.375)**2)/(.015**2) - ((yvalues - .695)**2)/(.015**2))\
         + .443507 * exp(-((xvalues-0.465)**2)/(.025**2) - ((yvalues - .675)**2)/(.005**2))\
         + .432523 * exp(-((xvalues-0.6)**2)/(.18**2) - ((yvalues - .755)**2)/(.115**2))\
         + .399770 * exp(-((xvalues-0.445)**2)/(.075**2) - ((yvalues - .915)**2)/(.045**2))
    c2 = .606893 * exp(-((xvalues-0.405)**2)/(.135**2) - ((yvalues - .33)**2)/(.24**2))\
         + .429825 * exp(-((xvalues-0.575)**2)/(.035**2) - ((yvalues - .39)**2)/(.05**2))\
         + .618103 * exp(-((xvalues-0.65)**2)/(.04**2) - ((yvalues - .385)**2)/(.085**2))\
         + .625116 * exp(-((xvalues-0.57)**2)/(.02**2) - ((yvalues - .53)**2)/(.02**2))\
         + .502399 * exp(-((xvalues-0.6)**2)/(.01**2) - ((yvalues - .515)**2)/(.015**2))\
         + .316215 * exp(-((xvalues-0.31)**2)/(.03**2) - ((yvalues - .595)**2)/(.025**2))\
         + .743633 * exp(-((xvalues-0.45)**2)/(.03**2) - ((yvalues - .625)**2)/(.035**2))\
         + .782661 * exp(-((xvalues-0.375)**2)/(.015**2) - ((yvalues - .695)**2)/(.015**2))\
         + .708838 * exp(-((xvalues-0.465)**2)/(.025**2) - ((yvalues - .675)**2)/(.005**2))\
         + .329171 * exp(-((xvalues-0.6)**2)/(.18**2) - ((yvalues - .755)**2)/(.115**2))\
         + .567646 * exp(-((xvalues-0.445)**2)/(.075**2) - ((yvalues - .915)**2)/(.045**2))
    c3 = .331261 * exp(-((xvalues-0.405)**2)/(.135**2) - ((yvalues - .33)**2)/(.24**2))\
         + .297917 * exp(-((xvalues-0.575)**2)/(.035**2) - ((yvalues - .39)**2)/(.05**2))\
         + .303581 * exp(-((xvalues-0.65)**2)/(.04**2) - ((yvalues - .385)**2)/(.085**2))\
         + .413251 * exp(-((xvalues-0.57)**2)/(.02**2) - ((yvalues - .53)**2)/(.02**2))\
         + .333327 * exp(-((xvalues-0.6)**2)/(.01**2) - ((yvalues - .515)**2)/(.015**2))\
         + .249389 * exp(-((xvalues-0.31)**2)/(.03**2) - ((yvalues - .595)**2)/(.025**2))\
         + .343188 * exp(-((xvalues-0.45)**2)/(.03**2) - ((yvalues - .625)**2)/(.035**2))\
         + .243203 * exp(-((xvalues-0.375)**2)/(.015**2) - ((yvalues - .695)**2)/(.015**2))\
         + .401408 * exp(-((xvalues-0.465)**2)/(.025**2) - ((yvalues - .675)**2)/(.005**2))\
         + .263288 * exp(-((xvalues-0.6)**2)/(.18**2) - ((yvalues - .755)**2)/(.115**2))\
         + .244822 * exp(-((xvalues-0.445)**2)/(.075**2) - ((yvalues - .915)**2)/(.045**2))
    for y in range(n):
        for x in range(n):
            ctotal[x, y] = c0[x, y] + c1[x, y] + c2[x, y] + c3[x, y] + cco[x, y]
    t = q * exp(-((xvalues - xcenter)**2 + (yvalues - ycenter)**2) / (.05**2))
    nextt = zeros([n, n])#must define the "next" variables for the space.
    nextc0 = zeros([n, n])
    nextc1 = zeros([n, n])
    nextc2 = zeros([n, n])
    nextc3 = zeros([n, n])
    nextcco = zeros([n, n])
    nextctotal = zeros([n, n])
    tplot = zeros([n-2, n-2])#Only plotting 100x100, so this is what we need to plot
    c0plot = zeros([n-2, n-2])
    c1plot = zeros([n-2, n-2])
    c2plot = zeros([n-2, n-2])
    c3plot = zeros([n-2, n-2])
    ccoplot = zeros([n-2, n-2])
    ttotal = 0#total populations
    c0total = 0
    c1total = 0
    c2total = 0
    c3total = 0
    ccototal = 0
    check = 50
    step = 0
    for x in range(1, n-1):#sum up the total each cell type.
        for y in range(1, n-1):
            ttotal = ttotal + t[x, y]
            c0total = c0total + c0[x, y]
            c1total = c1total + c1[x, y]
            c2total = c2total + c2[x, y]
            c3total = c3total + c3[x, y]
            ccototal = ccototal + cco[x, y]
    for x in range(1, n-1):#values that are actually plotted
        for y in range(1, n-1):
            tplot[x-1, y-1] = t[x, y]
            c0plot[x-1, y-1] = c0[x, y]
            c1plot[x-1, y-1] = c1[x, y]
            c2plot[x-1, y-1] = c2[x, y]
            c3plot[x-1, y-1] = c3[x, y]
            ccoplot[x-1, y-1] = cco[x, y]
    tmax = 0
    c0max = 0
    c1max = 0
    c2max = 0
    c3max = 0
    ccomax = 0
    for x in range(100):#find the max value to use for the colorbar
        for y in range(100):
            if tplot[x, y] > tmax:
                tmax = tplot[x, y]
            if c0plot[x, y] > c0max:
                c0max = c0plot[x, y]
            if c1plot[x, y] > c1max:
                c1max = c1plot[x, y]
            if c2plot[x, y] > c2max:
                c2max = c2plot[x, y]
            if c3plot[x, y] > c3max:
                c3max = c3plot[x, y]
            if ccoplot[x, y] > ccomax:
                ccomax = ccoplot[x, y]
    tpend = [ttotal]
    c0pend = [c0total]
    c1pend = [c1total]
    c2pend = [c2total]
    c3pend = [c3total]
    ccopend = [ccototal]
    steppend = [step]
def observe():
    global t, c0, c1, c2, c3, cco, ctotal, nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal, ttotal, c0total, c1total, c2total, c3total, ccototal, tplot, c0plot, c1plot, c2plot, c3plot, ccoplot, step, check, tmax, c0max, c1max, c2max, c3max, ccomax, tpend, c0pend, c1pend, c2pend, c3pend, ccopend, steppend
    clf()
    subplot(5, 5, 1)
    imshow(tplot, vmin=0, vmax = tmax)#show t-cells. 
    title('CAR T-cells')
    colorbar()
    subplot(5, 5, 3)
    imshow(c1plot, vmin = 0, vmax = c1max)#show HER2
    title('HER2')
    colorbar()
    subplot(5, 5, 5)
    imshow(c2plot, vmin = 0, vmax = c2max)#show EGFR
    title('EGFR')
    colorbar()
    subplot(5, 5, 11)
    imshow(c3plot, vmin = 0, vmax = c3max)#show IL13ra2
    title('IL13ra2')
    colorbar()
    subplot(5, 5, 21)
    imshow(c0plot, vmin = 0, vmax = c0max)#show no expression
    title('No expression')
    colorbar()
    subplot(5, 5, 13)
    imshow(ccoplot, vmin = 0, vmax = ccomax)#show coexpression
    title('coexpression')
    colorbar()
    subplot(5, 5, 23)
    plot(steppend, tpend, 'b.', label = 'CAR T-cells')    
    plot(steppend, c0pend, 'r.', label ='No expression')
    plot(steppend, c1pend, 'g.', label ='HER2')
    plot(steppend, c2pend, 'y.', label ='EGFR')
    plot(steppend, c3pend, 'm.', label ='IL13ra2')
    plot(steppend, ccopend, 'k.', label = 'coexpression')
    legend(bbox_to_anchor=(1.05, 1), prop = {'size':8}, loc=2, borderaxespad=0.)
    title('total population sizes')
    xlabel('Time')
    ylabel('population size')    
def update():
    global t, c0, c1, c2, c3, cco, ctotal, nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal, ttotal, c0total, c1total, c2total, c3total, ccototal, tplot, c0plot, c1plot, c2plot, c3plot, ccoplot, step, check, tmax, c0max, c1max, c2max, c3max, ccomax, tpend, c0pend, c1pend, c2pend, c3pend, ccopend, steppend  
    xcenter = random()#set random x and y for where the new dose is placed
    ycenter = random()
    if step<7.001 and step > 6.999 or step <14.001 and step > 13.999 or step < 21.001 and step > 20.999:
        t = t + q * exp(-((xvalues - xcenter)**2 + (yvalues - ycenter)**2) / (.05**2))#add a new dose every 7 days
    for x in range(1, n-1):
        for y in range(1, n-1):
            tC, tR, tL, tU, tD = t[x, y], t[(x + 1), y], t[(x - 1), y], \
                                 t[x, (y + 1)], t[x, (y - 1)]#setting values for neighbor cells
            c0C, c0R, c0L, c0U, c0D = c0[x, y], c0[(x + 1), y], c0[(x - 1), y], \
                                      c0[x, (y + 1)], c0[x, (y - 1)]            
            c1C, c1R, c1L, c1U, c1D = c1[x, y], c1[(x + 1), y], c1[(x - 1), y], \
                                      c1[x, (y + 1)], c1[x, (y - 1)]
            c2C, c2R, c2L, c2U, c2D = c2[x, y], c2[(x + 1), y], c2[(x - 1), y], \
                                      c2[x, (y + 1)], c2[x, (y - 1)]      
            c3C, c3R, c3L, c3U, c3D = c3[x, y], c3[(x + 1), y], c3[(x - 1), y], \
                                      c3[x, (y + 1)], c3[x, (y - 1)]
            ccoC, ccoR, ccoL, ccoU, ccoD = cco[x, y], cco[(x + 1), y], cco[(x - 1), y], \
                                           cco[x, (y + 1)], cco[x, (y - 1)]                          
            ctotalC = ctotal[x, y]                                                   
            tLapNum = tU + tD + tL + tR - 4 * tC#values of Laplacians
            c0LapNum = c0U + c0D + c0R + c0L - 4 * c0C
            c1LapNum = c1U + c1D + c1R + c1L - 4 * c1C
            c2LapNum = c2U + c2D + c2R + c2L - 4 * c2C
            c3LapNum = c3U + c3D + c3R + c3L - 4 * c3C
            ccoLapNum = ccoU + ccoD + ccoR + ccoL - 4 * ccoC
            #Be warned; the next few equations are HORRIFYINGLY LONG. Have a good day :) ‘’’the equations have the same comments as the samelocbiCAR, so check that out. No I am not being lazy, what are you talking about?’’’
            nextt[x, y] = tC + (-d * tC - mu1 * (((tR - tL) / (2 * dh))*((c1R - c1L) / (2 * dh)) + ((tU - tD) / (2 * dh)) * ((c1U - c1D) / (2 * dh)) + tC * (c1LapNum / dh**2)) - mu2 * (((tR - tL) / (2 * dh))*((c2R - c2L) / (2 * dh)) + ((tU - tD) / (2 * dh)) * ((c2U - c2D) / (2 * dh)) + tC * (c2LapNum / dh**2)) - mu3 * (((tR - tL) / (2 * dh))*((c3R - c3L) / (2 * dh)) + ((tU - tD) / (2 * dh)) * ((c3U - c3D) / (2 * dh)) + tC * (c3LapNum / dh**2))- muco * (((tR - tL) / (2 * dh))*((ccoR - ccoL) / (2 * dh)) + ((tU - tD) / (2 * dh)) * ((ccoU - ccoD) / (2 * dh)) + tC * (ccoLapNum / dh**2)) + Dt * (tLapNum/dh**2)) * dt
            nextc0[x, y] = c0C + (a * c0C * (1 - (ctotalC)/K) + Dc * (c0LapNum / dh**2)) * dt
            nextc1[x, y] = c1C + (a * c1C * (1 - (ctotalC)/K) - b1 * c1C * tC + Dc * (c1LapNum / dh**2)) * dt
            nextc2[x, y] = c2C + (a * c2C * (1 - (ctotalC)/K) - b2 * c2C * tC + Dc * (c2LapNum / dh**2)) * dt
            nextc3[x, y] = c3C + (a * c3C * (1 - (ctotalC)/K) - b3 * c3C * tC + Dc * (c3LapNum / dh**2)) * dt
            nextcco[x, y] = ccoC + (a * ccoC * (1 - (ctotalC)/K) - bd * ccoC * tC + Dc * (ccoLapNum/dh**2)) * dt
            nextctotal[x, y] = c0C + c1C + c2C + c3C + ccoC
    ttotal = 0
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
    for x in range(1, n-1):
        for y in range(1, n-1):#update plotted variables
            tplot[x-1, y-1] = t[x, y]
            c0plot[x-1, y-1] = c0[x, y]
            c1plot[x-1, y-1] = c1[x, y]
            c2plot[x-1, y-1] = c2[x, y]
            c3plot[x-1, y-1] = c3[x, y]
            ccoplot[x-1, y-1] = cco[x, y]
    tmax = 0
    c0max = 0
    c1max = 0
    c2max = 0
    c3max = 0
    ccomax = 0
    for x in range(100):
        for y in range(100):
            if tplot[x, y] > tmax:
                tmax = tplot[x, y]
            if c0plot[x, y] > c0max:
                c0max = c0plot[x, y]
            if c1plot[x, y] > c1max:
                c1max = c1plot[x, y]
            if c2plot[x, y] > c2max:
                c2max = c2plot[x, y]
            if c3plot[x, y] > c3max:
                c3max = c3plot[x, y]
            if ccoplot[x, y] > ccomax:
                ccomax = ccoplot[x, y]
    t, c0,c1, c2, c3, cco, ctotal, nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal = nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal, t, c0, c1, c2, c3, cco, ctotal
    step = step + dt
    check = check + 1#increase image counter thing. 
    if check == 50:
        tpend.append(ttotal)
        c0pend.append(c0total)
        c1pend.append(c1total)
        c2pend.append(c2total)
        c3pend.append(c3total)
        ccopend.append(ccototal)
        steppend.append(step)
def takepic():#defining when to take the pic, and what the pic is saved as.
    global t, c0, c1, c2, c3, cco, ctotal, nextt, nextc0, nextc1, nextc2, nextc3, nextcco, nextctotal, ttotal, c0total, c1total, c2total, c3total, ccototal, tplot, c0plot, c1plot, c2plot, c3plot, ccoplot, step, check
    if check == 50:
        observe()
        savefig(str(step) + 'random injections.png')
        check = 0

def picmode():#the mode when we want to save a bunch of figures
    initialize()
    takepic()
    while step < 28:
        update()
        takepic()
    figtext(0.8, 0.4, 'T =' + str(ttotal))#when done, display the total population size on the figure.
    figtext(0.8, 0.35, 'None =' + str(c0total))
    figtext(0.8, 0.3, 'HER2 =' + str(c1total))
    figtext(0.8, 0.25, 'EGFR =' + str(c2total))
    figtext(0.8, 0.2, 'IL13ra2 =' + str(c3total))
    figtext(0.8, 0.15, 'coexpression =' + str(ccototal))
    savefig('finalrandominjections.png')
picmode()