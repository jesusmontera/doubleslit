# shows how to calculate and plot double slit:
    # 1. intensities( slit2Intensity )
    # 2. fringes ( dslit2fringes )
    # 3. waves and their intersections( drawwaves )    
    # from wave length, slit separation, slit width and screen distance
    # so u can see that the 3 methods plot in the same place ok, the common
    # NOTE : common aproximation is not good enougth when the screen is near,
    # but becuase I use slit width and angles it works fine.
    # to get information on the concepts visit the page and view the video:
    # https://en.universaldenker.org/formulas/1119
    
     
    
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['toolbar'] = 'None'
#from numba import jit
from curvesintersection import curveintersec # to get intersections points on 2 curves
from math import asin,sin,degrees,tan,atan,cos

HPLANK=6.626E-34
MASSELECTRON=9.1e-31
MASSPROTON=1.67262192e-27

def broglieWL(m, v):    
    return HPLANK/(m*v)

def dslit2fringes(g,a,wavelength,bConstructive=True):
    # calculate fringes distances from center
    # same as dslit2fringesaprox but more precisse
     #g slit separtaion
     # a screen dfistance
    fringesdist=[]
    for m in range(1,6):         
        if bConstructive:
           sinalfa=(wavelength*m)/g           
        else:
           sinalfa=(wavelength*(m-0.5))/g
            
        alfa =  asin(sinalfa)        
        x = tan(alfa)*a            
        fringesdist.append(x)
    return fringesdist   

##def dslit2fringesaprox(g,a,wavelength,bConstructive=True):
##    # calculate fringes distances from center
##    # this aproximation is valid only if the screen dist(a) >> slit separation (g)
##     #g slit separtaion
##     # a screen dfistance
##    fringesdist=[]
##    for m in range(1,6):         
##        if bConstructive:
##           ds=m * wavelength
##        else:
##           ds = (m-0.5)* wavelength
##        x = (a *ds) /g# fringe separation
##        
##        print(x)
##        fringesdist.append(x)
##    return fringesdist   
        
def calcPointIntensity(puntoY,wl,L,d,w):
    # L screen distance
    # d slit separation
    # w slit width
    IMAXATCENTER=1.
    angulo=atan(puntoY/L)
    a=cos( (np.pi * d * sin(angulo)) / wl)
    bn = sin(np.pi * w * sin(angulo)/wl)
    bd= np.pi * w * sin(angulo)/wl
    return IMAXATCENTER * a * a * pow(bn/bd,2.0)
def slit2Intensity(wl,L,d,w,NPOINTS,screenmaxY):
    # intensity from center to right with slit width care
    # L screen distance
    # d slit separation
    # w slit width
    # wl wave length
    intensity=np.zeros(NPOINTS)
    y = np.linspace(0, screenmaxY, NPOINTS)    
    intensity[0]=1.
    for i in range(1,NPOINTS):
        intensity[i] = calcPointIntensity(y[i],wl,L,d,w)      
    return intensity

##def dslit1fringes(w,a,wavelength):
##     #w = slit width 
##     # a  = screen dfistance
##     # x = distance between screen center max and minimun
##     for m in range(1,6):
##         x  = m * wavelength * a /w
##         print("distance from center to minimun",m, " = ",x)
         


def intersecPoints2(w1X,w1Y,w2X,w2Y):
    ip=[]
    n1 = len(w1X)
    n2 = len(w2X)
    for i in range(n1):        
        for l in range(n2):
            x,y=curveintersec( w1X[i], w1Y[i],w2X[l],w2Y[l])
            for k in range(len(x)):
                ip.append((x[k],y[k]))            
    return ip
            
    
            
        
def drawwaves(wl, n,slitseparation,slitwidth,screendist,screenmaxY):    

    # calculate waves curves at wave maximum(up) 
    
    NPOINTSCURVE=360
    a = np.linspace(0, np.pi, NPOINTSCURVE) - np.pi/2
    
    wXmax = np.zeros([n,NPOINTSCURVE]) # slit 1 wave coordinates at max
    wYmax = np.zeros([n,NPOINTSCURVE])
    wXmin = np.zeros([n,NPOINTSCURVE]) # slit 1 wave coordinates at min
    wYmin = np.zeros([n,NPOINTSCURVE])
    
    rMax=wl/2 # the first wave after slit is at maximun(up)
    rMin=wl # the first wave minimun after slit is 1 wl
    for i in range(n):        
        wXmax[i] =  rMax * np.cos(a)
        wYmax[i] =  rMax * np.sin(a) +slitseparation/2 

        wXmin[i] =  rMin * np.cos(a)
        wYmin[i] =  rMin * np.sin(a) +slitseparation/2         
        rMax+=wl
        rMin+=wl

    
    
    # plot slits    
    
    ymin = -screenmaxY
    ymax = screenmaxY
    
    ay=slitseparation/2+slitwidth/2.
    by=-slitseparation/2- slitwidth/2 
    plt.plot([0,0], [ymax,ay], color="black",linewidth=4)
    plt.plot([0,0], [ymin,by], color="black",linewidth=4)
    
    plt.plot([0,0], [ay-slitwidth,by+slitwidth],linewidth=4, color="black")
    
    
    # plot waves
    for i in range(n):        
        #plot waves at maximun
        plt.plot(wXmax[i],wYmax[i], color="red")  # wave slit 1 at max     
        plt.plot(wXmax[i],wYmax[i]-slitseparation, color="green") # wave slit 2 at max     
        #plot waves at minimun
        plt.plot(wXmin[i],wYmin[i], color="red",linestyle='dashed',alpha=0.6)  # wave slit 1 at max     
        plt.plot(wXmin[i],wYmin[i]-slitseparation, color="green",linestyle='dashed',alpha=0.6) # wave slit 2 at max
        
        
    # plot points at waves instersection 
    ipmax= intersecPoints2(wXmax,wYmax,wXmax,wYmax-slitseparation )
    for i in range(len(ipmax)):
        if ipmax[i][0] <= screendist:
            plt.scatter(ipmax[i][0],ipmax[i][1], color="blue",s=20)

    ipmin= intersecPoints2(wXmin,wYmin,wXmin,wYmin-slitseparation )
    for i in range(len(ipmin)):
        if ipmin[i][0] <= screendist:
            plt.scatter(ipmin[i][0],ipmin[i][1], color="black",s=20)

                      
# here starts the program, the screen distance is very near(at 9 wave lengths)
# so I just have to plot 14 waves to see the interfencere
wl=4.8e-7
screendist=wl*9
slitseparation=wl*6.1
slitwidth=slitseparation/10.
numwaves=14
screenmaxY = wl*numwaves/1.5
SCREEN_Y_RES=1000



fig = plt.figure(figsize = (7, 7))
fig.tight_layout()

plt.gca().set_aspect('equal')
plt.xlim((-wl, wl*numwaves + wl*2 ))
plt.ylim(( -screenmaxY, screenmaxY))

Yintensity=slit2Intensity(wl,screendist,slitseparation,slitwidth,SCREEN_Y_RES,screenmaxY)
Yintensity = np.concatenate((Yintensity[::-1], Yintensity)) * 1.5e-6 + screendist



drawwaves(wl,numwaves,slitseparation,slitwidth,screendist,screenmaxY) # plot waves

# plot intensities at screen
Xintensity=np.linspace(-screenmaxY,screenmaxY,SCREEN_Y_RES*2)
plt.plot(Yintensity,Xintensity,color="black",linewidth=2)



# plot fringes in screen
fringesdist= dslit2fringes(g = slitseparation,a = screendist, wavelength=wl,bConstructive=True)
plt.plot([screendist,screendist], [-wl/2,wl/2], color="black",linewidth=4)
plt.plot([0,screendist], [0,0], color="black",linewidth=1)
for i in range(4):
    ymin= fringesdist[i]-wl/2
    ymax= fringesdist[i]+wl/2
    plt.plot([screendist,screendist], [ymin,ymax], color="black",linewidth=3)
    plt.plot([screendist,screendist], [-ymin,-ymax], color="black",linewidth=3)    
plt.show()        
        
