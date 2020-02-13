import numpy as np 
from numpy import linalg as LA
import math
             
bodies = [[           0      ,          0    ,              0       ,         0           ,    1.5707963      ,    0            , 0, 0 ,0  ,  0    ,   0     ,       0         ,     0       ,     0 ,0, 0],
          [-0.00000000000000  , 0.03907757912148 ,    0.44992952528334  ,-0.00000000000003 ,  1.59273072796596 ,  0.00000000000000 ,  0, 0, 0  ,  1400 ,  2235  ,      2269 ,          379  ,         0, 0 ,0],
          [ 0.51527369815677 ,  1.14433409974204  , 0.22508473490622 ,  0.00033816594169 ,  1.59272979230103 ,  0.00553507938926 ,  0, 0, 0  ,    3.4405 ,3.48E-004 ,  3.7933E-002 ,   3.7856E-002 ,  0, 0, 0],
          [ 0.58030888627625 ,  1.13500527938571 , 0.64496567637505,  -0.00053196546317 ,  1.59273303719436,   0.00868390769008 ,  0, 0, 0  ,     2.187,  4.071E-003 , 1.5218E-002 ,   1.152E-002  ,  0, 0 ,0],
          [-0.51527369815681 ,  1.14433409974207,   0.22508473490617 , -0.00033816594173 ,  1.59272979230102,  -0.00553507938972 ,  0 ,0 ,0,      3.4405, 3.48E-004,   3.7933E-002,    3.7856E-002,   0, 0, 0],
          [-0.58030888627619 ,  1.13500527938579,  0.64496567637501 ,  0.00053196546308 ,  1.59273303719436 , -0.00868390768979,   0, 0, 0  ,     2.187 , 4.071E-003,  1.5218E-002  ,  1.152E-002 ,   0 ,0 ,0],
          [ 0.69018601168402 ,  1.14008763903236 ,  0.41599454115878,  -0.00197637945134 ,   1.59322298488907 ,  0.00124933871770,   0, 0 , 0 ,     11.678, 0.120541 ,   7.7691E-002  ,  6.5647E-002 ,  0 ,0 ,0],
          [-0.69018601168390 ,  1.14008763903241,  0.41599454115874 ,  0.00197637944848 ,  1.59322298488898 , -0.00124933871848 ,  0 ,0 ,0 ,      11.678, 0.120541  ,  7.7691E-002 ,   6.5647E-002 ,  0 ,0 ,0],
          [ 0.00000000000003 ,  1.00173276295906 , 0.34284774902796 , -0.00000000000001 ,  1.57079632679490 ,  0.00000000000028 ,  0 ,0 ,0   ,    0.575  ,2.3E-005  ,  3.876E-003   ,  3.876E-003  ,  0 ,0, 0],
          [ 0.54034604459266 , -1.05866849342128 ,  0.33512395377186 , -0.00000000000003 ,  1.59273072796596 ,  0.12123884615529 ,  0, 0, 0 ,      3.4405, 3.48E-004 ,  3.7933E-002 ,   3.7856E-002 , 0 ,0 ,0],
          [-0.54034604459270 , -1.05866849342124 ,  0.33512395377193 , -0.00000000000003 , 1.59273072796596 , -0.12123884615548 , 0 ,0, 0 ,        3.4405 ,3.48E-004 ,  3.7933E-002 ,   3.7856E-002 , 0, 0, 0] ]   

ground = bodies[0]
    

# revolute joint matrix
# [body1 body2 locinbdy1 locinbdy2 jtframe_angles1 jtframe_angles2]

RevJts= [[1, 2 ,  0.34  , -0.2, -1.1,  -0.175 ,   0.05  ,    0 ,     0 ,   -0.061086523  , 0  , 0 , -0.061086523, 0],
         [1, 3 ,  0.4708, 0.16, -1.1, -0.1096 , -0.01   ,   0  ,    0  ,   0.061086523   ,0   , 0 ,  0.061086523, 0],
         [1, 4 , -0.34  , -0.2,-1.1 ,  0.175  ,  0.05   ,   0  ,    0  ,  -0.061086523   , 0  , 0 , -0.061086523, 0], 
         [1, 5 , -0.4708, 0.16, -1.1,   0.1096,  -0.01  ,    0  ,    0 ,    0.061086523  , 0  , 0 ,  0.061086523, 0],
         [1, 9,  0.34  , -0.1,  1.1,  -0.2   ,   0.015 ,    0  ,    0 ,    0            , 0  , 0 ,  0          , 0],
         [1, 10 ,-0.34  , -0.1,  1.1,  0.2    ,  0.015  ,   0  ,    0  ,   0             , 0  , 0 ,  0          , 0]]

# spherical joint matrix
# [body1 body2 locinbdy1 locinbdy2]

SphJts =[[ 2, 6,  0.175 , -0.05, -0.022,     0    ,  -0.24 ,  -0.022],
	 [ 3, 6,  0.1096,  0.01,  0.04 ,     0    ,  0.24  ,  0.04 ],
	 [ 4, 7, -0.175 , -0.05, -0.022,     0    ,  -0.24 ,  -0.022],
	 [ 5, 7, -0.1096,  0.01,  0.04 ,     0    ,   0.24 ,   0.04 ],
	 [ 6, 8, -0.02  , -0.07,  0.14 ,     0.67 ,   0    ,   0    ],
	 [ 7, 8,  0.02  , -0.07,  0.14 ,    -0.67 ,   0    ,   0    ]]

	
# Revolute joint parameters
JtGeomPar1=0.1 #used in r=JtGeomPar(1)*lrgdim
JtGeomPar2=0.3 #used in re=JtGeomPar(2)*r
JtGeomPar3=1   #used in l=JtGeomPar(3)*r
JtGeomPar4=0.3 #used in le=JtGeomPar(4)*l
JtGeomPar5=0.3 #used in lh=JtGeomPar(5)*l
JtGeomPar6=0.1 #used in gap=JtGeomPar(6)*r
  
# Spherical joint parameters
JtGeomPar7=60 #used in ang=JtGeomPar(7)
  
JtGeomPar8=24 #used in ncptspts=JtGeomPar(8)
JtGeomPar9=24 #used in ninoutsgg=JtGeomPar(9)
JtGeomPar10=12#used in noutsg=JtGeomPar(10)
JtGeomPar11=10#used in ngrd=JtGeomPar(11)


RJScl=0.1*np.ones((6,1))
SJScl=0.1*np.ones((6,1))
BScl=0.5*np.ones((11,1))
BClrs=2*np.ones((11,1))

nb=len(bodies)
nrj=len(RevJts)
noutsgj=len(SphJts)
  
  # Creating bodies
lrgdim=0
I=np.zeros((3,3))
for i in range(1,nb):
    m=bodies[i][9]
    I[0][0]=bodies[i][10]
    I[1][1]=bodies[i][11]
    I[2][2]=bodies[i][12]
    I[0][1]=bodies[i][13]
    I[1][0]=bodies[i][13]
    I[1][2]=bodies[i][14]
    I[2][1]=bodies[i][14]
    I[2][0]=bodies[i][15]
    I[0][2]=bodies[i][15]
    if (m <= 0):
      bx=[ 0 , 0 ]
      bx=np.reshape(bx,(3,1))
      by=[ 0 , 0 ]
      by=np.reshape(by,(3,1))
      bz=[ 0 , 0 ]
      bz=np.reshape(bz,(3,1))
    else:
      [eigval,eigvec]=LA.eig(I)              #check syntax
      e1=eigval[0] 
      e2=eigval[1]
      e3=eigval[2]
      v1=eigvec[0]
      v2=eigvec[1]
      v3=eigvec[2]
      v1=v1/(LA.norm(v1))
      v2=v2/(LA.norm(v2))
      v3=v3/(LA.norm(v3))
  
      s1=0.5*np.math.sqrt(6*(e2+e3-e1)/m)
      s2=0.5*np.math.sqrt(6*(e3+e1-e2)/m)
      s3=0.5*np.math.sqrt(6*(e1+e2-e3)/m)
  
      bx=[[s1,s1,s1,s1,s1,s1,s1,s1,-s1,-s1,-s1,-s1],[-s1,-s1,-s1,-s1,s1,s1,s1,s1,-s1,-s1,-s1,-s1]]
      by=[[s2,-s2,-s2,s2,s2,-s2,-s2,s2,s2,-s2,-s2,s2],[s2,-s2,-s2,s2,-s2,-s2,s2,s2,-s2,-s2,s2,s2]]
      bz=[[s3,s3,-s3,-s3,s3,s3,-s3,-s3,s3,s3,-s3,-s3],[s3,s3,-s3,-s3,s3,-s3,-s3,s3,s3,-s3,-s3,s3]]
      bx=BScl[i]*bx
      by=BScl[i]*by
      bz=BScl[i]*bz
      R = np.zeros(shape=(3,3))
      R=[v1,v2,v3]
      n=len(bx[1])  # syntax
      r=[[bx[0], bx[1]], [by[0], by[1]], [bz[0], bz[1]]]
      r=np.reshape(r,(3,24))
      xyz=np.matmul(R,r)
      bx=np.reshape(xyz[0],(2,12))
      by=np.reshape(xyz[1],(2,12))
      bz=np.reshape(xyz[2],(2,12))
      if (BScl[i] > 0):
          a=np.max([bx ,by ,bz]/BScl[i])
          
          lrgdim=a
  
      r=JtGeomPar1*lrgdim
      re=JtGeomPar2*r
      l=JtGeomPar3*r
      le=JtGeomPar4*l
      lh=JtGeomPar5*l
      gap=JtGeomPar6
      ang=JtGeomPar7
  
      ncptspts=JtGeomPar8
      ninoutsgg=JtGeomPar9
      noutsg=JtGeomPar10
      ngrd=JtGeomPar11
      if (nrj > 0):
        for irj in range(1,nrj):
          RJS=RJScl[irj]
          b1=RevJts[irj][0]
          pos11  =[0.34  ,- 0.2 , - 1.1 ]
          pos11=np.reshape(pos11,(3,1)) 
          phi1=RevJts[irj][8]
          th1=RevJts[irj][9]
          drcn1=[math.sin(th1)*math.cos(phi1),math.sin(th1)*math.sin(phi1),math.cos(th1)]
          if (LA.norm(drcn) == 0):
            print("Revelute jt direction vector is null.  Aborting.")
    
          else:
            czaxis=drcn/(LA.norm(drcn))
          if (LA.norm(pos1) != 0):
            cxdir=-pos1
          else:
            cxdir[0][0]=random.uniform(0,1)-0.5
            cxdir[1][0]=random.uniform(0,1)-0.5
            cxdir[2][0]=random.uniform(0,1)-0.5
  
         cxdir=cxdir/(LA.norm(cxdir))
  
         cyaxis[0][0]=czaxis[1]*cxdir[2]-czaxis[2]*cxdir[1]
         cyaxis[1][0]=czaxis[2]*cxdir[0]-czaxis[0]*cxdir[2]
         cyaxis[2][0]=czaxis[0]*cxdir[1]-czaxis[1]*cxdir[0]

         if (LA.norm(cyaxis) == 0):
           yover=0
          while (yover == 0):
            if (czaxis[0] != 0):
              cyaxis[1][0]=1
              cyaxis[2][0]=1
              cyaxis[0]=-(czaxis[1]*cyaxis[1]+czaxis[2]*cyaxis[2])/czaxis[0]
              yover=1
            elif (czaxis[1] != 0): 
              cyaxis[0][0]=1
              cyaxis[2][0]=1
              cyaxis[1]=-(czaxis[0]*cyaxis[0]+czaxis[2]*cyaxis[2])/czaxis[1]
              yover=1
            elif (czaxis[2] !=0): 
              cyaxis[0][0]=1
              cyaxis[1][0]=1
              cyaxis[2]=-(czaxis[0]*cyaxis[0]+czaxis[1]*cyaxis[1])/czaxis[2]
              yover=1
  
        cyaxis=cyaxis/norm(cyaxis)

        cxaxis[0][0]=cyaxis[1]*czaxis[2]-cyaxis[2]*czaxis[1]
        cxaxis[1][0]=cyaxis[2]*czaxis[0]-cyaxis[0]*czaxis[2]
        cxaxis[2][0]=cyaxis[0]*czaxis[1]-cyaxis[1]*czaxis[0]
  
        R=[cxaxis ,cyaxis ,czaxis]

  #Axial segments
       twopi=2*(22/7)
       thetas=[0:twopi-twopi/2/(noutsg-1):twopi/(noutsg-1)]
       lby2=(l*RJS)/2
       g=r*RJS*gap
       z1=lby2-g
       z2=-lby2+g
       cx=[] 
       cy=[]
       cz=[]
       rr=r+g
  
       for i in range(1,noutsg-1):
         x1=rr*np.cos(thetas(i))
         x2=x1
         y1=rr*np.sin(thetas(i))
         y2=y1
         cx=[cx ,[[x1],[x2]]]  #tranoutsgpos1e 
         cy=[cy ,[[y1],[y2]]]
         cz=[cz ,[[z1],[z2]]]
  

  #Circles
       phis=[0:twopi:twopi/(ncpts-1)]  #list creation   #check syntax
       np=len(phis,2)
       crcx=rr*np.cos(phis)
       crcx1=crcx(1:np-1)
       crcx2=crcx(2:np)
       crcy=rr*np.sin(phis)
       crcy1=crcy(1:np-1)
       crcy2=crcy(2:np)
       crcz1=z1*np.ones(crcx1)
       crcz2=-crcz1
  
       cx=[cx ,[[crcx1],[crcx2]], [[crcx1],[crcx2]]]     
       cx=[cy ,[[crcy1],[crcy2]], [[crcy1],[crcy2]]] 
       cz=[cz ,[[crcz1],[crcz1]] ,[[crcz2],[crcz2]]]
  
  #Handle
       hx=[[rr],[rr+lh]]
       hy=[[0],[0]]
       hz=[[0],[0]]
       cx=[cx , hx]
       cy=[cy , hy]
       cz=[cz , hz]
       cxn1=R(1,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos1[0]
       cxn2=R(1,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos1[0]
       cyn1=R(2,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos1[1]
       cyn2=R(2,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos1[1]
       czn1=R(3,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos1[2]
       czn2=R(3,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos1[2]
       n=len(cxn1[1])
       cx=[[cxn1, 0],[cxn2 , cxn2(n)]]
       cy=[[cyn1, 0],[cyn2 , cyn2(n)]]
       cz=[[czn1, 0],[czn2 , czn2(n)]]
 
  

  





      