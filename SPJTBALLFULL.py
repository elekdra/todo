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

RevJts= [[2, 3 ,  0.34  , -0.2, -1.1,  -0.175 ,   0.05  ,    0 ,     0 ,   -0.061086523  , 0  , 0 , -0.061086523, 0],
         [2, 4 ,  0.4708, 0.16, -1.1, -0.1096 , -0.01   ,   0  ,    0  ,   0.061086523   ,0   , 0 ,  0.061086523, 0],
         [2, 5 , -0.34  , -0.2,-1.1 ,  0.175  ,  0.05   ,   0  ,    0  ,  -0.061086523   , 0  , 0 , -0.061086523, 0], 
         [2, 6 , -0.4708, 0.16, -1.1,   0.1096,  -0.01  ,    0  ,    0 ,    0.061086523  , 0  , 0 ,  0.061086523, 0],
         [2, 10,  0.34  , -0.1,  1.1,  -0.2   ,   0.015 ,    0  ,    0 ,    0            , 0  , 0 ,  0          , 0],
         [2, 11 ,-0.34  , -0.1,  1.1,  0.2    ,  0.015  ,   0  ,    0  ,   0             , 0  , 0 ,  0          , 0]]

# spherical joint matrix
# [body1 body2 locinbdy1 locinbdy2]

SphJts =[[ 3, 7,  0.175 , -0.05, -0.022,     0    ,  -0.24 ,  -0.022],
	 [ 4, 7,  0.1096,  0.01,  0.04 ,     0    ,  0.24  ,  0.04 ],
	 [ 5, 8, -0.175 , -0.05, -0.022,     0    ,  -0.24 ,  -0.022],
	 [ 6, 8, -0.1096,  0.01,  0.04 ,     0    ,   0.24 ,   0.04 ],
	 [ 7, 9, -0.02  , -0.07,  0.14 ,     0.67 ,   0    ,   0    ],
	 [ 8, 9,  0.02  , -0.07,  0.14 ,    -0.67 ,   0    ,   0    ]]

	
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
for i in range(0,nb):
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
     
      bx=np.zeros((2,1))
      by=np.zeros((2,1))
      bz=np.zeros((2,1))
      
    else:
      [eigval,eigvec]=LA.eig(I)             #check syntax
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
          lrgdim=max(a,lrgdim)                      
r=JtGeomPar1*lrgdim
re=JtGeomPar2*r
l=JtGeomPar3*r
le=JtGeomPar4*l
lh=JtGeomPar5*l
gap=JtGeomPar6
ang=JtGeomPar7
ncpts=JtGeomPar8
ninsg=JtGeomPar9
noutsg=JtGeomPar10
ngrd=JtGeomPar11
 
if (nsj > 0):
    
    for isj in range(0:nsj):
      SJS=SJScl[isj]
      b2=SphJts[isj][1]
      SphJ=SphJts[isj]
      SphJts=SphJ[5:8]
      pos2=np.reshape(SphJts,(3,1)) 
      t=r*SJS
      if (LA.norm(pos2) == 0):
        bzaxis=[0;0;1]
      else:
        bzaxis=pos2
      yover=0
      while (yover == 0):
          if (bzaxis[0] != 0):
            byaxis[1][0]=1 
            byaxis[2][0]=1
            byaxis[0]=-(bzaxis[1]*byaxis[1]+bzaxis[2]*byaxis[2])/bzaxis[0]
            yover=1
          elif (bzaxis[1] != 0):
            byaxis[0][0]=1 
            byaxis[2][0]=1
            byaxis[1]=-(bzaxis[0]*byaxis[0]+bzaxis[2]*byaxis[2])/bzaxis[1]
            yover=1
          elif (bzaxis[2] !=0):
            byaxis[0][0]=1
            byaxis[1][0]=1
            byaxis[2]=-(bzaxis[0]*byaxis[0]+bzaxis[1]*byaxis[1])/bzaxis[2]
            yover=1
    
  
      bzaxis=bzaxis/LA.norm(bzaxis)
      byaxis=byaxis/LA.norm(byaxis)
      bxaxis[0][0]=byaxis[1]*bzaxis[2]-byaxis[2]*bzaxis[1]
      bxaxis[1][0]=byaxis[2]*bzaxis[0]-byaxis[0]*bzaxis[2]
      bxaxis[2][0]=byaxis[0]*bzaxis[1]-byaxis[1]*bzaxis[0]
  
      R=[bxaxis , byaxis ,bzaxis]

      twopi=2*3.1415
      phis=[0:twopi/(nc-1):twopi];
      thetas=[0:twopi/(ng-1):twopi-twopi/2/(ng-1)];
      z=r*np.cos(phis)
      z1=z(1:nc-1)
      z2=z(2:nc)
      rsphi=r*np.sin(phis)
      bx=[]
      by=[]
      bz=[]
  
      for i in range(1,ng-1):
        x=rsphi*np.cos(thetas(i))
        x1=x(1:nc-1)
        x2=x(2:nc)
        y=rsphi*np.sin(thetas(i))
        y1=y(1:nc-1)
        y2=y(2:nc)
        bx=[bx ,[x1,x2]]
        by=[by ,[y1,y2]]
        bz=[bz ,[z1,z2]]
  

      bxn1=R(1,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[0];
      bxn2=R(1,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[0];
      byn1=R(2,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[1];
      byn2=R(2,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[1];
      bzn1=R(3,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[2];
      bzn2=R(3,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[2];

      bx=[[bxn1 ,0],[ bxn2, pos[0]-r*bzaxis[0]]]
      by=[[byn1 ,0],[ byn2, pos[1]-r*bzaxis[1]]]
      bz=[[bzn1 ,0],[ bzn2, pos[2]-r*bzaxis[2]]]

  
