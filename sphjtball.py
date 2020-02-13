import numpy as np
from numpy import linalg as LA 
import random

def SphJtBall(pos,r,nc,ng):
  
  if (LA.norm(pos) == 0):
    bzaxis[0][0]=random.uniform(0,1)-0.5
    bzaxis[1][0]=random.uniform(0,1)-0.5
    bzaxis[2][0]=random.uniform(0,1)-0.5
  else:
    bzaxis=pos
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
  
  R=[bxaxis, byaxis ,bzaxis]

  twopi=2*3.1415
  #phis=[0:twopi/(nc-1):twopi];
  #thetas=[0:twopi/(ng-1):twopi-twopi/2/(ng-1)];
  z=r*np.cos(phis)
  #z1=z(1:nc-1)
  #z2=z(2:nc)
  rsphi=r*np.sin(phis)
  bx=[]
  by=[]
  bz=[]
  
  for i in range(1,ng-1):
    x=rsphi*np.cos(thetas(i))
    #x1=x(1:nc-1)
    #x2=x(2:nc)
    y=rsphi*np.sin(thetas(i))
    #y1=y(1:nc-1)
    #y2=y(2:nc)
    bx=[bx ,[x1,x2]]
    by=[by ,[y1,y2]]
    bz=[bz ,[z1,z2]]
  

  #bxn1=R(1,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[0];
  #bxn2=R(1,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[0];
  #byn1=R(2,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[1];
  #byn2=R(2,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[1];
  #bzn1=R(3,:)*[bx(1,:); by(1,:); bz(1,:)]+pos[2];
  #bzn2=R(3,:)*[bx(2,:); by(2,:); bz(2,:)]+pos[2];

  bx=[[bxn1 ,0],[ bxn2, pos[0]-r*bzaxis[0]]]
  by=[[byn1 ,0],[ byn2, pos[1]-r*bzaxis[1]]]
  bz=[[bzn1 ,0],[ bzn2, pos[2]-r*bzaxis[2]]]

  return bx,by,bz

