# input function
# bodies =
# [glbl_pos_of_com glbl_orntn_of_bdyframe loc_of_com_in_bdyframe 
# mass Ixx_about_com Iyy Yzz Ixy Iyz Izx]
import pygame
from pygame.locals import *
import random
import time
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np 
from numpy import linalg as LA
             
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


#t=0;
#BdyMn=read("q.txt");
#BdyMn=read("straight_const_vel_q.dat",40,132);
#t=read("straight_const_vel_t.dat",40,1);
#BdyMn=read("leftturn_const_vel_const_rad_q.dat",501,132);
#t=read("leftturn_const_vel_const_rad_t.dat",501,1);



BdyMn=open("double_lane_change_q.dat", "r")
if BdyMn.mode == 'r':
	BdyMn =BdyMn.read()	
t=open("double_lane_change_t.dat","r")
if t.mode == 'r':
	t =t.read()
	#print(t[1][0])

	
# Revolute joint parameters
JtGeomPar1=0.1 #used in r=JtGeomPar(1)*lrgdim
JtGeomPar2=0.3 #used in re=JtGeomPar(2)*r
JtGeomPar3=1   #used in l=JtGeomPar(3)*r
JtGeomPar4=0.3 #used in le=JtGeomPar(4)*l
JtGeomPar5=0.3 #used in lh=JtGeomPar(5)*l
JtGeomPar6=0.1 #used in gap=JtGeomPar(6)*r
  
# Spherical joint parameters
JtGeomPar7=60 #used in ang=JtGeomPar(7)
  
JtGeomPar8=24 #used in ncpts=JtGeomPar(8)
JtGeomPar9=24 #used in ninsg=JtGeomPar(9)
JtGeomPar10=12#used in noutsg=JtGeomPar(10)
JtGeomPar11=10#used in ngrd=JtGeomPar(11)


RJScl=0.1*np.ones((6,1))
SJScl=0.1*np.ones((6,1))
BScl=0.5*np.ones((11,1))
BClrs=2*np.ones((11,1))

#function [bx,by,bz] = RightCuboid(m,I,scl)
  
  # To return the body edges
  #
  # Inputs:
  # m   - mass in Kg
  # I   - inertia tensor about CoM
  # scl - scale to be applied (positive number)
  #
  # Outputs:
  # bx  - x coordinates of edges of right cuboid body
  # by  - y coordinates of edges of right cuboid body
  # bz  - z coordinates of edges of right cuboid body
def RightCuboid(m,I,scl):  
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
  
  s1=0.5*math.sqrt(6*(e2+e3-e1)/m)
  s2=0.5*math.sqrt(6*(e3+e1-e2)/m)
  s3=0.5*math.sqrt(6*(e1+e2-e3)/m)
  
  bx=[[s1,s1,s1,s1,s1,s1,s1,s1,-s1,-s1,-s1,-s1],[-s1,-s1,-s1,-s1,s1,s1,s1,s1,-s1,-s1,-s1,-s1]]
  by=[[s2,-s2,-s2,s2,s2,-s2,-s2,s2,s2,-s2,-s2,s2],[s2,-s2,-s2,s2,-s2,-s2,s2,s2,-s2,-s2,s2,s2]]
  bz=[[s3,s3,-s3,-s3,s3,s3,-s3,-s3,s3,s3,-s3,-s3],[s3,s3,-s3,-s3,s3,-s3,-s3,s3,s3,-s3,-s3,s3]]
  bx=scl[i]*bx
  by=scl[i]*by
  bz=scl[i]*bz
  R = np.zeros(shape=(3,3))
  R=[v1,v2,v3]
  n=len(bx[1])  # syntax
  r=[[bx[0], bx[1]], [by[0], by[1]], [bz[0], bz[1]]]
  r=np.reshape(r,(3,24))
  xyz=np.matmul(R,r)
  bx=np.reshape(xyz[0],(2,12))
  by=np.reshape(xyz[1],(2,12))
  bz=np.reshape(xyz[2],(2,12))
  return bx,by,bz

#RevJtRings();
#RevJtCyl();
#SphJtRings();
#SphJtBall();
#Rot_mat();
#function [wndow] = RBSanim(bodies,RevJts,SphJts,t,BdyMn,JtGeomPar,RJScl,SJScl,BScl,BClrs)
  
  # To animate the given mechanism


def main():
    pygame.init()
    display = (800,600)
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
    gluPerspective(45, (display[0]/display[1]), .1 , 50.0)
    glTranslatef(0.0,0.0, -5)
    temp=0
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        RBSanim(bodies,RevJts,SphJts,t,BdyMn,JtGeomPar,RJScl,SJScl,BScl,BClrs);
        pygame.display.flip()
        pygame.time.wait(10)    
main()


