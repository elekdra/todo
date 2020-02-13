#function [cx,cy,cz]=RevJtRings(pos,drcn,r,l,gap,lh,nc,ns)

  # Inputs:
  # pos  - position of center of jt in com frame - col vector
  # drcn - direction of joint z axis
  # r    - radius of inner cylinder
  # l    - length of inner cylinder
  # gap  - fraction of r as radial and axial gaps
  # lh   - length of handle on the outer cylinders
  # nc   - number of points on full circles
  # ns   - number of longitudinal lines to be drawn
  #
  # Outputs:
  # rx  - x coordinates of set of cylinder segments 
  # ry  - y coordinates of set of cylinder segments 
  # rz  - z coordinates of set of cylinder segments 

  # Finding rotation matrix
def RevJtRings(pos,drcn,r,l,gap,lh,nc,ns):  
  if (LA.norm(drcn) == 0):
    print("Revelute jt direction vector is null.  Aborting.")
    
  else:
    czaxis=drcn/(LA.norm(drcn))
  

  #Choosing x axis to be towards the CoM of body
  if (LA.norm(pos) != 0):
    cxdir=-pos
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
  #thetas=[0:twopi/(ns-1):twopi-twopi/2/(ns-1)]; # check syntax
  lby2=l/2
  g=r*gap
  z1=lby2-g
  z2=-lby2+g
  cx=[] 
  cy=[]
  cz=[]
  rr=r+g
  
  for i in range(1,ns-1):
    x1=rr*np.cos(thetas(i))
    x2=x1
    y1=rr*np.sin(thetas(i))
    y2=y1
    cx=[cx ,[[x1],[x2]]]  #transpose 
    cy=[cy ,[[y1],[y2]]]
    cz=[cz ,[[z1],[z2]]]
  

  #Circles
  phis=[0:twopi:twopi/(nc-1)]  #list creation   #check syntax
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
  cxn1=R(1,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[0]
  cxn2=R(1,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[0]
  cyn1=R(2,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[1]
  cyn2=R(2,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[1]
  czn1=R(3,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[2]
  czn2=R(3,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[2]
  n=len(cxn1[1])
  cx=[[cxn1, 0],[cxn2 , cxn2(n)]]
  cy=[[cyn1, 0],[cyn2 , cyn2(n)]]
  cz=[[czn1, 0],[czn2 , czn2(n)]]
  return cx,cy,cz
  

  




