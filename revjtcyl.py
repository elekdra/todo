#function [cx,cy,cz]=RevJtCyl(pos,drcn,r,re,l,le,lh,nc,ns)

  # Inputs:
  # pos  - position of center of jt in com frame - col vector
  # drcn - direction of joint z axis
  # r    - radius of inner cylinder
  # l    - length of inner cylinder
  # re   - extra radius of  cylinders
  # le   - length of  cylinders
  # lh   - length of handle on the outer cylinders
  # nc   - number of points on full circles
  # ng   - number of longitudinal lines to be drawn
  #
  # Outputs:
  # cx  - x coordinates of set of cylinder segments 
  # cy  - y coordinates of set of cylinder segments 
  # cz  - z coordinates of set of cylinder segments 

  # Finding rotation matrix
def  RevJtCyl(pos,drcn,r,re,l,le,lh,nc,ns):
  if (LA.norm(drcn) == 0):
    print("Revelute jt direction vector is null.  Aborting.")
    
  else:
    czaxis=drcn/LA.norm(drcn)
  

  #Choosing x axis to be towards the CoM of body
  if (LA.norm(pos) != 0):
    cxdir=-pos
  else:
    cxdir[0][0]=random.uniform(0,1)-0.5
    cxdir[1][0]=random.uniform(0,1)-0.5
    cxdir[2][0]=random.uniform(0,1)-0.5
  
  cxdir=cxdir/LA.norm(cxdir);
  
  cyaxis[0][0]=czaxis[1]*cxdir[2]-czaxis[2]*cxdir[1];
  cyaxis[1][0]=czaxis[2]*cxdir[0]-czaxis[0]*cxdir[2];
  cyaxis[2][0]=czaxis[0]*cxdir[1]-czaxis[1]*cxdir[0];

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
      
    
  

  cyaxis=cyaxis/LA.norm(cyaxis)

  cxaxis[0][0]=cyaxis[1]*czaxis[2]-cyaxis[2]*czaxis[1]
  cxaxis[1][0]=cyaxis[2]*czaxis[0]-cyaxis[0]*czaxis[2]
  cxaxis[2][0]=cyaxis[0]*czaxis[1]-cyaxis[1]*czaxis[0]
  
  R=[cxaxis ,cyaxis, czaxis]

  #Axial segments
  twopi=2*3.14159265359
  #thetas=list(range(0,twopi-twopi/2/(ns-1),twopi/(ns-1))
  lby2=l/2
  z1=lby2
  z2=-lby2
  cx=[]
  cy=[]
  cz=[]
  ro=r+re
  zo1=lby2+le
  zo2=lby2
  zo3=-lby2
  zo4=-lby2-le
  
  for i in range (1,ns-1):
    x1=r*np.cos(thetas(i))
    x2=x1
    y1=r*np.sin(thetas(i))
    y2=y1
    cx=[cx ,[[x1],[x2]]]
    cy=[cy ,[[y1],[y2]]]
    cz=[cz ,[[z1],[z2]]]
    xo1=ro*np.cos(thetas(i))
    xo2=xo1
    yo1=ro*np.sin(thetas(i))
    yo2=yo1
    cx=[cx ,[[xo1, xo1],[xo2, xo2]]]
    cy=[cy ,[[yo1, yo1],[yo2, yo2]]]
    cz=[cz ,[[zo1, zo1],[zo2, zo2]]]
  

  #Circles
  #phis=list(range(0,twopi,twopi/(nc-1))
  np=len(phis[1])
  incrcx=r*np.cos(phis)
  incrcx1=incrcx(1:np-1)
  incrcx2=incrcx(2:np)
  incrcy=r*np.sin(phis)
  incrcy1=incrcy(1:np-1)
  incrcy2=incrcy(2:np)
  incrcz1=lby2*np.ones(incrcx1)
  incrcz2=-incrcz1
  
  cx=[cx ,[[incrcx1],[incrcx2]], [[incrcx1], [incrcx2]]]
  cy=[cy ,[[incrcy1],[incrcy2]], [[incrcy1], [incrcy2]]]
  cz=[cz ,[[incrcz1],[incrcz1]], [[incrcz2], [incrcz2]]]
  
  outcrcx=ro*np.cos(phis)
  outcrcx1=outcrcx(1:np-1)
  outcrcx2=outcrcx(2:np)
  outcrcy=ro*np.sin(phis)
  outcrcy1=outcrcy(1:np-1)
  outcrcy2=outcrcy(2:np)
  incrcz1=(lby2+le)*np.ones(incrcx1)
  incrcz2=incrcz1-le
  incrcz3=incrcz2-l
  incrcz4=incrcz3-le
  
  cx=[cx, [outcrcx1, outcrcx2] ,[outcrcx1, outcrcx2]]
  cy=[cy, [outcrcy1, outcrcy2], [outcrcy1, outcrcy2]]
  cz=[cz, [incrcz1, incrcz1] ,[incrcz2, incrcz2]]
  cx=[cx, [outcrcx1, outcrcx2] ,[outcrcx1, outcrcx2]]
  cy=[cy, [outcrcy, outcrcy2] ,[outcrcy1, outcrcy2]]
  cz=[cz, [incrcz3,incrcz3] ,[incrcz4, incrcz4]]

  #Handle
  hx=[[ro, ro],[ ro, ro+lh]]
  hy=[[0 ,0],[ 0, 0]]
  hz=[[lby2, l/5], [l/5, l/5]]
  cx=[cx, hx]
  cy=[cy, hy]
  cz=[cz, hz]

  cxn1=R(1,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[0];
  cxn2=R(1,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[0];
  cyn1=R(2,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[1];
  cyn2=R(2,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[1];
  czn1=R(3,:)*[cx(1,:); cy(1,:); cz(1,:)]+pos[2];
  czn2=R(3,:)*[cx(2,:); cy(2,:); cz(2,:)]+pos[2];

  
  n=len(cxn1[1])

  cx=[[cxn1, 0],[ cxn2 ,cxn2(n)]]
  cy=[[cyn1, 0],[ cyn2 ,cyn2(n)]]
  cz=[[czn1, 0],[ czn2 ,czn2(n)]]
  return cx,cy,cz


