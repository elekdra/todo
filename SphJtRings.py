#function [bx,by,bz]=SphJtRings(pos,r,ang,gap,nc,ns)

  # Inputs:
  # pos - position of center of jt in CoM frame - col vector
  # r   - ball radius
  # ang - angle in deg at which rings are to be drawn
  #       0 is equator and 90 is pole
  # gap - gap between rings and ball in fraction of radius
  # nc  - number of points on full circles
  # ng  - number of circlular segments
  #
  # Outputs:
  # bx  - x coordinates of set of ball segments 
  # by  - y coordinates of set of ball segments 
  # bz  - z coordinates of set of ball segments 

  # Finding rotation matrix to set up the ball points
def SphJtRings(pos,r,ang,gap,nc,ns):

  rl=r*(1+gap)
  
  if (LA.norm(pos) == 0):
    bzaxis[0][0]=random.uniform(0,1)-0.5
    bzaxis[1][0]=random.uniform(0,1)-0.5
    bzaxis[2][0]=random.uniform(0,1)-0.5
  else:
    bzaxis=pos

  yover=0
  byaxis[0][0]=random.uniform(0,1)
  byaxis[1][0]=random.uniform(0,1)
  byaxis[2][0]=random.uniform(0,1)
  while (yover <= 0):
    if (bzaxis[0][0] != 0):
      byaxis[0][0]=-(bzaxis[1][0]*byaxis[1][0]+bzaxis[2][0]*byaxis[2][0])/bzaxis[0][0]
      yover=1
    elif (bzaxis[1][0] != 0):
      byaxis[1][0]=-(bzaxis[0][0]*byaxis[0][0]+bzaxis[2][0]*byaxis[2][0])/bzaxis[1][0]
      yover=1
    elif (bzaxis[2][0] !=0):
      byaxis[2][0]=-(bzaxis[0][0]*byaxis[0][0]+bzaxis[1][0]*byaxis[1][0])/bzaxis[2][0]
      yover=1
 
  bzaxis=bzaxis/(LA.norm(bzaxis))
  byaxis=byaxis/(LA.norm(byaxis))
  bxaxis[0][0]=byaxis[1][0]*bzaxis[2][0]-byaxis[2][0]*bzaxis[1][0]
  bxaxis[1][0]=byaxis[2][0]*bzaxis[0][0]-byaxis[0][0]*bzaxis[2][0]
  bxaxis[2][0]=byaxis[0][0]*bzaxis[1][0]-byaxis[1][0]*bzaxis[0][0]
  
  R=[bxaxis, byaxis, bzaxis]
    
  twopi=2*(22/7)
  angr=ang*(22/7)/180
  phis=list(range(angr,(22/7)-angr,twopi/(nc-1))) 
  np=len(phis[1])
  if (phis(np) < (22/7)-angr): 
    np=np+1


  
  thetas=list(range(0,twopi-twopi/2/(ns-1),twopi/(ns-1)))#list
  z=rl*(np.cos(phis))
  z1=z[0:np-1]
  z2=z[1:np] 
  rsphi=rl*(np.sin(phis))
  bx=[]
  by=[]
  bz=[]
  rthetas=list(range(0,twopi,twopi/(nc-1)))
  nrt=len(rthetas[1])
  
  for i in range(1,ns-1):
    x=rsphi*np.cos(thetas(i))
    x1=x[0:np-1]
    x2=x[1:np]
    y=rsphi*sin(thetas(i))
    y1=y[0:np-1] ##
    y2=y[1:np]  ##
    bx=[bx ,[[x1],[x2]]]
    by=[by ,[[y1],[y2]]]
    bz=[bz ,[[z1],[z2]]]
  
  #Constructing the two rings
  ringx=rl*np.sin(angr)*np.cos(rthetas)
  ringy=rl*np.sin(angr)*np.sin(rthetas)
  ringzu=rl*np.cos(angr)*np.ones(ringx)
  ringzl=-ringzu
  
  bx=[bx ,[[ringx[0:nrt-1], ringx[0:nrt-1],[ ringx[1:nrt], ringx[1:nrt]]]]]  
  by=[by ,[[ringy[0:nrt-1], ringy[0:nrt-1],[ ringy[1:nrt], ringy[1:nrt]]]]]
  bz=[bz ,[[ringzu[0:nrt-1], ringzl[0:nrt-1],[ ringzu[1:nrt], ringzl[1:nrt]]]]]

  #Constructing the connection of lower ring to stem
  cphis=list(range(%pi-angr,22/7,twopi/(nc-1)))
  ncp=len(cphis[2])
  if (cphis(ncp) < (22/7)):
    cphis=[cphis ,(22/7)] # ensuring that cphis go upto %pi
    ncp=ncp+1
  
  
  ncon=len(cphis[2])
  conz=rl*np.cos(cphis)
  conx=rl*np.sin(cphis)
  cony=np.zeros(conx)
  
  bx=[bx ,[conx[0:ncon-1],conx[1:ncon]]]
  by=[by ,[cony[0:ncon-1],cony[1:ncon]]]
  bz=[bz ,[conz[0:ncon-1],conz[1:ncon]]]
  
  bxn1=R[0]*[bx[0], by[0], bz[0]]+pos[0]
  bxn2=R[0]*[bx[1], by[1], bz[1]]+pos[0]
  byn1=R[1]*[bx[0], by[0], bz[0]]+pos[1]
  byn2=R[1]*[bx[1], by[1], bz[1]]+pos[1]
  bzn1=R[2]*[bx[0], by[0], bz[0]]+pos[2]
  bzn2=R[2]*[bx[1], by[1], bz[1]]+pos[2]

  bx=[[bxn1, 0] ,[bxn2 ,pos[0]-rl*bzaxis[0]]]
  by=[[byn1, 0] ,[byn2 ,pos[1]-rl*bzaxis[1]]]
  bz=[[bzn1, 0] ,[bzn2 ,pos[3]-rl*bzaxis[3]]]

  return bx,by,bz
