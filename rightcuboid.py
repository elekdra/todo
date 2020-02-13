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
