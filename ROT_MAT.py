#function R = Rot_mat(orntn)
def Rot_mat(orntn):
  phi=orntn[0]
  cph=np.cos(phi)
  sph=np.sin(phi)
  theta=orntn[1]
  cth=np.cos(theta)
  sth=np.sin(theta)
  psi=orntn[2]
  cps=np.cos(psi)
  sps=np.sin(psi)
    
  R=[[cps*cph-cth*sph*sps , -sps*cph-cth*sph*cps, sth*sph],
     [ cps*sph+cth*cph*sps, -sps*sph+cth*cph*cps,-sth*cph],
     [sps*sth ,cps*sth,cth]]

  return R


