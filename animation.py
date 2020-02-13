#function [wndow] = RBSanim(bodies,RevJts,SphJts,t,BdyMn,JtGeomPar,RJScl,SJScl,BScl,BClrs)
  
  # To animate the given mechanism
def RBSanim(bodies,RevJts,SphJts,t,BdyMn,JtGeomPar,RJScl,SJScl,BScl,BClrs):  
  nb=len(bodies)
  nrj=len(RevJts)
  nsj=len(SphJts)
  
  # Creating bodies
  
  lrgdim=0
  
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
      bx.transpose()=[ 0 , 0 ]
      bx=bx.transpose()
      by.transpose()=[ 0 , 0 ]
      by.transpose()
      bz.transpose()=[ 0 , 0 ]
      bz.transpose()
    else:
      [bx,by,bz]=RightCuboid(m,I,BScl[i]);
    

    #Finding the largest body dimension in order to size the joints
    
    if (BScl[i] > 0):
      lrgdim=max(lrgdim,max([bx by bz]/BScl[i]))
    
    
    x[i]=x+string[i]
    y[i]=y+string[i]
    z[i]=z+string[i]

    execstr(strcat([x(i) "=bx"]));
    execstr(strcat([y(i) "=by"]));
    execstr(strcat([z(i) "=bz"]));

  

  # Revolute joint parameters
  r=JtGeomPar[1]*lrgdim
  re=JtGeomPar[2]*r
  l=JtGeomPar[3]*r
  le=JtGeomPar[4]*l
  lh=JtGeomPar[5]*l
  gap=JtGeomPar[6]
  
  # Spherical joint parameters
  ang=JtGeomPar[7]
  
  ncpts=JtGeomPar[8]
  ninsg=JtGeomPar[9]
  noutsg=JtGeomPar[10]
  ngrd=JtGeomPar[11]

  # Constructing Revolute joints
  
  if (nrj > 0):

    for irj in range(1,nrj):
 
      RJS=RJScl[irj]
      
      b1=RevJts[irj][1]
      pos1=RevJts(irj,3:5)'-bodies(b1,7:9)'
      phi1=RevJts[irj][9]
      th1=RevJts[irj][10]
      drcn1=[[[np.sin(th1)*np.cos(phi1)],[np.sin(th1)*np.sin(phi1)]], np.cos(th1)]
      #pause
      [bx,by,bz]=RevJtRings(pos1,drcn1,r*RJS,l*RJS,gap,lh*RJS,ncpts,noutsg)
      execstr(strcat([x(b1) "=[bx " x(b1) "]"]))
      execstr(strcat([y(b1) "=[by " y(b1) "]"]))
      execstr(strcat([z(b1) "=[bz " z(b1) "]"]))
      
      b2=RevJts[irj][2]
      pos2=RevJts(irj,6:8)'-bodies(b2,7:9)'
      phi2=RevJts[irj][12]
      th2=RevJts[irj][13]
      drcn2=[np.sin(th2)*np.cos(phi2) , np.sin(th2)*np.sin(phi2),np.cos(th2)]
      #pause
      [bx,by,bz]=RevJtCyl(pos2,drcn2,r*RJS,re*RJS,l*RJS,le*RJS,lh*RJS,ncpts,ninsg);
      execstr(strcat([x(b2) "=[bx " x(b2) "]"]))
      execstr(strcat([y(b2) "=[by " y(b2) "]"]))
      execstr(strcat([z(b2) "=[bz " z(b2) "]"]))
      
    
    
  

  # Constructing Spherical joints

  if (nsj > 0):
    
    for isj in range(1,nsj):

      SJS=SJScl[isj]
      
      b1=SphJts[isj[1]]
      pos1=SphJts(isj,3:5)'-bodies(b1,7:9)';
      #pause
      [bx,by,bz]=SphJtRings(pos1,r*SJS,ang,gap,ncpts,noutsg);
      execstr(strcat([x(b1) "=[bx " x(b1) "]"]));
      execstr(strcat([y(b1) "=[by " y(b1) "]"]));
      execstr(strcat([z(b1) "=[bz " z(b1) "]"]));
      
      b2=SphJts(isj,2);
      pos2=SphJts(isj,6:8)'-bodies(b2,7:9)';
      #pause
      [bx,by,bz]=SphJtBall(pos2,r*SJS,ncpts,ninsg);
      execstr(strcat([x(b2) "=[bx " x(b2) "]"]));
      execstr(strcat([y(b2) "=[by " y(b2) "]"]));
      execstr(strcat([z(b2) "=[bz " z(b2) "]"]));
      
  

  // Finding the extent of ground to be drawn and creating grid for that
  for i in range(2,nb):
    execstr(strcat(["xyz=[" x(i) " " y(i) " " z(i) "]"]))
    mxdim(i)=max(xyz)
    mnxs(i)=min(BdyMn(:,(i-1)*6+1))-mxdim(i)
    mxxs(i)=max(BdyMn(:,(i-1)*6+1))+mxdim(i)
    mnys(i)=min(BdyMn(:,(i-1)*6+2))-mxdim(i)
    mxys(i)=max(BdyMn(:,(i-1)*6+2))+mxdim(i)
    mnzs(i)=min(BdyMn(:,(i-1)*6+3))-mxdim(i)
    mxzs(i)=max(BdyMn(:,(i-1)*6+3))+mxdim(i)
  
  
  mnx=min(mnxs) 
  mxx=max(mxxs)
  mny=min(mnys)
  mxy=max(mxys)
  mnz=min(mnzs) 
  mxz=max(mxzs)

  grdx=[mnx:(mxx-mnx)/(ngrd-1):mxx]
  grdy=[mny:(mxy-mny)/(ngrd-1):mxy]
  x0=[mnx*np.ones(1,ngrd) grdx; mxx*np.ones(1,ngrd) grdx]
  y0=[grdy mny*np.ones(1,ngrd); grdy mxy*np.ones(1,ngrd)]
  z0=np.zeros(x0)
  
  # Animating
  
  mhndle=scf(100001)
  
  nq=len(BdyMn[1])

  loopout=0
  
  while (loopout == 0):
  
    clf(mhndle,"reset");

    a=gca()
    mhndle.immediate_drawing="on"

    realtimeinit(0.001)

    for i in range(1,1):

      realtime(i)
      mhndle.immediate_drawing="off"

      for j in range(1,nb):
      
        orntn=BdyMn(i,(j-1)*6+4:(j-1)*6+6)
//        if (j ==1):
//          orntn=[0 0 0]
//        
        R=Rot_mat(orntn);
      
        execstr(strcat(["xyz1=[" x(j) "(1,:);" y(j) "(1,:);" z(j) "(1,:)]"]))
        execstr(strcat(["xyz2=[" x(j) "(2,:);" y(j) "(2,:);" z(j) "(2,:)]"]))
        nxyz1=R*xyz1; nxyz2=R*xyz2

        posn=BdyMn(i,(j-1)*6+1:(j-1)*6+3)
        nx=[nxyz1(1,:); nxyz2(1,:)]+posn(1)
        ny=[nxyz1(2,:); nxyz2(2,:)]+posn(2)
        nz=[nxyz1(3,:); nxyz2(3,:)]+posn(3)
      
        plot3d(nx,ny,nz);
        p(j)=strcat(["p" string(j)]);
        execstr(strcat([p(j) "=gce()"]));
        execstr(strcat([p(j) ".thickness=2"]));
        execstr(strcat([p(j) ".foreground=BClrs(j)"]));

      end
      
      plot3d(x0,y0,z0);
      p0=gce()
      p0.thickness=2
      p0.foreground=BClrs(1);
      
      mhndle.immediate_drawing="on";

    end
    
    a.isoview="on";
    a.data_bounds=[mnx mny mnz; mxx mxy mxz];

    disp("If you want to change orientation do that now and then enter return.")
    pause
  
    for i=2:10:nq

      realtime(i)
      mhndle.immediate_drawing="off";

      for j=1:nb
      
        orntn=BdyMn(i,(j-1)*6+4:(j-1)*6+6);
        if (j == 1)
          orntn=[0 0 0];
        end
        R=Rot_mat(orntn);
      
        execstr(strcat(["xyz1=[" x(j) "(1,:);" y(j) "(1,:);" z(j) "(1,:)]"]));
        execstr(strcat(["xyz2=[" x(j) "(2,:);" y(j) "(2,:);" z(j) "(2,:)]"]));
        nxyz1=R*xyz1; nxyz2=R*xyz2;
  
        posn=BdyMn(i,(j-1)*6+1:(j-1)*6+3);
        nx=[nxyz1(1,:); nxyz2(1,:)]+posn(1);
        ny=[nxyz1(2,:); nxyz2(2,:)]+posn(2);
        nz=[nxyz1(3,:); nxyz2(3,:)]+posn(3);
      
        execstr(strcat([p(j) ".data.x=nx"]));
        execstr(strcat([p(j) ".data.y=ny"]));
        execstr(strcat([p(j) ".data.z=nz"]));
  
      
      end

      p0.data.x=x0;
      p0.data.y=y0;
      p0.data.z=z0;
      
      mhndle.immediate_drawing="on";

    end

    loopout=input("Enter 0 to contune and 1 to stop animation.")
    
  end
  

  
  wndow=[];
  
endfunction



