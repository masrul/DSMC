module coreKernel
   use dataTypes
   use omp_lib


contains 
double precision function myrand()
   implicit none
   call random_number(myrand)
end function myrand

double precision function randVel()
   implicit none 
   double precision::t
   t=max(myrand(),1e-200)
   randVel=inputs%vtemp*sqrt(-log(t))
end function randVel

function Get_randomDir()result(randDir)
   implicit none 
   double precision::randDir(3)
   double precision::B,A,theta
   B=2.d0*myrand()-1.d0
   A=dsqrt(1.d0-B*B)
   theta=myrand()*2.d0*PI

   randDir(1)=B
   randDir(2)=A*dcos(theta)
   randDir(3)=A*dsin(theta)
end function Get_randomDir


subroutine sorter ()
   implicit none

   integer(8)::         &
      iPar,          &
      iCell,         &
      i,j,k,         &
      cellIDx,start, &
      npCount(inputs%nCells),&
      refIDx,idxmax

   double precision::beginTime,endTime

   

   !$omp do  private(i,j,k,iPar,iCell)
   do iPar=1,ParticleList%nPars
      i=min(int(floor((ParticleList%x(iPar)+1.d0)/inputs%dx)),inputs%ni)
      j=min(int(floor((ParticleList%y(iPar)+1.d0)/inputs%dy)),inputs%nj)
      k=min(int(floor((ParticleList%z(iPar)+1.d0)/inputs%dz)),inputs%nk)
      ParticleList%index(iPar)=i*inputs%nj*inputs%nk+j*inputs%nk+k+1
   enddo
   !$omp end do 


   !$omp do private(iCell) 
   do iCell=1,inputs%nCells
      cellInfo%nPars(iCell)=0
   end do 
   !$omp end do 

   !$omp single 
   do iPar=1,ParticleList%nPars
      CellIDx=ParticleList%index(iPar)
      cellInfo%nPars(CellIDx)=cellInfo%nPars(CellIDx)+1
   enddo

   
   start=1
   do iCell=1,inputs%nCells
      cellInfo%startIDx(iCell)=start
      start=start+cellInfo%nPars(iCell)
   enddo

   npCount(:)=0
   do iPar=1,ParticleList%nPars
      CellIDx=ParticleList%index(iPar)
      refIDx=npCount(CellIdx)+cellInfo%startIDx(CellIDx)
      cellInfo%xref(refIDx)=iPar
      npCount(CellIDx)=npCount(CellIDx)+1
   enddo
   !$omp end single 

end subroutine  sorter

subroutine initBoundaries()

   implicit none 

   type(dataParticle)::p

   double precision::&
      cx,cy,cz,      &
      vel(3)

      
   integer(8)::j,k,m,nPars
   double precision::beginTime,endTime
   

   !$omp  do private(j,k,m,vel,p,cx,cy,cz)
   alongY:do j=1,inputs%nj
      alongZ:do k=1,inputs%nk
         cx= -1.d0-inputs%dx
         cy= -1.d0+dble(j-1)*inputs%dy
         cz= -1.d0+dble(k-1)*inputs%dz
         cellLoop:do m=1,inputs%mppc
            p%x=cx+myrand()*inputs%dx
            p%y=cy+myrand()*inputs%dy
            p%z=cz+myrand()*inputs%dz
            vel(:)=Get_randomDir()
            
            vel=vel*randVel()
            p%vx=vel(1)+inputs%vmean
            p%vy=vel(2)
            p%vz=vel(3)

            p%type=0
            p%index=0
            
            !$omp critical 
            call ParticleList%add(p)
            !$omp end critical
         end do cellLoop
      end do alongZ
   end do alongY
   !$omp end do


end subroutine initBoundaries

subroutine moveParticlesWithBCs()
   implicit none 

   type(dataParticle)::pt
   integer(8)::iPar
   double precision::&
      x,y,z,         &
      xt,yt,zt,      &
      t
   double precision::beginTime,endTime


   !$omp  do private(iPar,x,y,z,xt,yt,zt,t,pt)
   do iPar=1,ParticleList%nPars
      x=ParticleList%x(iPar)
      y=ParticleList%y(iPar)
      z=ParticleList%z(iPar)

      xt=x+ParticleList%vx(iPar)*inputs%deltaT
      yt=y+ParticleList%vy(iPar)*inputs%deltaT
      zt=z+ParticleList%vz(iPar)*inputs%deltaT


      if((x<inputs%plate_x .and. xt>inputs%plate_x) .or. &
         (x>inputs%plate_x .and. xt<inputs%plate_x))then

         t=(x-inputs%plate_x)/(x-xt)
         pt%x=x*(1.d0-t)+xt*t
         pt%y=y*(1.d0-t)+yt*t
         pt%z=z*(1.d0-t)+zt*t

         if((pt%y<inputs%plate_dy .and. pt%y> -inputs%plate_dy) .and. &
            (pt%z<inputs%plate_dz .and. pt%z> -inputs%plate_dz))then

            xt=xt-2.d0*(xt-inputs%plate_x)
            ParticleList%vx(iPar)=-ParticleList%vx(iPar)
            ParticleList%type(iPar)=colWithPlate
         end if 
      end if 
   
      if(yt>+1.d0)yt=yt-2.d0
      if(yt<-1.d0)yt=yt+2.d0

      if(zt>+1.d0)zt=zt-2.d0
      if(zt<-1.d0)zt=zt+2.d0

      ParticleList%x(iPar)=xt
      ParticleList%y(iPar)=yt
      ParticleList%z(iPar)=zt
   end do 
   !$omp end  do


end subroutine moveParticlesWithBCs 

subroutine removeOutsideParticles()
   implicit none 
   integer(8)::iPar
   integer::nPars,nremove
   logical::keepDelting
   double precision::beginTime,endTime
   

   

   nPars=particleList%nPars
   nremove=0

   do iPar=1,particleList%nPars
      if(iPar==particleList%nPars)then
         particleList%nPars=particleList%nPars-1
         exit
      end if 
      if (ParticleList%x(iPar)<-1.d0 .or. ParticleList%x(iPar)>1.d0)then
         nremove=nremove+1
         call ParticleList%del(iPar,keepDelting)
         if (.not.keepDelting)exit  
      endif
   end do 
   particleList%nPars=particleList%nPars-1

end subroutine removeOutsideParticles

subroutine collideParticles(nsample)
   
   implicit none 
   integer(8)::nsample

   integer(8)::          &
      iCell,nInstant, &
      nselect,iCol,   &
      iType,jType,    &
      cType,iPar,jPar,&
      nCells,startIDx 
      


   double precision::&
      nMean,aselect, &
      cmax,vrm,      &
      crate,         & 
      vcm(3),vp(3),  &
      dir(3)
   double precision::beginTime,endTime

   

   !$omp do private(iCell, nInstant,nMean,aselect,nselect,iPar,jPar,vrm,vcm,cmax,itype,jtype,ctype,crate,vp,dir,startidx)
   CellLoop:do iCell=1,inputs%nCells
      nMean=dble(sampleInfo%nPars(iCell))/dble(nSample)
      nInstant=cellInfo%nPars(iCell)

      aselect=(                                   &
         (nInstant)*nMean*inputs%pnum*            &
         collisionInfo%maxCollisionRate(iCell)*   &
         (inputs%deltaT/inputs%cellVol) +         &
         collisionInfo%collisionRemainder(iCell)  &
      )
         

      nselect=int(aselect)
      collisionInfo%collisionRemainder(iCell)=aselect-dble(nselect)
      if(nselect>0)then
         if(nInstant<2)then
            collisionInfo%collisionRemainder(iCell)=&
            collisionInfo%collisionRemainder(iCell)+nselect
         else
            cmax=collisionInfo%maxCollisionRate(iCell)
            ColsnLoop:do iCol=1,nselect
               iPar=int(myrand()*nInstant)
               jPar=int(myrand()*nInstant)

               do while (iPar==jPar)
                  jPar=int(myrand()*nInstant)
               enddo 


               startIDx=cellInfo%startIDx(iCell)
               iPar=cellInfo%xref(startIDx+iPar)
               jPar=cellInfo%xref(startIDx+jPar)
               
               vrm=sqrt(                                          &
                  (ParticleList%vx(iPar)-ParticleList%vx(jPar))**2+ &
                  (ParticleList%vy(iPar)-ParticleList%vy(jPar))**2+ &
                  (ParticleList%vz(iPar)-ParticleList%vz(jPar))**2  &
               )  
               crate=inputs%sigmak*vrm
               if(crate>cmax)cmax=crate

               !! Check if these particle actually collide
               if(myrand()<crate/collisionInfo%maxCollisionRate(iCell))then

                  vcm(1)=0.5d0*(ParticleList%vx(iPar)+particleList%vx(jpar))
                  vcm(2)=0.5d0*(ParticleList%vy(iPar)+particleList%vy(jpar))
                  vcm(3)=0.5d0*(ParticleList%vz(iPar)+particleList%vz(jpar))

                  dir=Get_randomDir()
                  vp(1)=dir(1)*vrm
                  vp(2)=dir(2)*vrm
                  vp(3)=dir(3)*vrm

                  !! Adjust particle velocities to reflect collision
                  ParticleList%vx(iPar)=vcm(1)+0.5d0*vp(1)
                  ParticleList%vy(iPar)=vcm(2)+0.5d0*vp(2)
                  ParticleList%vz(iPar)=vcm(3)+0.5d0*vp(3)

                  ParticleList%vx(jPar)=vcm(1)-0.5d0*vp(1)
                  ParticleList%vy(jPar)=vcm(2)-0.5d0*vp(2)
                  ParticleList%vz(jPar)=vcm(3)-0.5d0*vp(3)

                  !! BookKeeping to track particle interactions
                  iType=Particlelist%type(iPar)
                  jType=Particlelist%type(jPar)
                  
                  cType=inflow
                  if((iType+jType)>0)cType=colWithPar

                  ParticleList%type(iPar)=max(cType,iType)
                  ParticleList%type(jPar)=max(cType,jType)
               endif 
            enddo ColsnLoop
            collisionInfo%maxCollisionRate(iCell)=cmax
         endif 
      endif 
   enddo CellLoop
   !$omp end do
   
end subroutine collideParticles


end module coreKernel
