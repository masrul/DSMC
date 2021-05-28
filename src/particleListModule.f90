submodule (dataTypes) particleListModule
contains
   module subroutine addParticle(self,particle)
      class(dataParticleList)::self
      type(dataParticle)::particle
      integer(8)::ii
      
      ii=self%nPars+1

      self%x(ii)=particle%x
      self%y(ii)=particle%y
      self%z(ii)=particle%z

      self%vx(ii)=particle%vx
      self%vy(ii)=particle%vy
      self%vz(ii)=particle%vz
      
      self%type(ii)=particle%type
      self%index(ii)=particle%index

      self%nPars=self%nPars+1
   end subroutine addParticle

   module subroutine delParticle(self,ParIDx,keepDelting)
      implicit none 
      class(dataparticleList)::self
      integer(8)::ParIDx
      logical::keepDelting
      integer(8)::ii,jj,j

      ii=ParIDx
      jj=ii
      
      keepDelting=.False.
      if(ii==self%nPars)then
         self%nPars=self%nPars-1
         return
      end if 

      do j=self%nPars,ii+1,-1
         if(self%x(j)>-1.d0 .and. self%x(j)<1.d0)then
            jj=j
            self%nPars=self%nPars-1
            keepDelting=.True.
            exit
         else
            self%nPars=self%nPars-1
         end if 
      enddo


      self%x(ii)=self%x(jj)
      self%y(ii)=self%y(jj)
      self%z(ii)=self%z(jj)

      self%vx(ii)=self%vx(jj)
      self%vy(ii)=self%vy(jj)
      self%vz(ii)=self%vz(jj)
      self%type(ii)=self%type(jj)
      self%index(ii)=self%index(jj)


   end subroutine delParticle
end submodule particleListModule
