submodule (dataTypes) allocateKernel

contains

module subroutine allocateCollisionInfo(self)
   implicit none 
   class(dataCollisionInfo)::self
   integer(8)::nCells
   integer(8)::ierr(2)

   nCells=inputs%nCells
   
   allocate(self%maxCollisionRate(nCells),stat=ierr(1))
   allocate(self%collisionRemainder(nCells),stat=ierr(2))

   if(any(ierr>0))then
      write(*,*)'Catastrophic Error: Insufficient memory, decrease grid number'
      stop
   endif
   self%maxCollisionRate(:)=0.d0
   self%collisionRemainder(:)=0.d0

end subroutine allocateCollisionInfo

module subroutine allocateCellInfo(self)
   class(dataCellInfo)::self
   integer(8)::nCells
   integer(8)::ierr(3)

   nCells=inputs%nCells
   allocate(self%xref(MAX_PARTICLE),stat=ierr(1))
   allocate(self%nPars(nCells),stat=ierr(2))
   allocate(self%startIDx(nCells),stat=ierr(3))

   if(any(ierr>0))then
      write(*,*)'Catastrophic Error: Insufficient memory, decrease grid number'
      stop
   endif

   self%xref(:)=0
   self%nPars(:)=0
   self%startIDx(:)=0
end subroutine allocateCellInfo 

module subroutine allocateSampleInfo(self)
   class(dataSampleInfo)::self
   integer(8)::nCells
   integer(8)::ierr(5)

   nCells=inputs%nCells
   
   allocate(self%nPars(nCells),stat=ierr(1))
   allocate(self%vx(nCells),stat=ierr(2))
   allocate(self%vy(nCells),stat=ierr(3))
   allocate(self%vz(nCells),stat=ierr(4))
   allocate(self%energy(nCells),stat=ierr(5))

   if(any(ierr>0))then
      write(*,*)'Catastrophic Error: Insufficient memory, decrease grid number'
      stop
   endif
   self%nPars(:)=0
   self%vx(:)=0.d0
   self%vy(:)=0.d0
   self%vz(:)=0.d0
   self%energy(:)=0.d0
end subroutine allocateSampleInfo

module subroutine allocateParticleList(self)
   class(dataParticleList)::self
   integer(8)::ierr(8)
   self%nPars=0
   allocate(self%x(MAX_PARTICLE),stat=ierr(1))
   allocate(self%y(MAX_PARTICLE),stat=ierr(2))
   allocate(self%z(MAX_PARTICLE),stat=ierr(3))
   allocate(self%vx(MAX_PARTICLE),stat=ierr(4))
   allocate(self%vy(MAX_PARTICLE),stat=ierr(5))
   allocate(self%vz(MAX_PARTICLE),stat=ierr(6))
   allocate(self%type(MAX_PARTICLE),stat=ierr(7))
   allocate(self%index(MAX_PARTICLE),stat=ierr(8))

   if(any(ierr>0))then
      write(*,*)'Catastrophic Error: Insufficient memory, decrease MAX_PARTICLE in dataTypes.f90'
      stop
   endif

end subroutine allocateParticleList
end submodule allocateKernel
