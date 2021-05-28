submodule (dataTypes) sampleKernel
   use omp_lib

contains

module subroutine sampleParticles(self)
   implicit none 
   class(dataSampleInfo)::self
   integer(8)::ii,iPar
   double precision::beginTime,endTime

 

   do iPar=1,particleList%nPars
      ii=ParticleList%index(iPar) 
      self%nPars(ii)=self%nPars(ii)+1
      self%vx(ii)=self%vx(ii)+particleList%vx(iPar)
      self%vy(ii)=self%vy(ii)+particleList%vy(iPar)
      self%vz(ii)=self%vz(ii)+particleList%vz(iPar)
      self%energy(ii)=self%energy(ii)+0.5*(&
         particleList%vx(ii)**2+&
         particleList%vy(ii)**2+&
         particleList%vz(ii)**2 &
      )
           
   enddo

end subroutine sampleParticles

module subroutine initSample(self)
   implicit none 
   class(dataSampleInfo)::self
   integer(8)::ii,iPar
   self%nsample=0
   do iPar=1,particleList%nPars
      ii=ParticleList%index(iPar) 
      sampleInfo%nPars(ii)=0
      sampleInfo%vx(ii)=0.d0
      sampleInfo%vy(ii)=0.d0
      sampleInfo%vz(ii)=0.d0
      sampleInfo%energy(ii)=sampleInfo%energy(ii)+0.5*(&
         particleList%vx(ii)**2+&
         particleList%vy(ii)**2+&
         particleList%vz(ii)**2 &
      )
           
   enddo
end subroutine initSample

end submodule sampleKernel 
