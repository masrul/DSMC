program main
   use dataTypes
   use coreKernel
   use omp_lib
   
   implicit none 
    
   integer(8)::iStep,nsample,i,iPar
   integer(8)::sample_reset,iCell
   double precision::time


   call inputs%parse()
   call inputs%initDefault(   &
      ni=64,nj=64,nk=64,      &
      mppc=10,                &
      mach=20.d0,             &
      density=1d30,           &
      plate_x=-0.25d0,          &
      plate_dy=0.25d0,          &
      plate_dz=0.5d0,           &
      time=-1.d0                &
   )


   call inputs%process()
   open(unit=100,file='particles.dat')
   print*,"Input parameters:"
   print*,"ni,nj,nk: ",inputs%ni,inputs%nj,inputs%nk
   print*,"mppc: ",inputs%mppc
   print*,"mach: ",inputs%mach
   print*,"density: ",inputs%density
   print*,"plate_x: ",inputs%plate_x
   print*,"plate_dy: ",inputs%plate_dy
   print*,"plate: ",inputs%plate_dz
   print*,"time: ",inputs%time
   print*,"sigmak: ",inputs%sigmak
   print*,"pnum: ",inputs%pnum
   print*,"vmean: ",inputs%vmean
   print*,"nCells: ",inputs%nCells
   print*,"nSteps: ",inputs%nSteps
   print*,"dx: ",inputs%dx
   print*,"dy: ",inputs%dy
   print*,"dz: ",inputs%dz
   print*,"cellVol: ",inputs%cellVol
   print*,"vtemp: ",inputs%vtemp
   print*,"deltaT: ",inputs%deltaT

   call particleList%malloc()
   call collisionInfo%malloc()
   call cellInfo%malloc()
   call sampleInfo%malloc()

   

   write(*,*)"time= ",inputs%time," ,nsteps = ",inputs%nSteps

   sampleInfo%nsample=0

   !! re-sample 4 time during simulation 
   sample_reset=inputs%nSteps/4

   !! Begin simulation. Initialize collision data
   call collisionInfo%init()

   !$omp parallel default(shared) private(myid)
   myid=omp_get_thread_num()


   !$omp single 
   if(omp_get_thread_num()==0)then
      print*,'Number of OMP thread in work: ', omp_get_num_threads()
   endif 
   !$omp end single 

   !$omp end  parallel  

   do iStep=1,inputs%nSteps
      
      !$omp parallel default(shared) private(myid)

       
      call initBoundaries()
      
      call moveParticlesWithBCs()

      !$omp single 
      call removeOutsideParticles()
      !$omp end single
      

      call sorter()

      !$omp single 
      if(mod((iStep-1),sample_reset)==0)then
        
         call SampleInfo%init()
         nsample=0
      endif

      nsample=nsample+1
      
      call sampleInfo%sample()
      !$omp end single
      
      call collideParticles(nsample)


      !$omp single
      if(mod((iStep-1),15)==0)then
         print*,'iStep:',iStep,'  ','nParticles: ',particleList%nPars
      end if 
      !$omp end single

      !$omp end parallel 
   end do

   !$omp single
   do iPar=1,ParticleList%nPars
      write(100,'(3F10.5,I5)')ParticleList%x(iPar),ParticleList%y(iPar),ParticleList%z(iPar),particleList%type(iPar)
   end do 
   !$omp  end single 
end program main
