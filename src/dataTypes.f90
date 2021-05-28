module dataTypes
   integer(8),parameter::inflow=0
   integer(8),parameter::colWithPar=1
   integer(8),parameter::colWithPlate=2
   integer(8),parameter::MAX_PARTICLE=10000000
   double  precision::PI=3.14159265358979323
  
   type dataTimer
      double precision::initBoundaries=0.d0
      double precision::moveParticles=0.d0
      double precision::removeParticles=0.d0
      double precision::sorter=0.d0
      double precision::sampling=0.d0
      double precision::collision=0.d0
   end type dataTimer


   type dataParticle
      double precision::x,y,z
      double precision::vx,vy,vz
      integer(8)::type
      integer(8)::index
   end type dataParticle

   type dataSampleInfo
      integer(8),allocatable::nPars(:)
      double precision,allocatable,dimension(:)::vx,vy,vz
      double precision,allocatable,dimension(:):: energy
      integer(8)::nsample

      contains
         procedure::malloc=>allocateSampleInfo
         procedure::sample=>sampleParticles
         procedure::init=>initSample
   endtype dataSampleInfo

   type dataCellInfo
      integer(8),allocatable::xref(:)
      integer(8),allocatable::nPars(:) 
      integer(8),allocatable::startIDx(:)    
      contains
        procedure::malloc=>allocateCellInfo 
   end type dataCellInfo

   type dataCollisionInfo
      double precision,allocatable:: maxCollisionRate(:)
      double precision,allocatable:: collisionRemainder(:)
      contains
         procedure::malloc=>allocateCollisionInfo
         procedure::init=>initCollisionInfo
   endtype


   type dataParticleList
      double precision,dimension(:),allocatable::x,y,z
      double precision,dimension(:),allocatable::vx,vy,vz
      integer(8),dimension(:),allocatable::type
      integer(8),dimension(:),allocatable::index
      integer(8)::nPars

      contains
         procedure::add=>addParticle
         procedure::del=>delParticle
         procedure::malloc=>allocateParticleList
   end type dataParticleList


   type dataInParams
      integer(8)::ni=-1,nj=-1,nk=-1
      integer(8)::mppc=-1
      double precision::mach=-1.d0
      double precision::density=-1.d0
      double precision::plate_x=-1.d0
      double precision::plate_dy=-1.d0
      double precision::plate_dz=-1.d0
      double precision::time=-1.d0
      double precision::sigmak=1d-28
      double precision::pnum=-1.d0
      double precision::vmean=-1.d0
      integer(8)::nCells
      integer(8)::nSteps
      double precision::dx=-1.d0
      double precision::dy=-1.d0
      double precision::dz=-1.d0
      double precision::cellVol=-1.d0
      double precision::vtemp=-1.d0
      double precision::deltaT=-1.d0

      contains
         procedure,pass(self)::initDefault=>initDefaultParams
         procedure::parse=>parseInParams
         procedure::process=>processInParams
   end type dataInParams

   interface

      module subroutine initDefaultParams( &
         self,ni,nj,nk,             &
         mppc,mach,density,         &
         plate_x,plate_dy,plate_dz, &
         time                       &
         )

         class(dataInParams)::self
         integer(8)::ni,nj,nk
         integer(8)::mppc
         double precision::mach
         double precision::density
         double precision::plate_x
         double precision::plate_dy
         double precision::plate_dz
         double precision::time
      end subroutine initDefaultParams

      module subroutine parseInParams(self)
          class(dataInParams)::self
      end subroutine parseInParams

      module subroutine addParticle(self,particle)
         class(dataParticleList)::self
         type(dataParticle)::particle
      end subroutine addParticle

      module subroutine delParticle(self,ParIDx,keepDelting)
         class(dataparticleList)::self
         integer(8)::ParIDx
         logical::keepDelting
      end subroutine delParticle


      module subroutine allocateCellInfo(self)
         class(dataCellInfo)::self
      end subroutine allocateCellInfo 

      module subroutine allocateSampleInfo(self)
         class(dataSampleInfo)::self
      end subroutine allocateSampleInfo

      module subroutine allocateParticleList(self)
         class(dataParticleList)::self
      end subroutine allocateParticleList
      
      module subroutine allocateCollisionInfo(self)
         class(dataCollisionInfo)::self
      end subroutine allocateCollisionInfo
      
      module subroutine sampleParticles(self)
         class(dataSampleInfo)::self
      end subroutine sampleParticles

      module subroutine initSample(self)
         class(dataSampleInfo)::self
         integer(8)::nsample
      end subroutine initSample

      module subroutine initCollisionInfo(self)
         class(dataCollisionInfo)::self
      end subroutine initCollisionInfo

      module subroutine processInParams(self)
         class(dataInParams)::self
      end subroutine processInParams
   endinterface

   type(dataInParams)::inputs
   type(dataCellInfo)::cellInfo
   type(dataSampleInfo)::SampleInfo
   type(dataCollisionInfo)::collisionInfo
   type(dataParticleList)::ParticleList 
   type(dataTimer)::Timer
   integer(8)::myid
   double precision,allocatable::seeds(:)
   integer::nthreads
   integer::iseed
end module dataTypes
