submodule (dataTypes) paramKernel
   
contains
   
module subroutine initDefaultParams( &
   self,ni,nj,nk,             &
   mppc,mach,density,         &
   plate_x,plate_dy,plate_dz, &
   time                       &
   )

   implicit none 
   class(dataInParams)::self
   integer(8)::ni,nj,nk
   integer(8)::mppc
   double precision::mach
   double precision::density
   double precision::plate_x
   double precision::plate_dy
   double precision::plate_dz
   double precision::time

   
   if(self%ni<0)         self%ni           =ni
   if(self%nj<0)         self%nj           =nj
   if(self%nk<0)         self%nk           =nk
   if(self%mppc<0)       self%mppc         =mppc
   if(self%mach<0)       self%mach         =mach
   if(self%density<0)    self%density      =density
   if(self%plate_x<0)    self%plate_x      =plate_x
   if(self%plate_dy<0)   self%plate_dy     =plate_dy
   if(self%plate_dz<0)   self%plate_dz     =plate_dz
   if(self%time<0)       self%time         =time
   
end subroutine initDefaultParams

module subroutine processInParams(self)
   implicit none 
   class(dataInParams)::self
   
   double precision::deltax,tsteps
   integer(8)::ln2steps

   self%vmean=1.d0
   self%nCells=self%ni*self%nj*self%nk
   self%vtemp=self%vmean/self%mach

   self%dx=2.d0/dble(self%ni)   
   self%dy=2.d0/dble(self%nj)   
   self%dz=2.d0/dble(self%nk)
   self%cellVol=self%dx*self%dy*self%dz
   
   !! compute number of molecules a particle represents
   self%pnum=self%density*self%cellVol/dble(self%mppc) 
   print*,'pnum: ',self%pnum 
   
   !! Compute resonable timestep 
   deltax=2.d0/dble(max(max(self%ni,self%nj),self%nk))
   self%deltaT=0.10*deltax/(self%vmean+self%vtemp)

   !! If time duration  not given, simulate for 4 free-stream flow through times 
   if(self%time<0)self%time=8.d0/(self%vmean+self%vtemp)

   !! Compute nearest power of 2 timesteps 
   tsteps=self%time/self%deltaT
   ln2steps=int(ceiling(log(tsteps)/log(2.0)))
   self%nSteps=lshift(1,ln2steps)
end subroutine processInParams

module subroutine parseInParams(self)
    implicit none 
    class(dataInParams)::self
    integer(4)::nArgs,argID
    character(len=120)::arg

    nArgs=command_argument_count()

    argID=0
    MainLoop:do  
        argID=argID+1
        call get_command_argument(argID,Arg)
        select case(trim(adjustl(Arg)))
        case('-h','--help')
            if(argID>1)then
                print*,'-h/--help can only be  used as first argument'
                print*,'Try: dsmc -h/-help for details '
                stop
            else
                call Print_Help()
                stop
            end if 

        case('-2d')
           self%nk=1
           self%plate_dz=1

        case('-ni')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%ni

        case('-nj')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%nj

        case('-nk')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%nk

        case('-mppc')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%mppc

        case('-mach')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%mach

        case('-density')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%density

        case('-platex')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%plate_x

        case('-platedy')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%plate_dy

        case('-platedz')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%plate_dz

        case('-time')
           argID=argID+1
           call get_command_argument(argID,Arg)
           read(Arg,*)self%time
        case default
			write(*,'(A)')'Invalid argument('//trim(adjustl(arg))//')!'
            print*,'Try: dsmc -h/--help for details '
            stop
		end select 
        if(argID .eq. nArgs)exit MainLoop
    end do MainLoop
end subroutine parseInParams

subroutine Print_Help()

write(*,'(A)')" *****************************************************************************"
write(*,'(A)')" **"
write(*,'(A)')" ** Simple Direct Simulation Monte Carlo rarefied gas simulator"
write(*,'(A)')" ** Runs on cube with corners [1,1,1] and [-1,-1,-1], and a center at the"
write(*,'(A)')" ** origin.  x=-1 is an inflow and x=1 is an outflow, all other boundaries"
write(*,'(A)')" ** are periodic."
write(*,'(A)')" **"
write(*,'(A)')" ** Main program arguments:"
write(*,'(A)')" **"
write(*,'(A)')" ** -2d:"
write(*,'(A)')" **      Switch to two dimensional mode (only one cell in z directions)"
write(*,'(A)')" **"
write(*,'(A)')" ** -ni:"
write(*,'(A)')" **      Number of cells in x direction"
write(*,'(A)')" **"
write(*,'(A)')" ** -nj:"
write(*,'(A)')" **      Number of cells in y direction"
write(*,'(A)')" **"
write(*,'(A)')" ** -nk:"
write(*,'(A)')" **      Number of cells in z direction"
write(*,'(A)')" **"
write(*,'(A)')" ** -mppc:"
write(*,'(A)')" **      Mean Particles per Cell in simulation, adjust number of virtual"
write(*,'(A)')" **      particles to meet this target for the inflow"
write(*,'(A)')" **"
write(*,'(A)')" ** -mach:"
write(*,'(A)')" **      Mach number or ratio of mean atom velocity versus mean thermal velocity"
write(*,'(A)')" **"
write(*,'(A)')" ** -density:"
write(*,'(A)')" **      Density of the incoming flow"
write(*,'(A)')" **"
write(*,'(A)')" ** -platex"
write(*,'(A)')" **      x-location of plate"
write(*,'(A)')" **"
write(*,'(A)')" ** -platedy"
write(*,'(A)')" **      y height of plate"
write(*,'(A)')" **"
write(*,'(A)')" ** -platedz"
write(*,'(A)')" **      z width of plate"
write(*,'(A)')" **"
write(*,'(A)')" ** -time"
write(*,'(A)')" **      simulation time step size (usually computed from the above parameters"
write(*,'(A)')" **"
write(*,'(A)')" *****************************************************************************"
end subroutine Print_Help

end submodule paramKernel
