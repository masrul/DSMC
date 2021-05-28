submodule (dataTypes) collisionInfoKernel

contains

module subroutine initCollisionInfo(self)
   class(dataCollisionInfo)::self
   integer(8)::nCells,iCell
   double precision::rand


   do iCell=1,inputs%nCells
      self%maxCollisionRate(iCell)=inputs%sigmak*inputs%vtemp
      call random_number(rand)
      self%collisionRemainder(iCell)=rand
   enddo
end subroutine initCollisionInfo

end submodule collisionInfoKernel
