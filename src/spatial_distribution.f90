MODULE spatial_distribution
  !! @note Module with subroutines for generating the initial spatial distribution 
  !! of the different partciles' species in the simulation. @endnote
  USE types

  IMPLICIT NONE
  
  PUBLIC :: intitial_spatial_distribution

CONTAINS


subroutine intitial_spatial_distribution(params,spp,F)
  !! @note Subroutine that contains calls to the different subroutines 
  !! for initializing the simulated particles with various
  !! spatial distribution functions. @endnote
  TYPE(KORC_PARAMS), INTENT(INOUT) 			  :: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) :: spp
  !! An instance of the derived type SPECIES containing all the parameters and 
  !! simulation variables of the different species in the simulation.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  INTEGER :: ii,np
  REAL :: Rbuff
  
  SELECT CASE (TRIM(spp%spatial_distrib))
  CASE ('TRACER')
     spp%vars%Y(:,1)=spp%Xtrace(1)
     spp%vars%Y(:,2)=spp%Xtrace(2)
     spp%vars%Y(:,3)=spp%Xtrace(3)
  CASE ('HLINE')

     np=params%mpi_params%nmpi*spp%ppp
     Rbuff=(params%rmax-params%rmin)/20._rp
     
     do ii=1,np
        spp%vars%Y(ii,1)=(params%rmin+Rbuff)+ &
             (params%rmax-2._rp*Rbuff-params%rmin)*(ii-1)/(np-1)
     enddo
     spp%vars%Y(:,2)=0._rp
     spp%vars%Y(:,3)=0._rp
  CASE DEFAULT
     spp%vars%Y(:,1)=spp%Xtrace(1)
     spp%vars%Y(:,2)=spp%Xtrace(2)
     spp%vars%Y(:,3)=spp%Xtrace(3)
  END SELECT

    
end subroutine intitial_spatial_distribution


END MODULE spatial_distribution
