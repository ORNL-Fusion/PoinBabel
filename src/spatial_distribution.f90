MODULE spatial_distribution
  !! @note Module with subroutines for generating the initial spatial distribution 
  !! of the different partciles' species in the simulation. @endnote
  USE types
  USE PB_input

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
  INTEGER :: ii,np,num_part_in,ios
  REAL(rp) :: Rbuff,phi_start,R_tmp,Z_tmp,rmin,rmax,Ro
  INTEGER, PARAMETER :: pos_file = 102

  rmax=params%rmax
  rmin=params%rmin
  Ro=rmin+(rmax-rmin)/2
  Rbuff=(rmax-rmin)/20._rp
  np=params%mpi_params%nmpi*spp%ppp
  
  SELECT CASE (TRIM(spp%spatial_distrib))
  CASE ('TRACER')
     spp%vars%Y(:,1)=spp%Xtrace(1)
     spp%vars%Y(:,2)=spp%Xtrace(2)
     spp%vars%Y(:,3)=spp%Xtrace(3)
  CASE ('HLINE')
     if (np.eq.1) call PB_abort(20)

     do ii=1,np
        spp%vars%Y(ii,1)=(rmin+Rbuff)+ &
             (rmax-2._rp*Rbuff-rmin)*(ii-1)/(np-1)        
     enddo
                       
  CASE('RADII')
     if (np.eq.1) call PB_abort(20)
     
     do ii=1,np
        spp%vars%Y(ii,1)=(Ro+Rbuff)+ &
             (rmax-2._rp*Rbuff-Ro)*(ii-1)/(np-1)      
     enddo    
  CASE ('INPUT')

     OPEN(UNIT=pos_file,FILE=position_filename,STATUS='OLD')

     READ(UNIT=pos_file,FMT='(I4,F3.2)',IOSTAT=ios) num_part_in,phi_start
     if (ios/=0) call PB_abort(20)

!     write(6,'(I4,F3.2)') num_part_in,phi_start
     
     if (num_part_in.gt.np) call PB_abort(20)
     
     do ii=1,num_part_in
        READ(UNIT=pos_file,FMT='(F8.5,F8.5)',IOSTAT=ios) R_tmp,Z_tmp


!        write(6,'(F7.5,F7.5)') R_tmp,Z_tmp
        
        spp%vars%Y(ii,1)=R_tmp
        spp%vars%Y(ii,2)=phi_start
        spp%vars%Y(ii,3)=Z_tmp

!        write(6,'(I1)') ii
!        write(6,*) spp%vars%Y(ii,1)
!        write(6,'(F7.2)') spp%vars%Y(ii,2)
!        write(6,'(F8.5)') spp%vars%Y(ii,3)
     end do


     
     if (num_part_in.lt.np) then
        do ii=num_part_in+1,np
               spp%vars%flagCon(ii) = 0_is
        enddo        
     endif
     
  CASE DEFAULT
     spp%vars%Y(:,1)=spp%Xtrace(1)
     spp%vars%Y(:,2)=spp%Xtrace(2)
     spp%vars%Y(:,3)=spp%Xtrace(3)
  END SELECT

    
end subroutine intitial_spatial_distribution


END MODULE spatial_distribution
