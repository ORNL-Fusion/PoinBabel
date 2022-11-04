module PB_ppusher
  !! @note Module with subroutines for advancing the particles' position and
  !! velocity in the simulations. @endnote
  use PB_types
  use PB_fields
  use PB_interp
  use PB_hpc
  use PB_coords
  use PB_constants

  IMPLICIT NONE

contains

  subroutine adv_eqn_top(params,F,spp)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)    :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section
    INTEGER             :: thread_num

    Bo=F%Bo
    phi_section=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z,thread_num)
    do pp=1_idef,spp%ppp,pchunk

       thread_num = OMP_GET_THREAD_NUM()

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)

          end do

          !write(6,*) thread_num,'R0',Y_R0
          !write(6,*) thread_num,'PHI0',Y_PHI0
          !write(6,*) thread_num,'Z0',Y_Z0

          call advance_eqn_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          if ((pchunk.eq.1).and.(flagCon(1).eq.0)) then
             write(6,*) thread_num,'trace left domain'
             exit
          endif

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)> &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                !write(6,*) thread_num,'R1',Y_R
                !write(6,*) thread_num,'PHI1',Y_PHI1
                !write(6,*) thread_num,'Z1',Y_Z

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                tt(cc)=tt(cc)+1
             else if ((Bo<0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)< &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                if (phi_section==0._rp) phi_section=phi_section+2*C_PI

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1
             endif
          enddo

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_eqn_top

#ifdef PSPLINE
    subroutine adv_interp_psi_top(params,F,spp)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: con_len
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)      :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section

    Bo=F%Bo
    phi_section=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z,con_len)
    do pp=1_idef,spp%ppp,pchunk

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)

          con_len(cc)=spp%vars%con_len(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)
          end do

          !write(6,*) 'R0',Y_R0
          !write(6,*) 'PHI0',Y_PHI0
          !write(6,*) 'Z0',Y_Z0

          if(params%output_orbit) then
             write(orbit_unit_write,*) '(R,PHI,Z-0):',Y_R0,Y_PHI0,Y_Z0
          end if

          call advance_interp_psi_vars(params,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
               flagCon,con_len)

          if(params%output_orbit) then
             write(orbit_unit_write,*) '(BR,BPHI,BZ)',B_R,B_PHI,B_Z
          end if

          if ((pchunk.eq.1).and.(flagCon(1).eq.0)) then
             write(6,*) 'trace left domain'
             exit
          endif

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)> &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                !write(6,*) 'R1',Y_R
                !write(6,*) 'PHI1',Y_PHI1
                !write(6,*) 'Z1',Y_Z

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                tt(cc)=tt(cc)+1
             else if ((Bo<0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)< &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                if (phi_section==0._rp) phi_section=phi_section+2*C_PI

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1
             endif
          enddo

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
          spp%vars%con_len(pp-1+cc)=con_len(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_interp_psi_top

  subroutine adv_interp_mars_top(params,F,spp)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)      :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section
    INTEGER             :: thread_num

    Bo=F%Bo
    phi_section=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk
        tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z,thread_num)
    do pp=1_idef,spp%ppp,pchunk

       thread_num = OMP_GET_THREAD_NUM()

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)
          end do
          !$OMP END SIMD

          !write(6,*) thread_num,'R0',Y_R0
          !write(6,*) thread_num,'PHI0',Y_PHI0
          !write(6,*) thread_num,'Z0',Y_Z0

          call advance_interp_mars_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          if ((pchunk.eq.1).and.(flagCon(1).eq.0)) then
             write(6,*) thread_num,'trace left domain'
             exit
          endif

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo
          !$OMP END SIMD

          do cc=1_idef,pchunk
             if ((Bo>0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)> &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                !write(6,*) thread_num,'R1',Y_R
                !write(6,*) thread_num,'PHI1',Y_PHI1
                !write(6,*) thread_num,'Z1',Y_Z

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                tt(cc)=tt(cc)+1
             else if ((Bo<0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)< &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                if (phi_section==0._rp) phi_section=phi_section+2*C_PI

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1
             endif
          enddo

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_interp_mars_top

  subroutine adv_interp_2DB_top(params,F,spp)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk) :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section
    INTEGER             :: thread_num

    Bo=F%Bo
    phi_section=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z,thread_num)
    do pp=1_idef,spp%ppp,pchunk

       thread_num = OMP_GET_THREAD_NUM()

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)
          end do

          !write(6,*) thread_num,'R0',Y_R0
          !write(6,*) thread_num,'PHI0',Y_PHI0
          !write(6,*) thread_num,'Z0',Y_Z0

          call advance_interp_2DB_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          do cc=1_idef,pchunk
             if ((pchunk.eq.1).and.(flagCon(cc).eq.0)) then
                write(6,*) thread_num,'trace left domain'
                exit
             endif
          end do

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)> &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                !write(6,*) thread_num,'R1',Y_R
                !write(6,*) thread_num,'PHI1',Y_PHI1
                !write(6,*) thread_num,'Z1',Y_Z

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                tt(cc)=tt(cc)+1
             else if ((Bo<0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)< &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                if (phi_section==0._rp) phi_section=phi_section+2*C_PI

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1
             endif
          enddo

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_interp_2DB_top

  subroutine adv_interp_3DB_top(params,F,spp)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,B_PHI0
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)           :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section0,phi_section1
    INTEGER             :: thread_num

    Bo=F%Bo
    phi_section0=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section0,phi_section1) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_PHI0,B_Z,thread_num)
    do pp=1_idef,spp%ppp,pchunk

       thread_num = OMP_GET_THREAD_NUM()

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)
          end do

          if(params%output_orbit) then
             write(orbit_unit_write,*) thread_num,pp,'(R,PHI,Z-0):',Y_R0,Y_PHI0,Y_Z0
          end if

          call advance_interp_3DB_vars(params,F,Y_R,Y_PHI,Y_Z,&
               B_R,B_PHI,B_Z,flagCon)

          if(params%output_orbit) then
             !write(orbit_unit_write,*) thread_num,pp,'(R,PHI,Z-1):',Y_R,Y_PHI,Y_Z
             write(orbit_unit_write,*) thread_num,pp,'(BR,BPHI,BZ)',B_R,B_PHI,B_Z
          end if

          if ((pchunk.eq.1).and.(flagCon(1).eq.0)) then
             write(6,*) thread_num,'trace left domain'
             exit
          endif

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          if (F%stel_sym) then
             do cc=1_idef,pchunk
                if ((B_PHI(cc)<0.and.(.not.(B_PHI(cc)*B_PHI0(cc).lt.0))).and. &
                     (modulo(Y_PHI(cc)-phi_section0,2*C_PI)> &
                     modulo(Y_PHI0(cc)-phi_section0,2*C_PI))) then

                   !write(6,*) 'tracing clockwise'
                   !write(6,*) 'phi_section',phi_section0
                   !write(6,*) thread_num,'R0',Y_R0
                   !write(6,*) thread_num,'PHI0',Y_PHI0
                   !write(6,*) thread_num,'Z0',Y_Z0

                   !write(6,*) thread_num,'R1',Y_R
                   !write(6,*) thread_num,'PHI1',Y_PHI1
                   !write(6,*) thread_num,'Z1',Y_Z

                   spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                   spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                   !write(6,*) 'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                   !write(6,*) 'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                   tt(cc)=tt(cc)+1
                else if ((B_PHI(cc)>0.and.(.not.(B_PHI(cc)*B_PHI0(cc).lt.0))).and. &
                     (modulo(Y_PHI(cc)-phi_section0,2*C_PI)< &
                     modulo(Y_PHI0(cc)-phi_section0,2*C_PI))) then

                   if (phi_section0==0._rp) THEN
                      phi_section1=phi_section0+2*C_PI
                   else
                      phi_section1=phi_section0
                   ENDIF

                   !write(6,*) 'tracing counter-clockwise'
                   !write(6,*) 'phi_section0',phi_section0
                   !write(6,*) 'phi_section1',phi_section1
                   !write(6,*) thread_num,'R0',Y_R0
                   !write(6,*) thread_num,'PHI0',Y_PHI0
                   !write(6,*) thread_num,'Z0',Y_Z0

                   !write(6,*) thread_num,'R1',Y_R
                   !write(6,*) thread_num,'PHI1',Y_PHI1
                   !write(6,*) thread_num,'Z1',Y_Z

                   spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                        (phi_section1-Y_PHI0(cc))* &
                        (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                   spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                        (phi_section1-Y_PHI0(cc))* &
                        (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                   !write(6,*) 'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                   !write(6,*) 'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                   tt(cc)=tt(cc)+1
                endif
             enddo
          else
             do cc=1_idef,pchunk
                if ((Bo>0).and. &
                     (modulo(Y_PHI(cc)-phi_section0,2*C_PI)> &
                     modulo(Y_PHI0(cc)-phi_section0,2*C_PI))) then

                   !write(6,*) thread_num,'R1',Y_R
                   !write(6,*) thread_num,'PHI1',Y_PHI1
                   !write(6,*) thread_num,'Z1',Y_Z

                   spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                   spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                   !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                   !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                   tt(cc)=tt(cc)+1
                else if ((Bo<0).and. &
                     (modulo(Y_PHI(cc)-phi_section0,2*C_PI)< &
                     modulo(Y_PHI0(cc)-phi_section0,2*C_PI))) then

                   if (phi_section0==0._rp) phi_section0=phi_section0+2*C_PI

                   spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                   spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                        (phi_section0-Y_PHI0(cc))* &
                        (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                   tt(cc)=tt(cc)+1
                endif
             enddo
          endif

          B_PHI0=B_PHI

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_interp_3DB_top
#endif

#ifdef FIO
  subroutine adv_fio_top(params,F,spp)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R0,Y_PHI0,Y_Z0,Y_PHI1
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)        :: tt
    !! time iterator.
    REAL(rp)  :: Bo,phi_section
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint
    INTEGER             :: thread_num

    Bo=F%Bo
    phi_section=params%phi_section
    num_punct=params%num_punctures
    pchunk=params%pchunk
    tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct,phi_section) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z,thread_num,hint)
    do pp=1_idef,spp%ppp,pchunk

       thread_num = OMP_GET_THREAD_NUM()

       !$OMP SIMD
       do cc=1_idef,pchunk
          tt(cc)=1

          Y_R(cc)=spp%vars%Y(pp-1+cc,1)
          Y_PHI(cc)=spp%vars%Y(pp-1+cc,2)
          Y_Z(cc)=spp%vars%Y(pp-1+cc,3)

          hint(cc)=spp%vars%hint(pp-1+cc)

          flagCon(cc)=spp%vars%flagCon(pp-1+cc)
       end do
       !$OMP END SIMD

       do while (maxval(tt).le.num_punct)

          !do cc=1_idef,pchunk
          !   if (mod(tt(cc),num_punct/10).eq.0) then
          !      write(6,*) thread_num,tt,cc
          !   end if
          !end do

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R0(cc)=Y_R(cc)
             Y_PHI0(cc)=Y_PHI(cc)
             Y_Z0(cc)=Y_Z(cc)
          end do
          !$OMP END SIMD

          !write(6,*) thread_num,'R0',Y_R0
          !write(6,*) thread_num,'PHI0',Y_PHI0
          !write(6,*) thread_num,'Z0',Y_Z0

          !write(6,*) phi_section

          call advance_fio_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon,hint)

          if ((pchunk.eq.1).and.(flagCon(1).eq.0)) then
             write(6,*) thread_num,'trace left domain'
             exit
          endif

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo
          !$OMP END SIMD


          !write(6,*) Bo
          !write(6,*) modulo(Y_PHI0(1)-phi_section,2*C_PI)
          !write(6,*) modulo(Y_PHI0(1)-.1-phi_section,2*C_PI)

          do cc=1_idef,pchunk
             if ((Bo>0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)> &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                !write(6,*) thread_num,'R1',Y_R
                !write(6,*) thread_num,'PHI1',Y_PHI1
                !write(6,*) thread_num,'Z1',Y_Z

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                !write(6,*) thread_num,'R_punct',spp%vars%punct(tt(cc),pp-1+cc,1)
                !write(6,*) thread_num,'Z_punct',spp%vars%punct(tt(cc),pp-1+cc,2)

                tt(cc)=tt(cc)+1
             else if ((Bo<0).and. &
                  (modulo(Y_PHI(cc)-phi_section,2*C_PI)< &
                  modulo(Y_PHI0(cc)-phi_section,2*C_PI))) then

                if (phi_section==0._rp) phi_section=phi_section+2*C_PI

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (phi_section-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1
             endif
          enddo

       end do !timestep iterator

       !$OMP SIMD
       do cc=1_idef,pchunk
          spp%vars%flagCon(pp-1+cc)=flagCon(cc)
       end do
       !$OMP END SIMD


    end do !particle chunk iterator
    !$OMP END PARALLEL DO

  end subroutine adv_fio_top
#endif

  subroutine advance_eqn_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon

    REAL(rp) :: ar,R0

    ar=F%AB%a
    R0=F%AB%Ro

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    call cyl_check_if_confined_p(pchunk,ar,R0,Y_R,Y_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       end if

    end do
    !$OMP END SIMD

  end subroutine advance_eqn_vars

#ifdef PSPLINE
  subroutine advance_interp_psi_vars(params,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon,con_len)
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z,con_len
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       else
          con_len(cc)=con_len(cc)+dx
       end if

    end do
    !$OMP END SIMD

  end subroutine advance_interp_psi_vars

  subroutine advance_interp_mars_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    call interp_3Dmars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
         flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    !write(6,*) 'bmag',Bmag
    !write(6,*) 'RHS',B_R/Bmag,B_PHI/(Y_R*Bmag),B_Z/Bmag

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       end if

    end do
    !$OMP END SIMD

  end subroutine advance_interp_mars_vars

    subroutine advance_interp_2DB_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       end if

    end do
    !$OMP END SIMD

  end subroutine advance_interp_2DB_vars

  subroutine advance_interp_3DB_vars(params,F,Y_R,Y_PHI,Y_Z,&
       B_R,B_PHI,B_Z,flagCon)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_Z,B_PHI
    REAL(rp),DIMENSION(params%pchunk) :: Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon
    INTEGER             :: thread_num

    thread_num = OMP_GET_THREAD_NUM()

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !write(6,*) thread_num,'BR0',B_R
    !write(6,*) thread_num,'BPHI0',B_PHI
    !write(6,*) thread_num,'BZ0',B_Z

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       end if

    end do
    !$OMP END SIMD

    call interp_3DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

  end subroutine advance_interp_3DB_vars
#endif

#ifdef FIO
  subroutine advance_fio_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon,hint)

    use omp_lib

    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dx
    !! Time step used in the leapfrog step (\(\Delta t\)).
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,Bmag
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon
    TYPE(C_PTR),DIMENSION(params%pchunk),INTENT(INOUT) :: hint
    INTEGER             :: thread_num

    thread_num = OMP_GET_THREAD_NUM()

    pchunk=params%pchunk
    dx=params%dx

    !$OMP SIMD
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon,hint)

    !write(6,*) thread_num,'BR0',B_R
    !write(6,*) thread_num,'BPHI0',B_PHI
    !write(6,*) thread_num,'BZ0',B_Z

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k1_R(cc)=dx*B_R(cc)/Bmag(cc)
       k1_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k1_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k2_R(cc)=dx*B_R(cc)/Bmag(cc)
       k2_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k2_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k3_R(cc)=dx*B_R(cc)/Bmag(cc)
       k3_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k3_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k4_R(cc)=dx*B_R(cc)/Bmag(cc)
       k4_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k4_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
         B_R,B_PHI,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k5_R(cc)=dx*B_R(cc)/Bmag(cc)
       k5_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k5_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       k6_R(cc)=dx*B_R(cc)/Bmag(cc)
       k6_PHI(cc)=dx*B_PHI(cc)/(Y_R(cc)*Bmag(cc))
       k6_Z(cc)=dx*B_Z(cc)/Bmag(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    do cc=1_idef,pchunk

       if (flagCon(cc).eq.0_is) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
       end if

    end do
    !$OMP END SIMD

  end subroutine advance_fio_vars
#endif



end module PB_ppusher
