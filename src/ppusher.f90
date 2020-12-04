module ppusher
  !! @note Module with subroutines for advancing the particles' position and
  !! velocity in the simulations. @endnote
  use types
  use field
  use interp
  use hpc
  use coords
  use constants

  IMPLICIT NONE

  PUBLIC :: adv_eqn_top,&
       adv_interp_psi_top,&
       adv_interp_2DB_top,&
       adv_interp_3DB_top,&
#ifdef FIO
       adv_fio_top,&
       advance_fio_vars,&
#endif
       advance_eqn_vars,&
       advance_interp_psi_vars,&
       advance_interp_2DB_vars,&
       advance_interp_3DB_vars

contains

  subroutine adv_eqn_top(params,F,spp)

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
    REAL(rp)  :: Bo

    Bo=F%AB%Bo
    num_punct=params%num_punctures
    pchunk=params%pchunk

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z)
    do pp=1_idef,spp%ppp,pchunk

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
          
          call advance_eqn_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo
          
          do cc=1_idef,pchunk
             if ((Bo>0).and.(Y_PHI(cc)>Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                
                tt(cc)=tt(cc)+1   
             else if ((Bo<0).and.(Y_PHI(cc)<Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
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
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,num_punct
    !! Chunk iterator.
    INTEGER(ip),DIMENSION(params%pchunk)      :: tt
    !! time iterator.
    REAL(rp)  :: Bo

    Bo=F%AB%Bo
    num_punct=params%num_punctures
    pchunk=params%pchunk
    tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z)
    do pp=1_idef,spp%ppp,pchunk

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
          
          call advance_interp_psi_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and.(Y_PHI(cc)>Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1   
             else if ((Bo<0).and.(Y_PHI(cc)<Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
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

  end subroutine adv_interp_psi_top

  subroutine adv_interp_2DB_top(params,F,spp)

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
    REAL(rp)  :: Bo

    Bo=F%AB%Bo
    num_punct=params%num_punctures
    pchunk=params%pchunk
    tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z)
    do pp=1_idef,spp%ppp,pchunk

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
          
          call advance_interp_2DB_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and.(Y_PHI(cc)>Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1   
             else if ((Bo<0).and.(Y_PHI(cc)<Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
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
    INTEGER(ip),DIMENSION(params%pchunk)           :: tt
    !! time iterator.
    REAL(rp)  :: Bo

    Bo=F%AB%Bo
    num_punct=params%num_punctures
    pchunk=params%pchunk
    tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z)
    do pp=1_idef,spp%ppp,pchunk

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
          
          call advance_interp_3DB_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and.(Y_PHI(cc)>Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (-2*C_PI-Y_PHI0(cc))* &
                     (Y_Z(cc)-Y_Z0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))

                tt(cc)=tt(cc)+1   
             else if ((Bo<0).and.(Y_PHI(cc)<Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
                     (Y_R(cc)-Y_R0(cc))/(Y_PHI1(cc)-Y_PHI0(cc))
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+ &
                     (2*C_PI-Y_PHI0(cc))* &
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

  end subroutine adv_interp_3DB_top


#ifdef FIO
  subroutine adv_fio_top(params,F,spp)

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
    REAL(rp)  :: Bo
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint

    Bo=F%AB%Bo
    num_punct=params%num_punctures
    pchunk=params%pchunk
    tt=1

    !$OMP PARALLEL DO default(none) &
    !$OMP& FIRSTPRIVATE(Bo,pchunk,num_punct) &
    !$OMP& shared(F,params,spp) &
    !$OMP& PRIVATE(pp,tt,cc,Y_R,Y_PHI,Y_Z,Y_R0,Y_PHI0,Y_Z0,Y_PHI1, &
    !$OMP& flagCon,B_R,B_PHI,B_Z)
    do pp=1_idef,spp%ppp,pchunk

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
          
          call advance_fio_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon,hint)

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_PHI1(cc)=Y_PHI(cc)
             Y_PHI(cc)=modulo(Y_PHI(cc),2*C_PI)
          enddo

          do cc=1_idef,pchunk
             if ((Bo>0).and.(Y_PHI(cc)>Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+(-2*C_PI-Y_PHI0)* &
                     (Y_R-Y_R0)/(Y_PHI1-Y_PHI0)
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+(-2*C_PI-Y_PHI0)* &
                     (Y_Z-Y_Z0)/(Y_PHI1-Y_PHI0)

                tt(cc)=tt(cc)+1   
             else if ((Bo<0).and.(Y_PHI(cc)<Y_PHI0(cc))) then

                spp%vars%punct(tt(cc),pp-1+cc,1)=Y_R0(cc)+(2*C_PI-Y_PHI0)* &
                     (Y_R-Y_R0)/(Y_PHI1-Y_PHI0)
                spp%vars%punct(tt(cc),pp-1+cc,2)=Y_Z0(cc)+(2*C_PI-Y_PHI0)* &
                     (Y_Z-Y_Z0)/(Y_PHI1-Y_PHI0)

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
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
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
       k1_R(cc)=dx*B_R(cc)              
       k1_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k1_Z(cc)=dx*B_Z(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k2_R(cc)=dx*B_R(cc)              
       k2_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k2_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k3_R(cc)=dx*B_R(cc)              
       k3_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k3_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k4_R(cc)=dx*B_R(cc)              
       k4_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k4_Z(cc)=dx*B_Z(cc)

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
       k5_R(cc)=dx*B_R(cc)              
       k5_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k5_Z(cc)=dx*B_Z(cc)

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
       k6_R(cc)=dx*B_R(cc)              
       k6_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k6_Z(cc)=dx*B_Z(cc)

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

  subroutine advance_interp_psi_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)
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
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
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
       k1_R(cc)=dx*B_R(cc)              
       k1_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k1_Z(cc)=dx*B_Z(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k2_R(cc)=dx*B_R(cc)              
       k2_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k2_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)   

    !$OMP SIMD
    do cc=1_idef,pchunk
       k3_R(cc)=dx*B_R(cc)              
       k3_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k3_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call calculate_magnetic_field_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)   

    !$OMP SIMD
    do cc=1_idef,pchunk
       k4_R(cc)=dx*B_R(cc)              
       k4_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k4_Z(cc)=dx*B_Z(cc)

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
       k5_R(cc)=dx*B_R(cc)              
       k5_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k5_Z(cc)=dx*B_Z(cc)

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
       k6_R(cc)=dx*B_R(cc)              
       k6_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k6_Z(cc)=dx*B_Z(cc)

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

  end subroutine advance_interp_psi_vars

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
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
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
       k1_R(cc)=dx*B_R(cc)              
       k1_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k1_Z(cc)=dx*B_Z(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD
    
    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k2_R(cc)=dx*B_R(cc)              
       k2_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k2_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k3_R(cc)=dx*B_R(cc)              
       k3_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k3_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k4_R(cc)=dx*B_R(cc)              
       k4_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k4_Z(cc)=dx*B_Z(cc)

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
       k5_R(cc)=dx*B_R(cc)              
       k5_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k5_Z(cc)=dx*B_Z(cc)

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
       k6_R(cc)=dx*B_R(cc)              
       k6_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k6_Z(cc)=dx*B_Z(cc)

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
  
  subroutine advance_interp_3DB_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon)
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
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
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
       k1_R(cc)=dx*B_R(cc)              
       k1_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k1_Z(cc)=dx*B_Z(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD
    
    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k2_R(cc)=dx*B_R(cc)              
       k2_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k2_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k3_R(cc)=dx*B_R(cc)              
       k3_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k3_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call interp_2DB_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,flagCon)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k4_R(cc)=dx*B_R(cc)              
       k4_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k4_Z(cc)=dx*B_Z(cc)

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
       k5_R(cc)=dx*B_R(cc)              
       k5_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k5_Z(cc)=dx*B_Z(cc)

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
       k6_R(cc)=dx*B_R(cc)              
       k6_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k6_Z(cc)=dx*B_Z(cc)

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

  end subroutine advance_interp_3DB_vars

#ifdef FIO
  subroutine advance_fio_vars(params,F,Y_R,Y_PHI,Y_Z,flagCon,hint)
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
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon
    TYPE(C_PTR),DIMENSION(params%pchunk),INTENT(IN) :: hint

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
       B_X,B_Y,B_Z,flagCon,hint)  

    !$OMP SIMD
    do cc=1_idef,pchunk
       k1_R(cc)=dx*B_R(cc)              
       k1_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k1_Z(cc)=dx*B_Z(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_X,B_Y,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k2_R(cc)=dx*B_R(cc)              
       k2_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k2_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_X,B_Y,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k3_R(cc)=dx*B_R(cc)              
       k3_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k3_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_X,B_Y,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k4_R(cc)=dx*B_R(cc)              
       k4_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k4_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
         B_X,B_Y,B_Z,flagCon,hint)
        
    call calculate_magnetic_field
    !$OMP SIMD
    do cc=1_idef,pchunk
       k5_R(cc)=dx*B_R(cc)              
       k5_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k5_Z(cc)=dx*B_Z(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
    end do
    !$OMP END SIMD

    call get_fio_magnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_X,B_Y,B_Z,flagCon,hint)

    !$OMP SIMD
    do cc=1_idef,pchunk
       k6_R(cc)=dx*B_R(cc)              
       k6_PHI(cc)=dx*B_PHI(cc)/Y_R(cc)    
       k6_Z(cc)=dx*B_Z(cc)

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



end module ppusher
