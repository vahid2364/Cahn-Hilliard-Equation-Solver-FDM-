!! filename = CH.f90
!! Vahid Attari
!! Created: 30 Feb. 2016
!! Modified: ....
!! Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
!!
!! Acknowledgements:  Based on Cahn-Hilliard 1965 paper
!!
!! Purpose:
!!   - Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases
!!     to self consistantly model the Spinodal Composition Phenomenon
!!   
!! General Algorithm function:
!!
!!   1. Retrieve parameter data from file "parameters.dat"
!!   2. Assess thermodynamics of the associated system 
!!   3. Reads initial phase distribution from "phase.dat" file
!!   4. Calculate Phase Evolution with time integration
!!      -  Nucleate Phases
!!      -  Resolve boundary conditions 
!!         -- Periodic boundaries in all directions
!!      -  Solve differential equations via 9-stencil finite difference
!!      -  Update phase information and concentration data
!!
!! Compilation instructions: >> make
!!    - Manual: >>  ifort -o a.out CH.f90
!!
!! Execution: >> ./a.out 
!!                                     
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!====================================================================================

!   This code simulates the early stage of Spinodal Decomposition...



!   Cahn- Hilliard solver 
!   with periodic boundary conditions.
!   Nonlinear term: f(u) = u - u**3

module INITIALIZATION

    integer, parameter::nx=300, ny=300
    real, dimension(nx,ny) :: phi,lap_phi,f,df,d2f
    real, dimension(nx,ny) :: lap_chem_pot,chem_pot
    real,save :: phiold, dx, dy, dt,h,pi!,!mobility,kappa
    real :: phitot
    real :: Betta_C,Betta_Max,Landa_C,Landa_Max,R_Max,betta,R
    integer, save :: l1, m1


end module INITIALIZATION

!==============================================================

module mts_pars

	! Energy parameters
    real,parameter :: a2 = -1
    real,parameter :: a4 = 1;	
	
	! Kinetic parameters
    real,parameter :: mobility = 1.0D0
    real,parameter :: kappa = 0.03D0

end module mts_pars

!==============================================================

MODULE wrt_opts

    CHARACTER ch*9
    CHARACTER(LEN=100) :: VAL
    CHARACTER(len=255) :: cwd,fileplace

    INTEGER, PARAMETER :: wrt_cycle  = 5000            !100000
    INTEGER, PARAMETER :: file_cycle = 20*wrt_cycle    !20000
    INTEGER, PARAMETER :: stop_iter  = 1*file_cycle   !100000
    INTEGER :: NNN3      = 1000

    integer            :: IPC
    integer            :: itimes
    real               :: Time

CONTAINS

    SUBROUTINE mk_dir

        !!***** MAKE DIRs *****
        call system('rm -r parameters')
		call system('rm -r microstructure')
        call system('rm -r results')
        !call system('rm -r heat')
        call system('rm *.mod')
        call system('rm *.dat')
        call system('rm *.plt')

        call system('mkdir results')
        call system('mkdir microstructure')
        !call system('mkdir parameters')
        !call system('mkdir heat')

    END SUBROUTINE

    SUBROUTINE wrt_disk

        !!***** MAKE DIRs *****
        CALL getcwd(cwd)
        WRITE(fileplace,*) ADJUSTL(TRIM(cwd))

    END SUBROUTINE

END MODULE wrt_opts

!==============================================================

program main
    use INITIALIZATION
    use wrt_opts
    use mts_pars
    implicit none
    integer::i,j
			
			
    IPC    = 0
    itimes = 0
	Time   = 0		    
    !! ************* INITIALIZATION
	call mk_dir		
    call init

    !! *********************************************
	OPEN(unit=3,file='results/itimes.dat')
	WRITE(3,*) IPC,itimes,phitot,Betta_C,Betta_Max
	!! *********************************************
	    
    !! ************* CALCULATING
    do itimes=1,stop_iter

        call gradandlaplace
        call evolution

        IF ( itimes.EQ. 1 .or. MOD(itimes,wrt_cycle).EQ.0 ) THEN
            IPC = IPC + 1
            CALL analyze
            CALL printdata(itimes)
            !! *********************************************
			WRITE(3,*) IPC,itimes,phitot,Betta_C,Betta_Max
			!! *********************************************
        ENDIF
        
        Time = Time + dt
        
    end do

    print *, 'SIMULATION FINISHED... GO HOME!!!'

    stop

end program main

!=================================================================

subroutine init
    use INITIALIZATION
    implicit none
    integer::i,j
    real::ranum


    h= 0.3D0

    dx=h
    dy=h
    l1=nx
    m1=ny
        
    dt=(dx)**4/32    !0.0001

    phitot = 0.0D0
    
    !setting the domain
    do i=1,nx
        do j=1,ny

            call random_number(ranum)
            phi(i,j)= 0.0d0 + 0.02*ranum
            phitot = phitot + phi(i,j)
 
        end do
    end do
    
    
    !OPEN(unit=4,file='betta.dat')
     
!    OPEN(unit=1,file='INITIAL.dat')
!    WRITE(1,*) 'ZONE ','I=',nx,'J= ',ny

!    ! INITITAL.dat file
!    do i=1,nx ; do j=1,ny
!        WRITE(1,*) i,j,phi(i,j)
!    end do ; end do
    
end subroutine init

!=======================================================================================================================

subroutine freeenergy(ii,jj,FE,dfdc,d2fdc2)
    
    use INITIALIZATION
    use mts_pars

    implicit none
    integer :: ii,jj
    real, dimension(nx,ny) :: FE,dfdc,d2fdc2
    

    FE(ii,jj)   = (a2/2)*phi(ii,jj)**2 + (a4/4)*phi(ii,jj)**4
    dfdc(ii,jj)  = a4 * phi(ii,jj)**3 + a2 * phi(ii,jj)
    d2fdc2(ii,jj) = 3*phi(ii,jj)**2 - 1
    
end subroutine

!=======================================================================================================================

subroutine analyze
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,nnn 
        
    PI=4.D0*DATAN(1.D0)
    
    !---------------------------
    !! System Total Mass : For Controlling the mass conservation
    phitot = phitot/DFLOAT((nx)*(ny))
    !---------------------------
        
    !---------------------------
    do i=1,nx; do j=1,ny;
        
        call freeenergy(i,j,f,df,d2f)
        ! wavenumber
        Betta_Max  = 0.5*sqrt(-d2f(i,j)/kappa)
        Betta_C    = Betta_Max*sqrt(2.0D0)

        ! wavelength
        Landa_Max  = 2*PI/Betta_Max                 ! The fastest growing wavelength
        Landa_C    = 2*PI/Betta_C 
        
        ! amplification factor
        R_Max      = mobility*((d2f(i,j))**2/8*kappa)

    end do; end do;
    !---------------------------
    
    !---------------------------
    betta = 0
!    WRITE(4,*) 'ZONE'
    do i=0,100,1
        R          = -mobility*(d2f(10,10))*betta**2-2*mobility*kappa*betta**4
        betta      = betta + 0.02 
        !WRITE(4,*) betta_C/betta,R/R_Max
    enddo
        
            
end subroutine analyze   
    
!================================================================================

subroutine gradandlaplace
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,ip,im,jp,jm
    integer :: ipp,imm,jpp,jmm
    real :: result

    do i=1,nx
        do j=1,ny


            !finite difference step
            ip = i+1  ; im = i-1  ; jp = j+1  ; jm = j-1  ; 
            ipp = i+2 ; imm = i-2 ; jpp = j+2 ; jmm = j-2 ;
            
            !periodic boundary condition
            if (im==0) im=l1
            if (jm==0) jm=m1
            if (ip==l1+1) ip=1
            if (jp==m1+1) jp=1
            
            if (imm==-1) imm=l1-1
            if (jmm==-1) jmm=m1-1
            if (imm==0)  imm=l1
            if (jmm==0)  jmm=m1
            
            if (ipp==l1+1) ipp=1
            if (jpp==m1+1) jpp=1
            if (ipp==l1+2) ipp=2
            if (jpp==m1+2) jpp=2

            !---------------------------

            !laplacians (9-stencil)
            lap_phi(i,j) = (2.0*(phi(ip,j)+phi(im,j)+phi(i,jp)+ phi(i,jm)) + phi(ip,jp)+phi(im,jm)+phi(im,jp)+ &
            phi(ip,jm) - 12.0*phi(i,j))/(3.0*dx*dx)
            
            !laplacians (13-stencil)
!            lap_phi(i,j) = (20*phi(i,j) - 8*(phi(ip,j)+phi(im,j)+phi(i,jp)+phi(i,jm)) &
!            + 2* (phi(im,jm) + phi(im,jp) + phi(ip,jm) + phi(ip,jp)) &
!            + phi(imm,j) + phi(ipp,j) + phi(i,jmm) + phi(i,jpp))/(dx*dx*dx*dx)
            
            result = lap_phi(i,j)
            IF(ISNAN(result)) THEN
            WRITE(*,*) i,j
            WRITE(*,*) lap_phi(i,j)
            WRITE(*,*) phi(i,j),phi(ip,j),phi(im,j),phi(i,jp),phi(i,jm), &
            phi(im,jm), phi(im,jp), phi(ip,jm) , phi(ip,jp), &
            phi(imm,j) ,phi(ipp,j) , phi(i,jmm) , phi(i,jpp),dx
            WRITE(*,*) (20*phi(i,j) - 8*(phi(ip,j)+phi(im,j)+phi(i,jp)+phi(i,jm)) &
            + 2*(phi(im,jm) + phi(im,jp) + phi(ip,jm) + phi(ip,jp)) &
            + phi(imm,j) + phi(ipp,j) + phi(i,jmm) + phi(i,jpp))
            WRITE(*,*) (1/dx**4)
            pause
            ENDIF 
            
            call freeenergy(i,j,f,df,d2f)
                                    
            chem_pot(i,j) = df(i,j) - 2*kappa*lap_phi(i,j)
                        
        end do
    end do

end subroutine gradandlaplace
!===============================================================
subroutine evolution
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,ip,im,jp,jm
    integer :: ipp,imm,jpp,jmm
    phitot = 0.0D0
     
    do i=1,nx
        do j=1,ny
        
            !finite difference step
            ip = i+1  ; im = i-1  ; jp = j+1  ; jm = j-1  ; 
            ipp = i+2 ; imm = i-2 ; jpp = j+2 ; jmm = j-2 ;
            
            !periodic boundary condition
            if (im==0) im=l1
            if (jm==0) jm=m1
            if (ip==l1+1) ip=1
            if (jp==m1+1) jp=1
            
            if (imm==-1) imm=l1-1
            if (jmm==-1) jmm=m1-1
            if (imm==0)  imm=l1
            if (jmm==0)  jmm=m1
            
            if (ipp==l1+1) ipp=1
            if (jpp==m1+1) jpp=1
            if (ipp==l1+2) ipp=2
            if (jpp==m1+2) jpp=2
            
            ! Laplace of chemical potential (9-stencil) 
            lap_chem_pot(i,j) = (2.0*(chem_pot(ip,j)+chem_pot(im,j)+chem_pot(i,jp)+ &
            chem_pot(i,jm)) + chem_pot(ip,jp)+chem_pot(im,jm)+chem_pot(im,jp)+ &
            chem_pot(ip,jm) - 12.0*chem_pot(i,j))/(3.0*dx*dx)
            
            ! Laplace of chemical potential (13-stencil)
!            lap_chem_pot(i,j) = (20*chem_pot(i,j) - 8*(chem_pot(ip,j)+chem_pot(im,j)+chem_pot(i,jp)+chem_pot(i,jm)) &
!            + 2* (chem_pot(im,jm) + chem_pot(im,jp) + chem_pot(ip,jm) + chem_pot(ip,jp)) &
!            + chem_pot(imm,j) + chem_pot(ipp,j) + chem_pot(i,jmm) + chem_pot(i,jpp))/(dx*dx*dx*dx)

            phi(i,j) = phi(i,j) + mobility*lap_chem_pot(i,j)*dt
      
            phitot = phitot + phi(i,j)

        end do
    end do

end subroutine evolution
!==================================================================

subroutine printdata(k)
    use INITIALIZATION
	use wrt_opts
	use mts_pars
    IMPLICIT NONE
    INTEGER :: k,i,j
    CHARACTER(len=80) :: FileName2
    
    !! *********** WRITTING OUTPUT *************** 
    WRITE(FileName2,FMT='(A30,I6.6,A4)') "microstructure/phi_",k,".dat"
    FileName2 = TRIM(FileName2)
    
    OPEN(unit=52,file=FileName2)
    WRITE(52,*) 'ZONE ','I=',nx,'J= ',ny
  
    DO j=1,ny
        DO i=1,nx
            WRITE(52,771) i,j,phi(i,j)
        ENDDO
    END DO

    CLOSE(52)
771 Format(I3,1X,i3,1X,1(1X,F8.4))
    !! *********************************************
    
    
    !! ************* DISPLAY OUTPUT ****************
    WRITE(*,*) '******Calculating...******'
    WRITE(*,*) 'IPC=',IPC,'Time=',Time,'itimes=',itimes
    WRITE(*,*) 'PHITOT=',phitot
    WRITE(*,*) 'Max. W.#.=',Betta_Max,'C. W.#.=',Betta_C
    WRITE(*,*) 'Max. W.L.=',Landa_Max,'C. W.L.=',Landa_C
    WRITE(*,*) 'R_Max=',R_Max
    WRITE(*,*) 'kappa=',kappa,'M=',Mobility
    WRITE(*,*)   
    !! *********************************************
    
    RETURN

end subroutine printdata  
