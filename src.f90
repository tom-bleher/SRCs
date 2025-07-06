! src.f90
!
! Translated and optimized from the original C++ version.
!
MODULE PhysicalConstants
    USE, INTRINSIC :: iso_fortran_env, ONLY: REAL64
    IMPLICIT NONE
    REAL(KIND=REAL64), PARAMETER :: PI = 3.14159265358979323846
    REAL(KIND=REAL64), PARAMETER :: HBARC = 197.3269804  ! MeV·fm
    REAL(KIND=REAL64), PARAMETER :: M_N = 939.0           ! MeV
END MODULE PhysicalConstants

MODULE GlobalData
    USE, INTRINSIC :: iso_fortran_env, ONLY: REAL64, INT32
    IMPLICIT NONE
    
    ! Shell state definition (n_r, l)
    TYPE :: ShellState
        INTEGER(KIND=INT32) :: n_r
        INTEGER(KIND=INT32) :: l
    END TYPE ShellState
    
    ! Nucleon definition
    TYPE :: Nucleon
        INTEGER(KIND=INT32) :: type          ! 0 = proton, 1 = neutron
        INTEGER(KIND=INT32) :: shell_index   ! index into protons or neutrons vector
        INTEGER(KIND=INT32) :: center_index  ! which spatial center it's at
        REAL(KIND=REAL64)   :: probability   ! probability of being at this center
        REAL(KIND=REAL64)   :: transparency  ! transparency of nucleon
    END TYPE Nucleon

    ! Global arrays for results
    REAL(KIND=REAL64), ALLOCATABLE, SAVE :: np(:,:,:,:), pp(:,:,:,:), nn(:,:,:,:)
    
    ! Shell configurations
    TYPE(ShellState), ALLOCATABLE, SAVE :: protons(:), neutrons(:)

END MODULE GlobalData

MODULE Procedures
    USE GlobalData
    USE PhysicalConstants
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: R_radial, GenDensity, getTransparency, configure_nucleus

    CONTAINS

    ! Efficient three-term recurrence for Laguerre polynomials
    FUNCTION assoc_laguerre(n, alpha, x) RESULT(res)
        INTEGER(KIND=INT32), INTENT(IN) :: n
        REAL(KIND=REAL64), INTENT(IN) :: alpha, x
        REAL(KIND=REAL64) :: res
        REAL(KIND=REAL64) :: L0, L1, Ln
        INTEGER(KIND=INT32) :: i
        
        IF (n < 0) THEN
            res = 0.0_REAL64
            RETURN
        ELSE IF (n == 0) THEN
            res = 1.0_REAL64
            RETURN
        ELSE IF (n == 1) THEN
            res = 1.0_REAL64 + alpha - x
            RETURN
        END IF
        
        ! Three-term recurrence relation
        L0 = 1.0_REAL64
        L1 = 1.0_REAL64 + alpha - x
        
        DO i = 2, n
            Ln = ((2.0_REAL64*REAL(i-1,REAL64) + 1.0_REAL64 + alpha - x)*L1 - &
                  (REAL(i-1,REAL64) + alpha)*L0) / REAL(i,REAL64)
            L0 = L1
            L1 = Ln
        END DO
        
        res = Ln
    END FUNCTION assoc_laguerre

    ! Radial wavefunction R_{n_r l}(r)
    FUNCTION R_radial(n_r, l, b, r) RESULT(res)
        INTEGER(KIND=INT32), INTENT(IN) :: n_r, l
        REAL(KIND=REAL64), INTENT(IN) :: b, r
        REAL(KIND=REAL64) :: res
        REAL(KIND=REAL64) :: rho, rho2, alpha, N_norm, L_val
        INTRINSIC :: sqrt, gamma, exp
        
        rho = r / b
        rho2 = rho * rho
        alpha = REAL(l, REAL64) + 0.5_REAL64
        
        N_norm = sqrt(2.0_REAL64 * gamma(REAL(n_r + 1, REAL64)) / &
                      (b**3 * gamma(REAL(n_r + l, REAL64) + 1.5_REAL64)))
        
        L_val = assoc_laguerre(n_r, alpha, rho2)
        
        res = N_norm * (rho**l) * EXP(-rho2 / 2.0_REAL64) * L_val
    END FUNCTION R_radial

    ! Woods-Saxon density profile
    FUNCTION GenDensity(nucl, r) RESULT(dens)
        INTEGER(KIND=INT32), INTENT(IN) :: nucl
        REAL(KIND=REAL64), INTENT(IN) :: r
        REAL(KIND=REAL64) :: dens
        INTRINSIC :: exp
        
        dens = 0.0_REAL64
        SELECT CASE (nucl)
        CASE (1) ! Carbon
            dens = 0.0922_REAL64 / (1.0_REAL64 + EXP((r - 2.861_REAL64) / 0.52_REAL64))
        CASE (2) ! Al
            dens = 0.1027_REAL64 / (1.0_REAL64 + EXP((r - 3.75_REAL64) / 0.52_REAL64))
        CASE (3) ! Fe56
            dens = 0.1095_REAL64 / (1.0_REAL64 + EXP((r - 4.781_REAL64) / 0.52_REAL64))
        CASE (4) ! Pb
            dens = 0.1156_REAL64 / (1.0_REAL64 + EXP((r - 7.423_REAL64) / 0.52_REAL64))
        CASE (5) ! Ca40
            dens = 0.1066_REAL64 / (1.0_REAL64 + EXP((r - 4.274_REAL64) / 0.52_REAL64))
        CASE (6) ! Ca48
            dens = 0.1081_REAL64 / (1.0_REAL64 + EXP((r - 4.543_REAL64) / 0.52_REAL64))
        CASE (7) ! Fe54
            dens = 0.1090_REAL64 / (1.0_REAL64 + EXP((r - 4.724_REAL64) / 0.52_REAL64))
        END SELECT
    END FUNCTION GenDensity

    ! Transparency calculation
    FUNCTION getTransparency(nucl, x, y, z_start, rmax, sigmaNN) RESULT(transp)
        INTEGER(KIND=INT32), INTENT(IN) :: nucl
        REAL(KIND=REAL64), INTENT(IN) :: x, y, z_start, rmax, sigmaNN
        REAL(KIND=REAL64) :: transp
        REAL(KIND=REAL64) :: integral, z, rPrime, dStep, y_sq
        INTRINSIC :: sqrt, exp
        
        dStep = 0.01_REAL64
        z = z_start
        integral = 0.0_REAL64
        y_sq = y*y
        
        DO WHILE (z**2 < rmax**2 - x**2 - y_sq)
            rPrime = sqrt(x**2 + y_sq + z**2)
            integral = integral + GenDensity(nucl, rPrime) * dStep
            z = z + dStep
        END DO
        transp = EXP(-integral * sigmaNN)
    END FUNCTION getTransparency
    
    SUBROUTINE configure_nucleus(nuclType, A, Z, N, b_p, b_n)
        INTEGER(KIND=INT32), INTENT(IN) :: nuclType
        REAL(KIND=REAL64), INTENT(OUT) :: A, b_p, b_n
        INTEGER(KIND=INT32), INTENT(OUT) :: Z, N
        REAL(KIND=REAL64) :: hbar_omega_p, hbar_omega_n
        INTRINSIC :: allocated, sqrt

        SELECT CASE(nuclType)
            CASE(0) ! Deuteron
                A = 2.0; Z = 1; N = 1
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(1), neutrons(1))
                protons = [ShellState(0,0)]; neutrons = [ShellState(0,0)]
            CASE(1) ! 12C
                A = 12.0; Z = 6; N = 6
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(6), neutrons(6))
                protons = [ShellState(0,0), ShellState(0,0), ShellState(0,1), ShellState(0,1), ShellState(0,1), ShellState(0,1)]
                neutrons = protons
            CASE(2) ! 27Al
                A = 27.0; Z = 13; N = 14
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(13), neutrons(14))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2)]
                neutrons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                           ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2)]
            CASE(3) ! 56Fe
                A = 56.0; Z = 26; N = 30
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(26), neutrons(30))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0),&
                          ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                          ShellState(0,3)]
                neutrons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                           ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(1,1),ShellState(1,1)]
            CASE(4) ! 208Pb
                A = 208.0; Z = 82; N = 126
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(82), neutrons(126))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),&
                          ShellState(1,0),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                          ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(1,1),ShellState(1,1),&
                          ShellState(1,1),ShellState(1,1),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                          ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(1,1),ShellState(1,1),&
                          ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                          ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                          ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),&
                          ShellState(1,2),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                          ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(2,0),&
                          ShellState(2,0),ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),&
                          ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),&
                          ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),&
                          ShellState(0,5),ShellState(0,5)]
                neutrons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                           ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),&
                           ShellState(1,0),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(1,1),ShellState(1,1),&
                           ShellState(1,1),ShellState(1,1),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(1,1),ShellState(1,1),&
                           ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                           ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                           ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),&
                           ShellState(1,2),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),&
                           ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(0,4),ShellState(2,0),&
                           ShellState(2,0),ShellState(1,2),ShellState(1,2),ShellState(1,2),ShellState(1,2),&
                           ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),&
                           ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),ShellState(0,5),&
                           ShellState(0,5),ShellState(0,5),ShellState(1,3),ShellState(1,3),ShellState(1,3),&
                           ShellState(1,3),ShellState(1,3),ShellState(1,3),ShellState(1,3),ShellState(1,3),&
                           ShellState(2,1),ShellState(2,1),ShellState(2,1),ShellState(2,1),ShellState(1,3),&
                           ShellState(1,3),ShellState(1,3),ShellState(1,3),ShellState(1,3),ShellState(1,3),&
                           ShellState(2,1),ShellState(2,1),ShellState(0,6),ShellState(0,6),ShellState(0,6),&
                           ShellState(0,6),ShellState(0,6),ShellState(0,6),ShellState(0,6),ShellState(0,6),&
                           ShellState(0,6),ShellState(0,6),ShellState(0,6),ShellState(0,6),ShellState(0,6),&
                           ShellState(0,6),ShellState(2,4),ShellState(2,4),ShellState(2,4),ShellState(2,4),&
                           ShellState(2,4),ShellState(2,4),ShellState(2,4),ShellState(2,4),ShellState(2,4),&
                           ShellState(2,4)]
            CASE(5) ! 40Ca
                A = 40.0; Z = 20; N = 20
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(20), neutrons(20))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0)]
                neutrons = protons
            CASE(6) ! 48Ca
                A = 48.0; Z = 20; N = 28
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(20), neutrons(28))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0)]
                neutrons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                           ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3)]
            CASE(7) ! 54Fe
                A = 54.0; Z = 26; N = 28
                IF(ALLOCATED(protons)) DEALLOCATE(protons, neutrons)
                ALLOCATE(protons(26), neutrons(28))
                protons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                          ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                          ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0),&
                          ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                          ShellState(0,3)]
                neutrons = [ShellState(0,0),ShellState(0,0),ShellState(0,1),ShellState(0,1),ShellState(0,1),&
                           ShellState(0,1),ShellState(0,1),ShellState(0,1),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(0,2),&
                           ShellState(0,2),ShellState(0,2),ShellState(0,2),ShellState(1,0),ShellState(1,0),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),ShellState(0,3),&
                           ShellState(0,3),ShellState(0,3),ShellState(0,3)]
            CASE DEFAULT
                WRITE(*,*) "Unknown nucleus type:", nuclType
                STOP 1
        END SELECT

        hbar_omega_p = 41.0_REAL64 * (REAL(Z, REAL64)**(-1.0_REAL64/3.0_REAL64))
        hbar_omega_n = 41.0_REAL64 * (REAL(N, REAL64)**(-1.0_REAL64/3.0_REAL64))
        b_p = HBARC / sqrt(M_N * hbar_omega_p)
        b_n = HBARC / sqrt(M_N * hbar_omega_n)
    END SUBROUTINE configure_nucleus

END MODULE Procedures


PROGRAM src_main
    USE, INTRINSIC :: iso_fortran_env, ONLY: REAL64, INT32, INT64
    USE omp_lib
    USE GlobalData
    USE PhysicalConstants
    USE Procedures
    
    IMPLICIT NONE
    
    INTEGER(KIND=INT32) :: nuclType, Z, N
    REAL(KIND=REAL64) :: Rsrc, pairProb, A
    CHARACTER(LEN=10) :: arg1, arg2, arg3
    REAL(KIND=REAL64) :: r_max, spacing, sigmaNN, Rsrc_squared
    REAL(KIND=REAL64) :: b_p, b_n, cellVol
    INTEGER(KIND=INT32) :: n_steps_1d, i, j, k, l_loop, tid
    TYPE(Nucleon), ALLOCATABLE :: nucleons(:)
    REAL(KIND=REAL64), ALLOCATABLE :: centers(:,:), center_r(:), center_transparency(:)
    INTEGER(KIND=INT32) :: c_idx, p_count, n_count
    REAL(KIND=REAL64) :: prob, transp, r_val, Rp, Rn
    REAL(KIND=REAL64) :: prob_threshold
    INTEGER(KIND=INT32), PARAMETER :: n_r_max = 4, l_max = 8
    
    ! Performance-related variables
    INTEGER(KIND=INT64), ALLOCATABLE :: cell_offsets(:)
    INTEGER(KIND=INT32), ALLOCATABLE :: flat_nucleons_in_cell(:)
    
    ! Variables for main computation loop
    REAL(KIND=REAL64), ALLOCATABLE :: thread_np(:,:), thread_pp(:,:), thread_nn(:,:)
    REAL(KIND=REAL64), ALLOCATABLE :: thread_P_np(:), thread_P_pp(:), thread_P_nn(:), thread_Reep(:)
    INTEGER(KIND=INT32) :: num_threads, search_radius_cells
    INTEGER(KIND=INT32) :: ix1, iy1, iz1, ix2, iy2, iz2
    INTEGER(KIND=INT64) :: c1_idx, c2_idx, start1, end1, start2, end2
    INTEGER(KIND=INT64) :: loop_i, loop_j, n1_flat_idx, n2_flat_idx, flat_idx
    REAL(KIND=REAL64) :: dx, dy, dz, dist_squared
    
    ! Timing variables
    REAL(KIND=REAL64) :: t_start, t_end

    ! ---- Argument Parsing ----
    IF (command_argument_count() < 3) THEN
        WRITE(*,*) "Usage: ./program NucleusType Rsrc pairProb"
        STOP 1
    END IF
    CALL get_command_argument(1, arg1)
    CALL get_command_argument(2, arg2)
    CALL get_command_argument(3, arg3)
    READ(arg1, *) nuclType
    READ(arg2, *) Rsrc
    READ(arg3, *) pairProb
    
    ! ---- Initialization ----
    t_start = omp_get_wtime()

    ALLOCATE(np(n_r_max, l_max, n_r_max, l_max))
    ALLOCATE(pp(n_r_max, l_max, n_r_max, l_max))
    ALLOCATE(nn(n_r_max, l_max, n_r_max, l_max))
    np = 0.0_REAL64; pp = 0.0_REAL64; nn = 0.0_REAL64

    r_max = 10.0_REAL64
    spacing = 0.3_REAL64
    sigmaNN = 3.5_REAL64
    Rsrc_squared = Rsrc * Rsrc
    
    CALL configure_nucleus(nuclType, A, Z, N, b_p, b_n)
    
    ! Calculate dynamic probability threshold
    BLOCK
        REAL(KIND=REAL64) :: accuracy, total_volume, N_nucleons, N_pairs
        accuracy = 1e-3_REAL64
        total_volume = (2.0_REAL64 * r_max)**3
        N_nucleons = A * total_volume / (spacing**3)
        N_pairs = N_nucleons * (N_nucleons - 1.0_REAL64) / 2.0_REAL64
        prob_threshold = sqrt(accuracy / N_pairs)
        WRITE(*,'(A, ES12.5)') " Probability threshold: ", prob_threshold
    END BLOCK
    
    ! ---- Grid and Center Pre-computation ----
    n_steps_1d = INT(FLOOR(2 * r_max / spacing)) + 1
    c_idx = n_steps_1d**3
    ALLOCATE(centers(3, c_idx))
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    DO i = 1, n_steps_1d
        BLOCK
            REAL(KIND=REAL64) :: x, y, z
            x = -r_max + (i-1) * spacing
            DO j = 1, n_steps_1d
                y = -r_max + (j-1) * spacing
                DO k = 1, n_steps_1d
                    z = -r_max + (k-1) * spacing
                    centers(:, (i-1)*n_steps_1d**2 + (j-1)*n_steps_1d + k) = [x, y, z]
                END DO
            END DO
        END BLOCK
    END DO
    !$OMP END PARALLEL DO
    
    ALLOCATE(center_r(c_idx), center_transparency(c_idx))
    !$OMP PARALLEL DO
    DO i = 1, c_idx
        center_r(i) = NORM2(centers(:,i))
        center_transparency(i) = getTransparency(nuclType, centers(1,i), centers(2,i), centers(3,i), r_max, sigmaNN)
    END DO
    !$OMP END PARALLEL DO
    
    WRITE(*,*) "Creating nucleon probability distribution..."
    
    ! ---- Nucleon Creation ----
    cellVol = spacing**3
    p_count = SIZE(protons)
    n_count = SIZE(neutrons)
    
    BLOCK
        TYPE(Nucleon), ALLOCATABLE :: temp_nucleons(:)
        INTEGER(KIND=INT32) :: current_nucleon_count
        ALLOCATE(temp_nucleons( (p_count + n_count) * c_idx ))
        current_nucleon_count = 0

        !$OMP PARALLEL
        BLOCK
            TYPE(Nucleon), ALLOCATABLE :: local_nucleons(:)
            INTEGER(KIND=INT32) :: local_count, p_idx_block, n_idx_block, c_idx_block
            ALLOCATE(local_nucleons( (p_count + n_count) * c_idx / omp_get_num_threads() + 1000))
            local_count = 0

            ! Protons
            !$OMP DO SCHEDULE(DYNAMIC, 1)
            DO p_idx_block = 1, p_count
                DO c_idx_block = 1, c_idx
                    r_val = center_r(c_idx_block)
                    Rp = R_radial(protons(p_idx_block)%n_r, protons(p_idx_block)%l, b_p, r_val)
                    prob = (Rp*Rp/(4*PI)) * cellVol
                    IF (prob > prob_threshold) THEN
                        local_count = local_count + 1
                        transp = center_transparency(c_idx_block)
                        local_nucleons(local_count) = Nucleon(0, p_idx_block, c_idx_block, prob, transp)
                    END IF
                END DO
            END DO

            ! Neutrons
            !$OMP DO SCHEDULE(DYNAMIC, 1)
            DO n_idx_block = 1, n_count
                DO c_idx_block = 1, c_idx
                    r_val = center_r(c_idx_block)
                    Rn = R_radial(neutrons(n_idx_block)%n_r, neutrons(n_idx_block)%l, b_n, r_val)
                    prob = (Rn*Rn/(4*PI)) * cellVol
                    IF (prob > prob_threshold) THEN
                        local_count = local_count + 1
                        transp = center_transparency(c_idx_block)
                        local_nucleons(local_count) = Nucleon(1, n_idx_block, c_idx_block, prob, transp)
                    END IF
                END DO
            END DO
            
            !$OMP CRITICAL
            temp_nucleons(current_nucleon_count+1 : current_nucleon_count+local_count) = local_nucleons(1:local_count)
            current_nucleon_count = current_nucleon_count + local_count
            !$OMP END CRITICAL
        END BLOCK
        !$OMP END PARALLEL
        ALLOCATE(nucleons(current_nucleon_count))
        nucleons = temp_nucleons(1:current_nucleon_count)
    END BLOCK
    
    WRITE(*,*) "Created ", SIZE(nucleons), " nucleon-position combinations"

    ! ---- Spatial Partitioning Setup ----
    BLOCK
        INTEGER(KIND=INT32), ALLOCATABLE :: temp_nucleon_indices(:)
        INTEGER(KIND=INT32) :: cell_idx, nucleon_idx
        
        ALLOCATE(cell_offsets(c_idx + 1))
        cell_offsets = 0
        
        DO nucleon_idx = 1, SIZE(nucleons)
            cell_idx = nucleons(nucleon_idx)%center_index
            cell_offsets(cell_idx + 1) = cell_offsets(cell_idx + 1) + 1
        END DO
        
        DO cell_idx = 2, c_idx + 1
            cell_offsets(cell_idx) = cell_offsets(cell_idx) + cell_offsets(cell_idx - 1)
        END DO
        
        ALLOCATE(flat_nucleons_in_cell(cell_offsets(c_idx+1)))
        ALLOCATE(temp_nucleon_indices(c_idx))
        temp_nucleon_indices = cell_offsets(1:c_idx)

        DO nucleon_idx = 1, SIZE(nucleons)
            cell_idx = nucleons(nucleon_idx)%center_index
            temp_nucleon_indices(cell_idx) = temp_nucleon_indices(cell_idx) + 1
            flat_nucleons_in_cell(temp_nucleon_indices(cell_idx)) = nucleon_idx
        END DO
    END BLOCK

    WRITE(*,*) "Computing pair probabilities..."

    ! ---- Main Computation Loop ----
    !$OMP PARALLEL
    num_threads = omp_get_num_threads()
    !$OMP SINGLE
        WRITE(*,*) "Running with ", num_threads, " threads."
        ALLOCATE(thread_np(n_r_max*l_max*n_r_max*l_max, num_threads))
        ALLOCATE(thread_pp(n_r_max*l_max*n_r_max*l_max, num_threads))
        ALLOCATE(thread_nn(n_r_max*l_max*n_r_max*l_max, num_threads))
        ALLOCATE(thread_P_np(num_threads), thread_P_pp(num_threads), thread_P_nn(num_threads), thread_Reep(num_threads))
        thread_np = 0.0; thread_pp = 0.0; thread_nn = 0.0
        thread_P_np = 0.0; thread_P_pp = 0.0; thread_P_nn = 0.0; thread_Reep = 0.0
    !$OMP END SINGLE
    
    search_radius_cells = INT(CEILING(Rsrc / spacing))

    !$OMP DO SCHEDULE(DYNAMIC, 1) COLLAPSE(3) &
    !$OMP PRIVATE(ix1,iy1,iz1,ix2,iy2,iz2,tid,c1_idx,c2_idx,start1,end1,start2,end2) &
    !$OMP PRIVATE(loop_i,loop_j,n1_flat_idx,n2_flat_idx,dx,dy,dz,dist_squared)
    DO i = 1, n_steps_1d
        DO j = 1, n_steps_1d
            DO k = 1, n_steps_1d
                tid = omp_get_thread_num()
                ix1 = i - 1; iy1 = j - 1; iz1 = k - 1
                c1_idx = ix1 * n_steps_1d**2 + iy1 * n_steps_1d + iz1 + 1
                
                start1 = cell_offsets(c1_idx) + 1
                end1 = cell_offsets(c1_idx + 1)
                IF (start1 > end1) CYCLE
                
                ! Pairs within the same cell
                DO loop_i = start1, end1
                    DO loop_j = loop_i + 1, end1
                        CALL process_pair(flat_nucleons_in_cell(loop_i), &
                                        flat_nucleons_in_cell(loop_j), tid)
                    END DO
                END DO

                ! Pairs with neighboring cells
                DO ix2 = ix1 - search_radius_cells, ix1 + search_radius_cells
                    DO iy2 = iy1 - search_radius_cells, iy1 + search_radius_cells
                        DO iz2 = iz1 - search_radius_cells, iz1 + search_radius_cells
                            IF (ix2 < 0 .OR. ix2 >= n_steps_1d .OR. &
                                iy2 < 0 .OR. iy2 >= n_steps_1d .OR. &
                                iz2 < 0 .OR. iz2 >= n_steps_1d) CYCLE

                            c2_idx = ix2 * n_steps_1d**2 + iy2 * n_steps_1d + iz2 + 1
                            IF (c2_idx <= c1_idx) CYCLE
                            
                            start2 = cell_offsets(c2_idx) + 1
                            end2 = cell_offsets(c2_idx + 1)
                            IF (start2 > end2) CYCLE

                            dx = centers(1, c1_idx) - centers(1, c2_idx)
                            dy = centers(2, c1_idx) - centers(2, c2_idx)
                            dz = centers(3, c1_idx) - centers(3, c2_idx)
                            dist_squared = dx*dx + dy*dy + dz*dz
                            
                            IF (dist_squared >= Rsrc_squared) CYCLE
                            
                            DO n1_flat_idx = start1, end1
                                DO n2_flat_idx = start2, end2
                                    CALL process_pair(flat_nucleons_in_cell(n1_flat_idx), &
                                                    flat_nucleons_in_cell(n2_flat_idx), tid)
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    !$OMP END DO
    
    ! Final combine of np,pp,nn arrays
    !$OMP DO PRIVATE(flat_idx,tid)
    DO i = 1, n_r_max
        DO j = 1, l_max
            DO k = 1, n_r_max
                DO l_loop = 1, l_max
                    DO tid = 1, num_threads
                        flat_idx = (i-1)*(l_max*n_r_max*l_max) + (j-1)*(n_r_max*l_max) + &
                                  (k-1)*l_max + l_loop-1
                        np(i,j,k,l_loop) = np(i,j,k,l_loop) + thread_np(flat_idx+1, tid)
                        pp(i,j,k,l_loop) = pp(i,j,k,l_loop) + thread_pp(flat_idx+1, tid)
                        nn(i,j,k,l_loop) = nn(i,j,k,l_loop) + thread_nn(flat_idx+1, tid)
                    END DO
                END DO
            END DO
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    t_end = omp_get_wtime()
    
    ! ---- Final Output ----
    BLOCK
        REAL(KIND=REAL64) :: P_np, P_pp, P_nn, Reep, Pp_tot, Pn_tot, src_fraction
        P_np = SUM(thread_P_np)
        P_pp = SUM(thread_P_pp)
        P_nn = SUM(thread_P_nn)
        Reep = SUM(thread_Reep)
        
        Pp_tot = 0.0; Pn_tot = 0.0
        DO i = 1, SIZE(nucleons)
            IF (nucleons(i)%type == 0) THEN
                Pp_tot = Pp_tot + nucleons(i)%probability
            ELSE
                Pn_tot = Pn_tot + nucleons(i)%probability
            END IF
        END DO
        
        src_fraction = (P_np + P_pp + P_nn) / A

        WRITE(*,'(A)') ""
        WRITE(*,'(A)') "Total close-proximity probabilities (small spheres):"
        WRITE(*,'(A, G0)') "  #np: ", P_np
        WRITE(*,'(A, G0)') "  #pp: ", P_pp
        WRITE(*,'(A, G0)') "  #nn: ", P_nn
        WRITE(*,'(A, G0)') "  R(e,e''p): ", Reep
        WRITE(*,'(A, G0)') "rmax: ", r_max
        WRITE(*,'(A, G0)') "SRC fraction = ", src_fraction
        IF (P_np > 0) THEN
            WRITE(*,'(A, G0, A, G0)') "pp/np = ", P_pp/P_np, "  nn/np = ", P_nn/P_np
        ELSE
            WRITE(*,'(A)') "pp/np and nn/np are undefined because P_np is zero."
        END IF
        WRITE(*,'(A, G0, A, G0)') "Pp = ", Pp_tot, " Pn = ", Pn_tot
        WRITE(*,'(A, G0)') "spacing: ", spacing
        WRITE(*,'(A, G0)') "Rsrc = ", Rsrc
        WRITE(*,'(A, G0)') "pairProb = ", pairProb
        WRITE(*,'(A, F8.3, A)') "Total execution time: ", t_end - t_start, " seconds."
    END BLOCK

CONTAINS
    SUBROUTINE process_pair(n1_idx, n2_idx, tid)
        USE GlobalData
        USE PhysicalConstants
        USE Procedures
        IMPLICIT NONE
        INTEGER(KIND=INT32), INTENT(IN) :: n1_idx, n2_idx, tid
        TYPE(Nucleon) :: n1, n2
        REAL(KIND=REAL64) :: Tnp_p, Tnp_n, Tpp1, Tpp2, Tnn1, Tnn2
        REAL(KIND=REAL64) :: np_prob, pp_prob, nn_prob, pair_prob
        INTEGER(KIND=INT32) :: p_idx, n_idx
        TYPE(ShellState) :: p_shell, n_shell, p1_shell, p2_shell, n1_shell, n2_shell
        REAL(KIND=REAL64), PARAMETER :: LimitAngular = 100.0
        INTEGER(KIND=INT64) :: flat_idx

        n1 = nucleons(n1_idx)
        n2 = nucleons(n2_idx)

        IF (n1%type == n2%type .AND. n1%shell_index == n2%shell_index) RETURN
        
        pair_prob = n1%probability * n2%probability
        Tnp_p=0; Tnp_n=0; Tpp1=0; Tpp2=0; Tnn1=0; Tnn2=0
        np_prob=0; pp_prob=0; nn_prob=0

        IF ((n1%type == 0 .AND. n2%type == 1) .OR. (n1%type == 1 .AND. n2%type == 0)) THEN
            IF (n1%type == 0) THEN
                p_idx = n1%shell_index; n_idx = n2%shell_index
                Tnp_p = n1%transparency; Tnp_n = n2%transparency
            ELSE
                p_idx = n2%shell_index; n_idx = n1%shell_index
                Tnp_p = n2%transparency; Tnp_n = n1%transparency
            END IF
            
            p_shell = protons(p_idx); n_shell = neutrons(n_idx)
            IF (ABS(p_shell%l - n_shell%l) <= LimitAngular) THEN
                np_prob = pair_prob
                thread_P_np(tid+1) = thread_P_np(tid+1) + np_prob
                flat_idx = (n_shell%n_r)*(l_max*n_r_max*l_max) + &
                          (n_shell%l)*(n_r_max*l_max) + &
                          (p_shell%n_r)*l_max + p_shell%l
                thread_np(flat_idx+1, tid+1) = thread_np(flat_idx+1, tid+1) + np_prob
            END IF
        ELSE IF (n1%type == 0 .AND. n2%type == 0) THEN
            p1_shell = protons(n1%shell_index); p2_shell = protons(n2%shell_index)
            Tpp1 = n1%transparency; Tpp2 = n2%transparency
            IF (ABS(p1_shell%l - p2_shell%l) <= LimitAngular) THEN
                pp_prob = pair_prob * pairProb
                thread_P_pp(tid+1) = thread_P_pp(tid+1) + pp_prob
                flat_idx = (p1_shell%n_r)*(l_max*n_r_max*l_max) + &
                          (p1_shell%l)*(n_r_max*l_max) + &
                          (p2_shell%n_r)*l_max + p2_shell%l
                thread_pp(flat_idx+1, tid+1) = thread_pp(flat_idx+1, tid+1) + pp_prob
            END IF
        ELSE IF (n1%type == 1 .AND. n2%type == 1) THEN
            n1_shell = neutrons(n1%shell_index); n2_shell = neutrons(n2%shell_index)
            Tnn1 = n1%transparency; Tnn2 = n2%transparency
            IF (ABS(n1_shell%l - n2_shell%l) <= LimitAngular) THEN
                nn_prob = pair_prob * pairProb
                thread_P_nn(tid+1) = thread_P_nn(tid+1) + nn_prob
                flat_idx = (n1_shell%n_r)*(l_max*n_r_max*l_max) + &
                          (n1_shell%l)*(n_r_max*l_max) + &
                          (n2_shell%n_r)*l_max + n2_shell%l
                thread_nn(flat_idx+1, tid+1) = thread_nn(flat_idx+1, tid+1) + nn_prob
            END IF
        END IF

        thread_Reep(tid+1) = thread_Reep(tid+1) + &
            (Tnp_p*np_prob+(Tpp1+Tpp2)*pp_prob)*2.4 + &
            (Tnp_n*np_prob+(Tnn1+Tnn2)*nn_prob)*1.0*0.05

    END SUBROUTINE process_pair

END PROGRAM src_main 