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
    
    ! Structure-of-Arrays for Nucleons for better cache performance
    INTEGER(KIND=INT32), ALLOCATABLE :: nucleon_type(:)
    INTEGER(KIND=INT32), ALLOCATABLE :: nucleon_shell_index(:)
    INTEGER(KIND=INT32), ALLOCATABLE :: nucleon_center_index(:)
    REAL(KIND=REAL64), ALLOCATABLE :: nucleon_probability(:)
    REAL(KIND=REAL64), ALLOCATABLE :: nucleon_transparency(:)
    ! Temporary sorted arrays for reordering
    INTEGER(KIND=INT32), ALLOCATABLE :: s_nucleon_type(:)
    INTEGER(KIND=INT32), ALLOCATABLE :: s_nucleon_shell_index(:)
    INTEGER(KIND=INT32), ALLOCATABLE :: s_nucleon_center_index(:)
    REAL(KIND=REAL64), ALLOCATABLE :: s_nucleon_probability(:)
    REAL(KIND=REAL64), ALLOCATABLE :: s_nucleon_transparency(:)
    
    REAL(KIND=REAL64), ALLOCATABLE :: centers(:,:), center_r(:), center_transparency(:)
    INTEGER(KIND=INT32) :: c_idx, p_count, n_count, total_nucleons
    REAL(KIND=REAL64) :: prob, transp, r_val, Rp, Rn
    REAL(KIND=REAL64), PARAMETER :: prob_threshold = 1e-6
    INTEGER(KIND=INT32), PARAMETER :: n_r_max = 4, l_max = 8
    
    ! Performance-related variables
    INTEGER(KIND=INT64), ALLOCATABLE :: cell_offsets(:)
    
    ! Neighbor list for efficient pair finding
    TYPE :: CellPair
        INTEGER(KIND=INT32) :: c1
        INTEGER(KIND=INT32) :: c2
    END TYPE CellPair
    TYPE(CellPair), ALLOCATABLE :: valid_pairs(:)

    ! Variables for main computation loop
    REAL(KIND=REAL64), ALLOCATABLE :: thread_np(:,:,:,:,:), thread_pp(:,:,:,:,:), thread_nn(:,:,:,:,:)
    INTEGER(KIND=INT32) :: num_threads, search_radius_cells
    INTEGER(KIND=INT32) :: ix1, iy1, iz1, ix2, iy2, iz2
    INTEGER(KIND=INT64) :: c1_idx, c2_idx, start1, end1, start2, end2
    INTEGER(KIND=INT64) :: loop_i, loop_j, n1_idx, n2_idx
    REAL(KIND=REAL64) :: dx, dy, dz, dist_squared
    REAL(KIND=REAL64) :: P_np_tot, P_pp_tot, P_nn_tot, Reep_tot
    
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
    
    ! ---- Grid and Center Pre-computation ----
    n_steps_1d = INT(FLOOR(2 * r_max / spacing)) + 1
    c_idx = n_steps_1d**3
    ALLOCATE(centers(3, c_idx))
    !$OMP PARALLEL DO PRIVATE(i, j, k) SCHEDULE(STATIC)
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
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, c_idx
        center_r(i) = NORM2(centers(:,i))
        center_transparency(i) = getTransparency(nuclType, centers(1,i), centers(2,i), centers(3,i), r_max, sigmaNN)
    END DO
    !$OMP END PARALLEL DO
    
    WRITE(*,*) "Creating nucleon probability distribution..."
    
    ! ---- Nucleon Creation (Optimized Two-Pass Approach) ----
    cellVol = spacing**3
    p_count = SIZE(protons)
    n_count = SIZE(neutrons)

    BLOCK
        INTEGER(KIND=INT32), ALLOCATABLE :: thread_nucleon_counts(:)
        INTEGER(KIND=INT64), ALLOCATABLE :: thread_offsets(:)
        INTEGER(KIND=INT32) :: p_idx, n_idx, c_idx_loop
        INTEGER(KIND=INT32) :: local_count
        
        ! Pass 1: Count nucleons in parallel to avoid critical sections
        !$OMP PARALLEL PRIVATE(tid, local_count, p_idx, c_idx_loop, r_val, Rp, prob, n_idx, Rn)
        !$OMP SINGLE
            num_threads = omp_get_num_threads()
            ALLOCATE(thread_nucleon_counts(num_threads))
            thread_nucleon_counts = 0
        !$OMP END SINGLE
        
        tid = omp_get_thread_num()
        local_count = 0

        ! Protons
        !$OMP DO SCHEDULE(DYNAMIC, 1)
        DO p_idx = 1, p_count
            DO c_idx_loop = 1, c_idx
                r_val = center_r(c_idx_loop)
                Rp = R_radial(protons(p_idx)%n_r, protons(p_idx)%l, b_p, r_val)
                prob = (Rp*Rp/(4*PI)) * cellVol
                IF (prob > prob_threshold) THEN
                    local_count = local_count + 1
                END IF
            END DO
        END DO

        ! Neutrons
        !$OMP DO SCHEDULE(DYNAMIC, 1)
        DO n_idx = 1, n_count
            DO c_idx_loop = 1, c_idx
                r_val = center_r(c_idx_loop)
                Rn = R_radial(neutrons(n_idx)%n_r, neutrons(n_idx)%l, b_n, r_val)
                prob = (Rn*Rn/(4*PI)) * cellVol
                IF (prob > prob_threshold) THEN
                    local_count = local_count + 1
                END IF
            END DO
        END DO
        thread_nucleon_counts(tid + 1) = local_count
        !$OMP END PARALLEL

        ! Serial prefix sum to determine total size and thread offsets
        total_nucleons = SUM(thread_nucleon_counts)
        ALLOCATE(thread_offsets(num_threads+1))
        thread_offsets(1) = 0
        DO i = 1, num_threads
            thread_offsets(i+1) = thread_offsets(i) + thread_nucleon_counts(i)
        END DO

        ALLOCATE(nucleon_type(total_nucleons), nucleon_shell_index(total_nucleons), &
                 nucleon_center_index(total_nucleons), nucleon_probability(total_nucleons), &
                 nucleon_transparency(total_nucleons))

        ! Pass 2: Fill SoA arrays in parallel
        !$OMP PARALLEL PRIVATE(tid, local_count, p_idx, c_idx_loop, r_val, Rp, prob, transp, n_idx, Rn)
        tid = omp_get_thread_num()
        local_count = 0

        ! Protons
        !$OMP DO SCHEDULE(DYNAMIC, 1)
        DO p_idx = 1, p_count
            DO c_idx_loop = 1, c_idx
                r_val = center_r(c_idx_loop)
                Rp = R_radial(protons(p_idx)%n_r, protons(p_idx)%l, b_p, r_val)
                prob = (Rp*Rp/(4*PI)) * cellVol
                IF (prob > prob_threshold) THEN
                    local_count = local_count + 1
                    transp = center_transparency(c_idx_loop)
                    nucleon_type(thread_offsets(tid+1) + local_count) = 0
                    nucleon_shell_index(thread_offsets(tid+1) + local_count) = p_idx
                    nucleon_center_index(thread_offsets(tid+1) + local_count) = c_idx_loop
                    nucleon_probability(thread_offsets(tid+1) + local_count) = prob
                    nucleon_transparency(thread_offsets(tid+1) + local_count) = transp
                END IF
            END DO
        END DO

        ! Neutrons
        !$OMP DO SCHEDULE(DYNAMIC, 1)
        DO n_idx = 1, n_count
            DO c_idx_loop = 1, c_idx
                r_val = center_r(c_idx_loop)
                Rn = R_radial(neutrons(n_idx)%n_r, neutrons(n_idx)%l, b_n, r_val)
                prob = (Rn*Rn/(4*PI)) * cellVol
                IF (prob > prob_threshold) THEN
                    local_count = local_count + 1
                    transp = center_transparency(c_idx_loop)
                    nucleon_type(thread_offsets(tid+1) + local_count) = 1
                    nucleon_shell_index(thread_offsets(tid+1) + local_count) = n_idx
                    nucleon_center_index(thread_offsets(tid+1) + local_count) = c_idx_loop
                    nucleon_probability(thread_offsets(tid+1) + local_count) = prob
                    nucleon_transparency(thread_offsets(tid+1) + local_count) = transp
                END IF
            END DO
        END DO
        !$OMP END PARALLEL
    END BLOCK
    
    WRITE(*,*) "Created ", total_nucleons, " nucleon-position combinations"

    ! ---- Data Reordering for Cache Locality ----
    BLOCK
        INTEGER(KIND=INT64), ALLOCATABLE :: insert_idx(:)
        INTEGER(KIND=INT32) :: nucleon_idx, cell_idx, insert_pos

        ALLOCATE(cell_offsets(c_idx + 1))
        cell_offsets = 0
        
        ! Create histogram of nucleons per cell
        DO nucleon_idx = 1, total_nucleons
            cell_idx = nucleon_center_index(nucleon_idx)
            cell_offsets(cell_idx + 1) = cell_offsets(cell_idx + 1) + 1
        END DO
        
        ! Convert to cumulative sum to get offsets
        DO cell_idx = 2, c_idx + 1
            cell_offsets(cell_idx) = cell_offsets(cell_idx) + cell_offsets(cell_idx - 1)
        END DO
        
        ! Allocate sorted arrays
        ALLOCATE(s_nucleon_type(total_nucleons), s_nucleon_shell_index(total_nucleons), &
                 s_nucleon_center_index(total_nucleons), s_nucleon_probability(total_nucleons), &
                 s_nucleon_transparency(total_nucleons))
        
        ALLOCATE(insert_idx(c_idx))
        insert_idx = cell_offsets(1:c_idx)

        ! Reorder the nucleons based on their cell index
        DO nucleon_idx = 1, total_nucleons
            cell_idx = nucleon_center_index(nucleon_idx)
            insert_pos = insert_idx(cell_idx) + 1
            
            s_nucleon_type(insert_pos) = nucleon_type(nucleon_idx)
            s_nucleon_shell_index(insert_pos) = nucleon_shell_index(nucleon_idx)
            s_nucleon_center_index(insert_pos) = nucleon_center_index(nucleon_idx)
            s_nucleon_probability(insert_pos) = nucleon_probability(nucleon_idx)
            s_nucleon_transparency(insert_pos) = nucleon_transparency(nucleon_idx)
            
            insert_idx(cell_idx) = insert_pos
        END DO
        DEALLOCATE(insert_idx)

        ! Replace original arrays with sorted ones
        CALL MOVE_ALLOC(s_nucleon_type, nucleon_type)
        CALL MOVE_ALLOC(s_nucleon_shell_index, nucleon_shell_index)
        CALL MOVE_ALLOC(s_nucleon_center_index, nucleon_center_index)
        CALL MOVE_ALLOC(s_nucleon_probability, nucleon_probability)
        CALL MOVE_ALLOC(s_nucleon_transparency, nucleon_transparency)
    END BLOCK

    ! ---- Pre-compute Neighbor List ----
    BLOCK
        INTEGER(KIND=INT32), ALLOCATABLE :: thread_pair_counts(:)
        INTEGER(KIND=INT64), ALLOCATABLE :: thread_pair_offsets(:)
        INTEGER(KIND=INT32) :: total_pairs, local_count, c1, c2, ix1, iy1, iz1, ix2, iy2, iz2
        
        search_radius_cells = INT(CEILING(Rsrc / spacing))
        
        !$OMP PARALLEL
        !$OMP SINGLE
            num_threads = omp_get_num_threads()
            ALLOCATE(thread_pair_counts(num_threads))
            thread_pair_counts = 0
        !$OMP END SINGLE
        
        tid = omp_get_thread_num()
        local_count = 0
        !$OMP DO SCHEDULE(DYNAMIC, 1) COLLAPSE(3)
        DO i = 1, n_steps_1d
            DO j = 1, n_steps_1d
                DO k = 1, n_steps_1d
                    ix1 = i - 1; iy1 = j - 1; iz1 = k - 1
                    c1 = ix1 * n_steps_1d**2 + iy1 * n_steps_1d + iz1 + 1
                    DO ix2 = ix1, ix1 + search_radius_cells
                        DO iy2 = iy1 - search_radius_cells, iy1 + search_radius_cells
                            DO iz2 = iz1 - search_radius_cells, iz1 + search_radius_cells
                                IF (ix2 < 0 .OR. ix2 >= n_steps_1d .OR. &
                                    iy2 < 0 .OR. iy2 >= n_steps_1d .OR. &
                                    iz2 < 0 .OR. iz2 >= n_steps_1d) CYCLE
                                
                                c2 = ix2 * n_steps_1d**2 + iy2 * n_steps_1d + iz2 + 1
                                IF (c2 <= c1) CYCLE

                                dx = centers(1, c1) - centers(1, c2)
                                dy = centers(2, c1) - centers(2, c2)
                                dz = centers(3, c1) - centers(3, c2)
                                dist_squared = dx*dx + dy*dy + dz*dz
                                IF (dist_squared < Rsrc_squared) THEN
                                    local_count = local_count + 1
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO
        thread_pair_counts(tid + 1) = local_count
        !$OMP END PARALLEL

        total_pairs = SUM(thread_pair_counts)
        ALLOCATE(valid_pairs(total_pairs))
        ALLOCATE(thread_pair_offsets(num_threads+1))
        thread_pair_offsets(1) = 0
        DO i = 1, num_threads
            thread_pair_offsets(i+1) = thread_pair_offsets(i) + thread_pair_counts(i)
        END DO

        !$OMP PARALLEL
        tid = omp_get_thread_num()
        local_count = 0
        !$OMP DO SCHEDULE(DYNAMIC, 1) COLLAPSE(3)
        DO i = 1, n_steps_1d
            DO j = 1, n_steps_1d
                DO k = 1, n_steps_1d
                    ix1 = i - 1; iy1 = j - 1; iz1 = k - 1
                    c1 = ix1 * n_steps_1d**2 + iy1 * n_steps_1d + iz1 + 1
                    DO ix2 = ix1, ix1 + search_radius_cells
                        DO iy2 = iy1 - search_radius_cells, iy1 + search_radius_cells
                            DO iz2 = iz1 - search_radius_cells, iz1 + search_radius_cells
                                IF (ix2 < 0 .OR. ix2 >= n_steps_1d .OR. &
                                    iy2 < 0 .OR. iy2 >= n_steps_1d .OR. &
                                    iz2 < 0 .OR. iz2 >= n_steps_1d) CYCLE
                                
                                c2 = ix2 * n_steps_1d**2 + iy2 * n_steps_1d + iz2 + 1
                                IF (c2 <= c1) CYCLE

                                dx = centers(1, c1) - centers(1, c2)
                                dy = centers(2, c1) - centers(2, c2)
                                dz = centers(3, c1) - centers(3, c2)
                                dist_squared = dx*dx + dy*dy + dz*dz
                                IF (dist_squared < Rsrc_squared) THEN
                                    local_count = local_count + 1
                                    valid_pairs(thread_pair_offsets(tid+1) + local_count) = CellPair(c1, c2)
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO
        !$OMP END PARALLEL
    END BLOCK

    WRITE(*,*) "Computing pair probabilities..."

    ! ---- Main Computation Loop (Optimized) ----
    P_np_tot = 0.0_REAL64
    P_pp_tot = 0.0_REAL64
    P_nn_tot = 0.0_REAL64
    Reep_tot = 0.0_REAL64

    !$OMP PARALLEL
    
    !$OMP SINGLE
        num_threads = omp_get_num_threads()
        WRITE(*,*) "Running with ", num_threads, " threads."
        ALLOCATE(thread_np(n_r_max, l_max, n_r_max, l_max, num_threads))
        ALLOCATE(thread_pp(n_r_max, l_max, n_r_max, l_max, num_threads))
        ALLOCATE(thread_nn(n_r_max, l_max, n_r_max, l_max, num_threads))
        thread_np = 0.0; thread_pp = 0.0; thread_nn = 0.0
    !$OMP END SINGLE
    
    !$OMP DO SCHEDULE(GUIDED, 8) PRIVATE(tid, c1_idx, start1, end1, loop_i, loop_j) &
    !$OMP REDUCTION(+:P_np_tot, P_pp_tot, P_nn_tot, Reep_tot)
    DO c1_idx = 1, c_idx
        tid = omp_get_thread_num()
        start1 = cell_offsets(c1_idx) + 1
        end1 = cell_offsets(c1_idx + 1)
        IF (start1 > end1) CYCLE

        ! Pairs within the same cell
        DO loop_i = start1, end1
            DO loop_j = loop_i + 1, end1
                CALL process_pair(loop_i, loop_j, tid, &
                                P_np_tot, P_pp_tot, P_nn_tot, Reep_tot)
            END DO
        END DO
    END DO
    !$OMP END DO

    !$OMP DO SCHEDULE(GUIDED, 8) &
    !$OMP PRIVATE(tid,c1_idx,c2_idx,start1,end1,start2,end2,n1_idx,n2_idx) &
    !$OMP REDUCTION(+:P_np_tot, P_pp_tot, P_nn_tot, Reep_tot)
    DO i = 1, SIZE(valid_pairs)
        tid = omp_get_thread_num()
        c1_idx = valid_pairs(i)%c1
        c2_idx = valid_pairs(i)%c2

        start1 = cell_offsets(c1_idx) + 1
        end1 = cell_offsets(c1_idx + 1)
        start2 = cell_offsets(c2_idx) + 1
        end2 = cell_offsets(c2_idx + 1)
        
        IF (start1 > end1 .OR. start2 > end2) CYCLE

        DO n1_idx = start1, end1
            DO n2_idx = start2, end2
                CALL process_pair(n1_idx, n2_idx, tid, &
                                P_np_tot, P_pp_tot, P_nn_tot, Reep_tot)
            END DO
        END DO
    END DO
    !$OMP END DO

    ! Final combine of np,pp,nn arrays
    !$OMP DO
    DO i = 1, n_r_max
        DO j = 1, l_max
            DO k = 1, n_r_max
                DO l_loop = 1, l_max
                    np(i,j,k,l_loop) = SUM(thread_np(i, j, k, l_loop, :))
                    pp(i,j,k,l_loop) = SUM(thread_pp(i, j, k, l_loop, :))
                    nn(i,j,k,l_loop) = SUM(thread_nn(i, j, k, l_loop, :))
                END DO
            END DO
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    t_end = omp_get_wtime()
    
    ! ---- Final Output ----
    BLOCK
        REAL(KIND=REAL64) :: Pp_tot_calc, Pn_tot_calc, src_fraction
        
        Pp_tot_calc = 0.0; Pn_tot_calc = 0.0
        DO i = 1, total_nucleons
            IF (nucleon_type(i) == 0) THEN
                Pp_tot_calc = Pp_tot_calc + nucleon_probability(i)
            ELSE
                Pn_tot_calc = Pn_tot_calc + nucleon_probability(i)
            END IF
        END DO
        
        src_fraction = (P_np_tot + P_pp_tot + P_nn_tot) / A

        WRITE(*,'(A)') ""
        WRITE(*,'(A)') "Total close-proximity probabilities (small spheres):"
        WRITE(*,'(A, G0)') "  #np: ", P_np_tot
        WRITE(*,'(A, G0)') "  #pp: ", P_pp_tot
        WRITE(*,'(A, G0)') "  #nn: ", P_nn_tot
        WRITE(*,'(A, G0)') "  R(e,e''p): ", Reep_tot
        WRITE(*,'(A, G0)') "rmax: ", r_max
        WRITE(*,'(A, G0)') "SRC fraction = ", src_fraction
        IF (P_np_tot > 0) THEN
            WRITE(*,'(A, G0, A, G0)') "pp/np = ", P_pp_tot/P_np_tot, "  nn/np = ", P_nn_tot/P_np_tot
        ELSE
            WRITE(*,'(A)') "pp/np and nn/np are undefined because P_np is zero."
        END IF
        WRITE(*,'(A, G0, A, G0)') "Pp = ", Pp_tot_calc, " Pn = ", Pn_tot_calc
        WRITE(*,'(A, G0)') "spacing: ", spacing
        WRITE(*,'(A, G0)') "Rsrc = ", Rsrc
        WRITE(*,'(A, G0)') "pairProb = ", pairProb
        WRITE(*,'(A, F8.3, A)') "Total execution time: ", t_end - t_start, " seconds."
    END BLOCK

CONTAINS
    SUBROUTINE process_pair(n1_idx, n2_idx, tid, P_np, P_pp, P_nn, Reep)
        USE GlobalData
        USE PhysicalConstants
        USE Procedures
        IMPLICIT NONE
        INTEGER(KIND=INT64), INTENT(IN) :: n1_idx, n2_idx
        INTEGER(KIND=INT32), INTENT(IN) :: tid
        REAL(KIND=REAL64), INTENT(INOUT) :: P_np, P_pp, P_nn, Reep
        
        INTEGER(KIND=INT32) :: n1_type, n2_type, n1_shell_idx, n2_shell_idx
        REAL(KIND=REAL64) :: n1_prob, n2_prob, n1_transp, n2_transp

        REAL(KIND=REAL64) :: Tnp_p, Tnp_n, Tpp1, Tpp2, Tnn1, Tnn2
        REAL(KIND=REAL64) :: np_prob, pp_prob, nn_prob, pair_prob
        INTEGER(KIND=INT32) :: p_idx, n_idx
        TYPE(ShellState) :: p_shell, n_shell, p1_shell, p2_shell, n1_shell, n2_shell
        REAL(KIND=REAL64), PARAMETER :: LimitAngular = 100.0

        n1_type = nucleon_type(n1_idx)
        n2_type = nucleon_type(n2_idx)
        n1_shell_idx = nucleon_shell_index(n1_idx)
        n2_shell_idx = nucleon_shell_index(n2_idx)

        IF (n1_type == n2_type .AND. n1_shell_idx == n2_shell_idx) RETURN
        
        n1_prob = nucleon_probability(n1_idx)
        n2_prob = nucleon_probability(n2_idx)
        n1_transp = nucleon_transparency(n1_idx)
        n2_transp = nucleon_transparency(n2_idx)

        pair_prob = n1_prob * n2_prob
        Tnp_p=0; Tnp_n=0; Tpp1=0; Tpp2=0; Tnn1=0; Tnn2=0
        np_prob=0; pp_prob=0; nn_prob=0

        IF ((n1_type == 0 .AND. n2_type == 1) .OR. (n1_type == 1 .AND. n2_type == 0)) THEN
            IF (n1_type == 0) THEN
                p_idx = n1_shell_idx; n_idx = n2_shell_idx
                Tnp_p = n1_transp; Tnp_n = n2_transp
            ELSE
                p_idx = n2_shell_idx; n_idx = n1_shell_idx
                Tnp_p = n2_transp; Tnp_n = n1_transp
            END IF
            
            p_shell = protons(p_idx); n_shell = neutrons(n_idx)
            IF (ABS(p_shell%l - n_shell%l) <= LimitAngular) THEN
                np_prob = pair_prob
                P_np = P_np + np_prob
                thread_np(n_shell%n_r+1, n_shell%l+1, p_shell%n_r+1, p_shell%l+1, tid+1) = &
                    thread_np(n_shell%n_r+1, n_shell%l+1, p_shell%n_r+1, p_shell%l+1, tid+1) + np_prob
            END IF
        ELSE IF (n1_type == 0 .AND. n2_type == 0) THEN
            p1_shell = protons(n1_shell_idx); p2_shell = protons(n2_shell_idx)
            Tpp1 = n1_transp; Tpp2 = n2_transp
            IF (ABS(p1_shell%l - p2_shell%l) <= LimitAngular) THEN
                pp_prob = pair_prob * pairProb
                P_pp = P_pp + pp_prob
                thread_pp(p1_shell%n_r+1, p1_shell%l+1, p2_shell%n_r+1, p2_shell%l+1, tid+1) = &
                    thread_pp(p1_shell%n_r+1, p1_shell%l+1, p2_shell%n_r+1, p2_shell%l+1, tid+1) + pp_prob
            END IF
        ELSE IF (n1_type == 1 .AND. n2_type == 1) THEN
            n1_shell = neutrons(n1_shell_idx); n2_shell = neutrons(n2_shell_idx)
            Tnn1 = n1_transp; Tnn2 = n2_transp
            IF (ABS(n1_shell%l - n2_shell%l) <= LimitAngular) THEN
                nn_prob = pair_prob * pairProb
                P_nn = P_nn + nn_prob
                thread_nn(n1_shell%n_r+1, n1_shell%l+1, n2_shell%n_r+1, n2_shell%l+1, tid+1) = &
                    thread_nn(n1_shell%n_r+1, n1_shell%l+1, n2_shell%n_r+1, n2_shell%l+1, tid+1) + nn_prob
            END IF
        END IF

        Reep = Reep + &
            (Tnp_p*np_prob+(Tpp1+Tpp2)*pp_prob)*2.4 + &
            (Tnp_n*np_prob+(Tnn1+Tnn2)*nn_prob)*1.0*0.05

    END SUBROUTINE process_pair

END PROGRAM src_main 