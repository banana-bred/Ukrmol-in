! Copyright 2019
!
! For a comprehensive list of the developers that contributed to these codes
! see the UK-AMOR website.
!
! This file is part of UKRmol-in (UKRmol+ suite).
!
!     UKRmol-in is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     UKRmol-in is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  UKRmol-in (in source/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.

!> \brief   Contracted Hamiltonian module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module handles the calculation of the contracted scattering hamiltonian
!> through the usage of the Contracted_Hamiltonian class. Eq. references refer to the paper by
!> J. Tennyson:  J. Phys B.: At. Mol. Opt. Phys. 29 (1996) 1817-1828
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial Revision
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Contracted_Hamiltonian_module

    use const_gbl,                  only: stdout
    use consts_mpi_ci,              only: NORMAL_PHASE, NO_PURE_TARGET_INTEGRAL_DIFF, MAT_SPARSE, NIDX
    use mpi_gbl,                    only: mpi_reduceall_max
    use precisn,                    only: longint, wp
    use BaseMatrix_module,          only: BaseMatrix
    use Contracted_Symbolic_module, only: ContractedSymbolicElementVector
    use Hamiltonian_module,         only: BaseHamiltonian
    use MemoryManager_module,       only: master_memory
    use Parallelization_module,     only: grid => process_grid
    use Symbolic_module,            only: SymbolicElementVector
    use Target_RMat_CI_module,      only: Target_Rmat_CI
    use Timing_Module,              only: master_timer
    use Utility_module,             only: triangular_index_to_ij, compute_total_triangular, compute_total_box, box_index_to_ij, &
                                          pack_ints, unpack_ints

    implicit none

    public Contracted_Hamiltonian

    private

    !> \brief   Computation of Hamiltonian
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This class computes the hamitonian using the contracted prototype scheme
    !> given by J. Tennyson:  J. Phys B.: At. Mol. Opt. Phys. 29 (1996) 1817-1828.
    !> The notation i1,n1,j1,i2,n2,j2 referes to the variables i,n,j,i',n',j' given in table 1.
    !> The class refer to the ones desribed in the Paper in Section 3.
    !>
    type, extends(BaseHamiltonian) :: Contracted_Hamiltonian

        class(Target_RMat_CI), pointer :: rmat_ci(:)    !< Our Target coefficients
        integer :: start_L2                             !< Start of the L2 CSFS
        integer :: total_num_csf                        !< Total number of CSFs we're dealing with
        integer :: num_L2_functions
        integer :: L2_mat_offset                        !< Matrix starting position of the L2 portion
        integer :: csf_skip                             !< Used to move to the next CSF (most continuum ones come in pairs so the next one is skipped usually)
        integer :: non_zero_prototype = 0               !< Counter for calculated prototypes that are non_zero

    contains

        procedure, public :: initialize        => initialize_contracted_hamiltonian
        procedure, public :: build_hamiltonian => build_contracted_hamiltonian

        !----------class routines----------------------!
        procedure, private :: evaluate_class_1
        procedure, private :: evaluate_class_2
        procedure, private :: evaluate_class_3
        procedure, private :: evaluate_class_4
        procedure, private :: evaluate_class_5
        procedure, private :: evaluate_class_6
        procedure, private :: evaluate_class_7
        procedure, private :: evaluate_class_8
        procedure, private :: evaluate_class_2_and_8
        !--------class routines------------------------!
       !procedure, private :: evaluate_class_7
        !---Contraction routines
        procedure, private :: fast_contract_class_1_3_diag
        procedure, private :: fast_contract_class_1_3_offdiag
        procedure, private :: fast_contract_class_567

        !---------Expansion routines-------------------!
        procedure, private :: expand_diagonal
        procedure, private :: expand_off_diagonal
        procedure, private :: expand_continuum_L2
        procedure, private :: expand_off_diagonal_general

        procedure, private :: expand_continuum_L2_eval
        procedure, private :: expand_diagonal_eval
        procedure, private :: contr_expand_diagonal_eval
        procedure, private :: contr_expand_off_diagonal_eval
        procedure, private :: contr_expand_continuum_L2_eval
        procedure, private :: expand_off_diagonal_eval
        procedure, private :: expand_off_diagonal_gen_eval

        !-------auxilary--routines-----------------!
        procedure, private :: compute_matrix_ij
        procedure, private :: get_target_coeff
        procedure, private :: reduce_prototypes
        procedure, private :: get_starting_index

    end type Contracted_Hamiltonian

contains

    !> \brief   Initializes the contracted hamiltonian and supplies target coefficients
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this  Hamiltonian object to update.
    !> \param[in] rmat_ci  Our Target_RMat_CI array containing the cofficients.
    !>
    subroutine initialize_contracted_hamiltonian (this, rmat_ci)
        class(Contracted_Hamiltonian)             :: this
        class(Target_Rmat_CI), target, intent(in) :: rmat_ci(:)

        !Assign the rmat pointer
        this % rmat_ci => rmat_ci
        !Our skip value is 2 for the moment
        this % csf_skip = 2
        !Where the L2 starts
        this % start_L2 = this % options % last_continuum + 1
        !total number of CSFS
        this % total_num_csf = this % options % num_csfs
        !Calculate where each L2 CSF belongs in the matrix
        this % L2_mat_offset = this % options % contracted_mat_size - this % total_num_csf
        this % num_L2_functions = this % options % num_L2_CSF

    end subroutine initialize_contracted_hamiltonian


    !> \brief   Builds the Contracted Hamiltonian class by class
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this             Hamiltonian object to update.
    !> \param[out]   matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine build_contracted_hamiltonian (this, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        integer :: num_target_sym, target_sym, num_targets_per_sym, i1, i2, n1, n2, l1, loop_skip, loop_ido, my_idx, ido
        integer :: num_csfs, m, n, num_continuum, total_l1_cont_elems, matrix_block_size, estimated_mat_size

        !The block size is the size of the continuum and continuum+L2 section of the matrix
        matrix_block_size = this % options % contracted_mat_size - this % total_num_csf + this % options % last_continuum

        write (stdout, "('continuum block_size =',i8)") matrix_block_size
        !We tell the matrix storage class the contracted size, the block size and the fact that it is sparse
        call matrix_elements % initialize_matrix_structure(this % options % contracted_mat_size, MAT_SPARSE, matrix_block_size)

        !Get the number of target symmetries
        num_target_sym = this % options % num_target_sym

        !------------------------------------------------------------------------------------------------------------------!
        !-----------------------------CONTINUUM CALCULATIONS---------------------------------------------------------------!
        !------------------------------------------------------------------------------------------------------------------!

        this % orbitals % MFLG = 0

        !-----------------------------------------------CLASS 1---------------------------------------------------------
        call master_timer % start_timer("Class 1")

        !Perform a class one ('diagonal same target symmetry') calculation on all target symmetries
        do i1 = 1, num_target_sym
            call this % evaluate_class_1(i1, this % options % num_target_state_sym(i1), matrix_elements)
        end do
        call master_timer % stop_timer("Class 1")
        write (stdout, "('Class 1 complete')")
        flush(stdout)
        call master_timer % report_timers
        !stop
        !-----------------------------------------------CLASS 3---------------------------------------------------------


        !-----------------------------------------------CLASSES 567---------------------------------------------------------
        !These classes deal with differing target symmetries
        !call master_timer%start_timer("Class 56")
        !Loop through every combination of target symmetry
        !do i1 = 1, num_target_sym-1
        !    do i2=i1+1,num_target_sym
        !        if(i1==i2) cycle !We've already done same target symmetries so ignore
                !Check if we have the same continuum symmetry
        !        if(this%options%lambda_continuum_orbitals_target(i1) == this%options%lambda_continuum_orbitals_target(i2)) then
        !            if(this%options%gerade_sym_continuum_orbital_target(i1) == this%options%gerade_sym_continuum_orbital_target(i2)) then

        !                cycle
        !            endif
        !        endif
        !    enddo
        !enddo
        !call master_timer%stop_timer("Class 56")

        call master_memory % print_memory_report

        !Before moving to the pure L2 calculations we need to clear out any other continuum calculations from the matrix format
        !call matrix_elements%update_continuum(.true.)

        !--------------------------------------------------------------------------------------------------------------------------!
        !----------------------------------------------PURE L2 Calculations--------------------------------------------------------!
        !--------------------------------------------------------------------------------------------------------------------------!

        this % orbitals % MFLG = 1

        call master_timer % start_timer("Class 3")
        !Perform a class three ('off-diagonal same target symmetry') calculation on all target symmetries
        do i1 = 1, num_target_sym
            call this % evaluate_class_3(i1, this % options % num_target_state_sym(i1), matrix_elements)
        end do
        call master_timer % stop_timer("Class 3")

        !call master_timer%report_timers

        write (stdout, "('Class 3 complete')")
        flush(stdout)

        call master_timer % report_timers
        call  master_memory % print_memory_report

        !-----------------------------------------------CLASS 4---------------------------------------------------------
        !Class four is a difficult one as we have to loop through the entire L2 functions which are quite large
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        call master_timer % start_timer("Class 4")
        do i1 = 1, num_target_sym
            total_l1_cont_elems = this % options % num_target_state_sym(i1) * this % options % num_continuum_orbitals_target(i1)
            do loop_ido = this % start_L2, this % total_num_csf, loop_skip
                l1 = my_idx + loop_ido

                if (l1 > this % total_num_csf) then
                    !This is dummy to ensure synchornization with other procs
                    do ido = 1, total_l1_cont_elems
                        call matrix_elements % update_pure_L2(.false.)
                    end do
                else
                    !Class four (continuum-L2 calculation)
                    call this % evaluate_class_4(i1, this % options % num_target_state_sym(i1), l1, matrix_elements)
                end if
            end do
        end do

        call master_timer % stop_timer("Class 4")
        write (stdout, "('Class 4 complete')")
        flush(stdout)

        call master_timer % report_timers
        !call matrix_elements%update_pure_L2(.true.)
        !call master_memory%print_memory_report

        !-----------------------------------------------CLASSES 2+8---------------------------------------------------------
        call master_timer % start_timer("Class 567")

        do i1 = 1, num_target_sym - 1
            do i2 = i1 + 1, num_target_sym

                if (i1 == i2) cycle !We've already done same target symmetries so ignore
                !write(stdout,"('class 567: ',2i8)") i1,i2
                !Check if we have the same continuum symmetry
                if (this % options % lambda_continuum_orbitals_target(i1) == &
                    this % options % lambda_continuum_orbitals_target(i2)) then
                    if (this % options % gerade_sym_continuum_orbital_target(i1) == &
                        this % options % gerade_sym_continuum_orbital_target(i2)) then
                        !If we do then do class five ('diagonal same continuum symmetry') calculation
                        !        write(stdout,"('class 5')")
                        call this % evaluate_class_5(i1, this % options % num_target_state_sym(i1), i2, &
                                                     this % options % num_target_state_sym(i2), matrix_elements)
                        !Then a class six ('off-diagonal same continuum symmetry') calculation
                        !        write(stdout,"('class 6')")
                        call this % evaluate_class_6 (i1, this % options % num_target_state_sym(i1), i2, &
                                                      this % options % num_target_state_sym(i2), matrix_elements)
                        !call master_timer%report_timers
                        !call master_memory%print_memory_report
                        cycle
                    end if
                end if

                !write(stdout,"('class 7')")
                !Otherwise do a class seven ('differing continuum symmetry') calculation
                call this % evaluate_class_7(i1, this % options % num_target_state_sym(i1), i2, &
                                             this % options % num_target_state_sym(i2), matrix_elements)
                !call master_timer%report_timers
                !call master_memory%print_memory_report

            end do
        end do

        call master_timer % stop_timer("Class 567")
        write (stdout, "('Class 567 complete')")
        flush(stdout)

        call master_timer % report_timers
        call master_memory % print_memory_report
        call master_timer % start_timer("Class 2+8")
        !Here I've combined both class two ('L2 diagonal') and class eight ('L2 offdiagonal') calculations
        !call this%evaluate_class_2(matrix_elements)
        call this % evaluate_class_2_and_8(matrix_elements)
        !call this%evaluate_class_8(matrix_elements)
        call master_timer % stop_timer("Class 2+8")

        write (stdout, "('Class 2+8 complete')")
        flush(stdout)

        call master_timer % report_timers
        call master_memory % print_memory_report

        !Clear up remaining L2 elements in a matrix formats stash
        call matrix_elements % update_pure_L2(.true.)

        !The matrix may have some last step thats required
        call  matrix_elements % finalize_matrix

        !Print out interesting information
        n = (this % options % num_csfs * (this % options % num_csfs + 1)) / 2
        m = (this % options % contracted_mat_size * (this % options % contracted_mat_size + 1)) / 2

        write (stdout, 6000) m, matrix_elements % get_size()
        6000 format(//,' TOTAL H Matrix     ELEMENTS      =',I15,/,' NON-ZERO ELEMENTS evaluated      =',I15)
        write (stdout, 6010) n, this % non_zero_prototype
        6010 format(' Total    prototype ELEMENTS      =',I15,/,' Non-zero prototype ELEMENTS      =',I15)
        ! WRITE(stdout,6020)nint0, nint00
        6020 format(/,' Number (prototype) integrals     =',I15,/,' Compressed (prototype) integrals =',I15)
        write (stdout, 6030) this % number_of_integrals
        6030 format(' Number integrals evaluated       =',I15)

    end subroutine build_contracted_hamiltonian


    !> \brief   Builds the class 1 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to the building of diagonal elements of the target states with the same target symmetries.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               The target symmetry.
    !> \param[in]  num_targets      The number of target states withing the symmetry.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine evaluate_class_1 (this, i1, num_targets, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        integer,           intent(in)    :: i1
        integer,           intent(in)    :: num_targets

        !>These are the contracted prototypes for each combination of target states
        type(ContractedSymbolicElementVector) :: master_prototype
        !>This stores te resultant slater rule calculation of prototypes before contracting for each state
        type(SymbolicElementVector)           :: temp_prototype
        !>The starting matrix index for each target state combination
        integer, allocatable                  :: matrix_index(:,:)

        !>The number of target state combinations
        integer               :: num_ci
        !------------------------------target state variables
        integer               :: n1, n2
        !----------------------------Looping variables------------------------------!
        integer               :: starting_index, ido, jdo, ii, jj, ij, err
        !----------------------------MPI looping variables--------------------------
        integer               :: loop_idx, loop_ido, my_idx, loop_skip
        !>The total number of configurations with this target symmetry
        integer               :: num_configurations
        !>CSF indicies
        integer               :: config_idx_a, config_idx_b
        !The target coefficient
        real(wp), allocatable :: alpha_coeffs(:,:)
        !The resultant matrix coefficient
        real(wp), allocatable :: mat_coeffs(:,:)
        !Number of continuum orbitals
        integer               :: continuum_orbital, num_continuum_orbitals

        allocate(matrix_index(2,num_targets*(num_targets+1)/2), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS1::MATRIX_INDEX')
        allocate(alpha_coeffs(num_targets*(num_targets+1)/2,1), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS1::alpha_coeffs')
        allocate(mat_coeffs(num_targets*(num_targets+1)/2,1), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS1::mat_coeffs')
        !Get number of continuum orbitals
        num_continuum_orbitals = this % options % num_continuum_orbitals_target(i1)

        !number of target state combinations (triangular form)
        num_ci = num_targets * (num_targets + 1) / 2

        ii = 1
        do n1 = 1, this % options % num_target_state_sym(i1)
            do n2 = n1, this % options % num_target_state_sym(i1)
                !Construct a master prototype for each target state combo
                !And get the starting matrix index
                call this % compute_matrix_ij(i1, n1, 1, i1, n2, 1, matrix_index(1,ii), matrix_index(2,ii))
                ii = ii + 1
            end do
        end do

        call master_prototype % construct(num_ci, 1)

        !Construct the temporary prototype
        call temp_prototype % construct

        !get the starting index of the configurations
        starting_index = this % get_starting_index(i1)

        !Get the number of configurations
        num_configurations = this % options % num_ci_target_sym(i1)

        !Set up our mpi variables
        loop_skip = max(1, grid % lprocs)
        my_idx = max(grid % lrank, 0)
        call master_timer % start_timer("Class 1 prototype")
        !This loop is commonly found in a lot of the classes. Basically it involves striding the loop across mpi processes
        !For example if nprocs is 5, process 0 will do indicies 1,6,11..etc while process 1 will do 2,7,12.. and so on
        do loop_ido = 1, num_configurations, loop_skip
            this % orbitals % MFLG = 0
            !Calculate our real index
            ido = loop_ido + my_idx
            !If we've gon past (which is possible in this setup then skip
            if (ido > num_configurations) cycle
            !Calculate our diagonal configuration index
            config_idx_a = starting_index + (ido - 1) * this % csf_skip
            !Slater it to get our temporary prototypes
            call this % slater_rules(this % csfs(config_idx_a), this % csfs(config_idx_a), temp_prototype, 0, .false.)
            alpha_coeffs = 0.0_wp
            !Loop through each target state
            ii = 1
            do n1 = 1, this % options % num_target_state_sym(i1)
                do n2 = n1, this % options % num_target_state_sym(i1)
                    ! Apply our coefficeints for each combination to the master, Left hand of Eq. 6
                    alpha_coeffs(ii,1) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i1, ido, n2)
                    ii = ii + 1
                end do
            end do
            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)
            !call this%fast_contract_class_1_3_diag(i1,num_targets,ido,temp_prototype,master_prototype)

            !Clear the temporary prototype for usage again
            call temp_prototype % clear

            !Now the off diagonal
            do jdo = 1, ido - 1

                !Get the index
                config_idx_b = starting_index + (jdo - 1) * this % csf_skip
                !Slater
                call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 1, .false.)

                !No prototypes? then we leave
                if (temp_prototype % is_empty()) cycle

                !Now contract the 'off diagonal' diagonal
                ii = 1
                do n1 = 1, this % options % num_target_state_sym(i1)
                    do n2 = n1, this % options % num_target_state_sym(i1)
                        ! (cimn*cim'n' + cim'n*cimn')H_imj,imj' Right hand of Eq. 6
                        alpha_coeffs(ii,1) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i1, jdo, n2) &
                                           + this % get_target_coeff(i1, jdo, n1) * this % get_target_coeff(i1, ido, n2)
                        ! call master_prototype(ii)%add_symbols(temp_prototype,alpha)
                        ii = ii + 1
                    end do
                end do
                !print *,temp_prototype%get_size()
                !print *,jdo
                call master_prototype % add_symbols(temp_prototype, alpha_coeffs)
                !call this%fast_contract_class_1_3_offdiag(i1,num_targets,ido,jdo,temp_prototype,master_prototype)
                !Clear the prototype for the next round
                call temp_prototype % clear
            end do
        end do
        call master_timer % stop_timer("Class 1 prototype")
        ij = 1

        call master_prototype % synchronize_symbols()

        !Without expansion, evaluate the first elements
        !do ii=1,num_targets
        !    do jj=ii,num_targets
        !        mat_coeff = this%evaluate_integrals(master_prototype(ij),NORMAL_PHASE)
                !Insert
                !these are special as they may add the eigenvalues after
                !We just let the master rank add them since they'll be combined later on
        !        if(myrank == master .and. this%diagonal_flag == NO_PURE_TARGET_INTEGRAL_DIFF) mat_coeff = mat_coeff + this%rmat_ci(i1)%eigenvalues(ii)
        !        call matrix_elements%insert_matrix_element(matrix_index(2,ij),matrix_index(1,ij),mat_coeff,1)
                !Possibly perform update (who knows really?)
        !        call matrix_elements%update_continuum(.false.)
        !        ij=ij+1
        !    enddo
        !enddo

        !now we have our prototype we can start the expansion
        config_idx_a = starting_index + (num_configurations - 1) * this % csf_skip
        !Get which continuum orbital this belongs to
        continuum_orbital = this % csfs(config_idx_a) % orbital_sequence_number

        !Set up our mpi variables
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        call master_timer % start_timer("Class 1 Expansion")
        !Now loop through our J
        do loop_ido = 1, num_continuum_orbitals, loop_skip
            ido = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false., num_targets * (num_targets - 1))

            if (ido > num_continuum_orbitals) cycle
            !Loop through each target combo
            mat_coeffs = 0
            !Immediately expand the correct prototype and then immediately evaluate the matrix element
            call this % contr_expand_diagonal_eval(continuum_orbital, master_prototype, ido, mat_coeffs)
            ij = 1
            do ii = 1, num_targets
                do jj = ii, num_targets
                    !Insert into our matrix format
                    if (this % diagonal_flag == NO_PURE_TARGET_INTEGRAL_DIFF .and. ii == jj) then
                        mat_coeffs(ij,1) = mat_coeffs(ij,1) + this % rmat_ci(i1) % eigenvalues(ii)
                    end if
                    call matrix_elements % insert_matrix_element(matrix_index(2,ij) + ido - 1, &
                                                                 matrix_index(1,ij) + ido - 1, mat_coeffs(ij,1), 8)
                    !An update may occur
                    ij = ij + 1
                end do
            end do
        end do

        call master_timer % stop_timer("Class 1 Expansion")
        !Destroy our temporary prototype
        call temp_prototype %destroy
        !Destroy all of the masters
        !do ii=1,num_ci
        call master_prototype % destroy
        !enddo
        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8, size(mat_coeffs))
        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

        this % non_zero_prototype = this % non_zero_prototype + 1

        !call master_timer%report_timers
    end subroutine evaluate_class_1


    !> \brief   Builds the class 2 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This has be deprecated, instead class 2 and 8 have been combined into a faster method.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine evaluate_class_2 (this, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        type(SymbolicElementVector)      :: symbolic_elements
        real(wp) :: mat_coeff
        integer  :: l1, mat_offset, loop_ido, loop_skip, my_idx

        this % diagonal_flag = 0

        !Construct our symbolic elements
        call symbolic_elements % construct

        loop_skip = max(1, grid % gprocs)
        my_idx = max(1, grid % grank)

        !Loop through the L2 functions
        do l1 = this % start_L2, this % total_num_csf
            call matrix_elements % update_pure_L2(.false.)
            !Slater them for the prototypes
            call this % slater_rules(this % csfs(l1), this % csfs(l1), symbolic_elements, 0)
            !Empty? then leave!
            if (symbolic_elements % is_empty()) cycle
            !Evaluate the integral
            mat_coeff = this % evaluate_integrals(symbolic_elements, NORMAL_PHASE)
            !Insert using an offset
            call matrix_elements % insert_matrix_element(l1 + this % L2_mat_offset, l1 + this % L2_mat_offset, mat_coeff, 2)
            !Clear the symbols for the next step
            call symbolic_elements % clear

            this % non_zero_prototype = this % non_zero_prototype + 1
        end do

        !Destroy the symbols
        call symbolic_elements % destroy

    end subroutine evaluate_class_2


    !> \brief   Builds the class 3 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to building the off-diagonal elements of the target states with the same target symmetries.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               The target symmetry.
    !> \param[in]  num_targets      The number of target states withing the symmetry.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine evaluate_class_3 (this, i1, num_targets, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        integer,           intent(in)    :: i1, num_targets

        !>These are the contracted prototypes for each combination of target states
        type(ContractedSymbolicElementVector) :: master_prototype
        !>This stores te resultant slater rule calculation of prototypes before contracting for each state
        type(SymbolicElementVector)           :: temp_prototype
        !>The starting matrix index for each target state combination
        integer, allocatable                  :: matrix_index(:,:)
        !>The number of target state combinations
        integer :: num_ci
        !------------------------------target state variables
        integer :: n1, n2
        !----------------------------Looping variables------------------------------!
        integer :: starting_index, ido, jdo, ii
        !----------------------------MPI looping variables--------------------------
        integer :: loop_idx, loop_ido, my_idx, loop_skip, err
        !>The total number of configurations with this target symmetry
        integer :: num_configurations
        !>CSF indicies
        integer :: config_idx_a, config_idx_b
        !The target coefficient
        real(wp), allocatable :: alpha_coeffs(:,:)
        !The resultant matrix coefficient
        real(wp), allocatable :: mat_coeffs(:,:)
        !Number of continuum orbitals
        integer :: continuum_orbital_a, continuum_orbital_b, num_continuum_orbitals, total_orbs

        allocate(matrix_index(2,num_targets*(num_targets+1)/2), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS3::MATRIX_INDEX')

        allocate(alpha_coeffs(num_targets*(num_targets+1)/2,1), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS3::alpha_coeffs')

        allocate(mat_coeffs(num_targets*(num_targets+1)/2,1), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS3::mat_coeffs')

        !Get number of continuum orbitals
        num_continuum_orbitals = this % options % num_continuum_orbitals_target(i1)

        !number of target state combinations (triangular form
        num_ci = num_targets * (num_targets + 1) / 2

        ii = 1
        do n1 = 1, this % options % num_target_state_sym(i1)
            do n2 = n1, this % options % num_target_state_sym(i1)
                !Construct a master prototype for each target state combo
                !And get the starting matrix index
                call this % compute_matrix_ij(i1, n1, 2, i1, n2, 1, matrix_index(1,ii), matrix_index(2,ii))
                ii = ii + 1
            end do
        end do

        call master_prototype % construct(num_targets * (num_targets + 1) / 2, 1)

        !Construct the temporary prototype
        call temp_prototype % construct

        !get the starting index of the configurations
        starting_index = this % get_starting_index(i1)

        !Get the number of configurations
        num_configurations = this % options % num_ci_target_sym(i1)

        !Set up our mpi variables
        loop_skip = max(1, grid % lprocs)
        my_idx = max(grid % lrank, 0)
        call master_timer % start_timer("Class 3 Prototype")

        !Do the diagonal
        do loop_ido = 1, num_configurations, loop_skip
            !Calculate our real index
            ido = loop_ido + my_idx
            !If we've gon past (which is possible in this setup then skip
            if (ido > num_configurations) cycle
            !Calculate our diagonal configuration index
            config_idx_a = starting_index + (ido - 1) * this % csf_skip
            !The second configuration is simply offset by one
            config_idx_b = config_idx_a + 1
            !Slater it to get our temporary prototypes
            call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 0, .false.)

            alpha_coeffs = 0.0_wp
            !Similar contraction with class 1 Left hand of Eq. 6
            ii = 1
            do n1 = 1, this % options % num_target_state_sym(i1)
                do n2 = n1, this % options % num_target_state_sym(i1)
                    alpha_coeffs(ii,1) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i1, ido, n2)
                    ii = ii + 1
                end do
            end do
            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)
            !call this%fast_contract_class_1_3_diag(i1,num_targets,ido,temp_prototype,master_prototype)
            !Clear temporary prototype
            call temp_prototype % clear

            !Now loop through offdiagonal
            do jdo = 1, num_configurations
                !Skip diagonal (already done above)
                if (jdo == ido) cycle
                !Get the second CSF index offset by one
                config_idx_b = starting_index + 1 + (jdo - 1) * this % csf_skip
                !Slater
                call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 1, .false.)

                !No prototypes? then we leave
                if (temp_prototype % is_empty()) cycle
                alpha_coeffs = 0.0_wp
                !Contract similary to class 1
                ! (cimn*cim'n' + cim'n*cimn')H_imj,imj' Right hand of Eq. 6
                ii = 1
                do n1 = 1, this % options % num_target_state_sym(i1)
                    do n2 = n1, this % options % num_target_state_sym(i1)
                        alpha_coeffs(ii,1) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i1, jdo, n2)
                        ii = ii + 1
                    end do
                end do
                call master_prototype % add_symbols(temp_prototype, alpha_coeffs)
                !call this%fast_contract_class_1_3_offdiag(i1,num_targets,ido,jdo,temp_prototype,master_prototype)
                call temp_prototype % clear
            end do
        end do
        call master_timer % stop_timer("Class 3 Prototype")
        !call master_timer % report_timers
        call master_prototype % synchronize_symbols()

        !Again without expansion, evaluate the first elements
        !ii=1
        !do n1=1,num_targets
        !    do n2=n1,num_targets
                ! cimn*cimn'

        !            call master_prototype(ii)%synchronize_symbols()

        !            if(myrank == master) then
        !                mat_coeffs(1) = this%evaluate_integrals(master_prototype(ii),NORMAL_PHASE)
        !
        !                call matrix_elements%insert_matrix_element(matrix_index(1,ii),matrix_index(2,ii),mat_coeffs(1),8)
        !                write(2014,*) matrix_index(1,ii),matrix_index(2,ii),mat_coeffs(1)
        !            endif
        !            call matrix_elements%update_pure_L2(.false.)
                    !this time it is possible that if the target states are not the same then we 'flip' the index over and insert it
        !            if (n1 /= n2) then
        !                if(myrank == master) then
        !                 call matrix_elements%insert_matrix_element(matrix_index(1,ii)-1,matrix_index(2,ii)+1,mat_coeffs(1),8)
        !                 write(2014,*) matrix_index(1,ii)-1,matrix_index(2,ii)+1,mat_coeffs(1)
        !                endif
        !                 call matrix_elements%update_pure_L2(.false.)
        !            endif
                    !Calaculate the next index
        !            matrix_index(1,ii) = matrix_index(1,ii) - 1

        !            ii=ii+1

        !    enddo
        !enddo

        !Now get expansion variable variables
        config_idx_a = starting_index + (num_configurations - 1) * this % csf_skip
        config_idx_b = config_idx_a + 1

        continuum_orbital_a = this % csfs(config_idx_a) % orbital_sequence_number
        continuum_orbital_b = this % csfs(config_idx_b) % orbital_sequence_number
        call master_timer % start_timer("Class 3 Expansion")

        !Set up our mpi variables
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        total_orbs = compute_total_box(num_continuum_orbitals, num_continuum_orbitals)

        do loop_ido = 1, total_orbs, loop_skip

            loop_idx = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false., num_targets * (num_targets - 1))
            !Skip if beyond total
            if (loop_idx > total_orbs) cycle
            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx, num_continuum_orbitals, ido, jdo)

            if (ido > num_continuum_orbitals) cycle
            if (jdo > num_continuum_orbitals) cycle
            if (jdo < ido + 1) cycle
            if (ido == jdo) cycle

            !Loop through each continuum orbital
            !do ido=1,num_continuum_orbitals-1
            !do jdo=max(3,ido+1),num_continuum_orbitals
            ii = 1
            !Loop through each target state
            !Expand and Evaluate the matrix element

            call this % contr_expand_off_diagonal_eval(continuum_orbital_a, continuum_orbital_b, &
                                                       master_prototype, ido, jdo, mat_coeffs)

            do n1 = 1, num_targets
                do n2 = n1, num_targets
                    !Insert into the matrix
                    call matrix_elements % insert_matrix_element(matrix_index(2,ii) + jdo - 1, &
                                                                 matrix_index(1,ii) + ido - 2, mat_coeffs(ii,1), 8)
                    !print *,matrix_index(2,ii)+jdo-2,matrix_index(1,ii)+ido-1
                    !write(2014,*) matrix_index(1,ii)+jdo-2,matrix_index(2,ii)+ido-1,mat_coeffs(ii)
                    !call matrix_elements%update_continuum(.false.)
                    !this time it is possible that if the target states are not the same then we 'flip' the index over and insert it
                    if (n1 /= n2) then
                        call matrix_elements % insert_matrix_element (matrix_index(2,ii) + ido - 1, &
                                                                      matrix_index(1,ii) + jdo - 2, mat_coeffs(ii,1), 8)
                    !write(2014,*) matrix_index(1,ii)+ido-2,matrix_index(2,ii)+jdo-1,mat_coeffs(ii)
                    !call matrix_elements%update_continuum(.false.)
                    end if
                    ii = ii + 1
                end do
            end do
            !enddo
        end do

        call master_timer % stop_timer("Class 3 Expansion")
        !call master_timer%report_timers

        !Clean up master
        call master_prototype % destroy

        !Clean up temporary
        call temp_prototype % destroy

        this % non_zero_prototype = this % non_zero_prototype + 1

        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8, size(mat_coeffs))
        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

    end subroutine evaluate_class_3


    !> \brief   Builds the class 4 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to building the target-L2 off-diagonal elements of the hamiltonian.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               The target symmetry
    !> \param[in]  num_states       The number of target states withing the symmetry
    !> \param[in]  l1               The index to the L2 function
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements
    !>
    subroutine evaluate_class_4 (this, i1, num_states, l1, matrix_elements)

        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        integer,           intent(in)    :: i1, num_states, l1

        integer :: config_idx_a             ! CSF index
        integer :: num_continuum_orbitals   ! Number of continuum orbitals
        integer :: num_configurations       ! Number of configurations within the target symmetry
        integer :: n1                       ! Target state index
        integer :: continuum_orbital        ! Continuum orbital index

        integer :: ido, jdo, starting_index, i, j   ! Loop variables
        integer :: loop_skip, loop_idx, my_idx, err ! MPI loop variables
        integer, allocatable :: matrix_index(:)     ! The starting matrix index for each state

        type(ContractedSymbolicElementVector) :: master_prototype   ! The contracted prototypes for each target state
        type(SymbolicElementVector)           :: temp_prototype     ! Resultant slater rule calculation of symbols before contracting for each state
        real(wp), allocatable :: alpha_coeffs(:,:), mat_coeffs(:)   ! Target coefficient and matrix element coefficient

        allocate(matrix_index(num_states), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS4::MATRIX_INDEX')

        allocate(alpha_coeffs(num_states,1), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS4::alpha_coeffs')

        allocate(mat_coeffs(num_states), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS4::mat_coeffs')

        !Get number of continuum orbitals
        num_continuum_orbitals = this % options % num_continuum_orbitals_target(i1)
        !Get number of configuration state functions for this symmtery
        num_configurations = this % options % num_ci_target_sym(i1)

        !get the starting index of the configurations
        starting_index = this % get_starting_index(i1)

        do n1 = 1, num_states
            !Construct and clear a symbol list for each state
            !Compute starting index for each state
            call this % compute_matrix_ij(i1, n1, 1, 1, 1, 1, matrix_index(n1), j)
        end do

        call master_prototype % construct(num_states, 1)

        !Construct our temporary symbols
        call temp_prototype % construct
        !Clear them just in case
        call temp_prototype % clear

        !Setup our loop variables
        !loop_skip = max(1,nprocs)
        !my_idx = max(myrank,0)

        !Again like class 1 and 3 we loop in a strided fashion across mpi nodes (or in a normal fashion without mpi)
        do ido = 1, num_configurations
            !Calculate our real index
            !loop_idx = ido + my_idx
            !If we've gone past it then leave
            !if(loop_idx > num_configurations) cycle
            !Caculate our CSf index
            config_idx_a = starting_index + (ido - 1) * this % csf_skip

            !Slater it with the L2 functions
            call this % slater_rules(this % csfs(l1), this % csfs(config_idx_a), temp_prototype, 1, .false.)
            !Nothing? Skip
            if (temp_prototype % is_empty()) cycle

            !Contract of the form cimn*Himj,l (Eq. 7)
            do n1 = 1, num_states
                alpha_coeffs(n1,1) = this % get_target_coeff(i1, ido, n1)
            end do

            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)
            call temp_prototype % clear
        end do

        call temp_prototype % clear

        !Get last CSF in symmetry
        config_idx_a = starting_index + (num_configurations - 1) * this % csf_skip
        !Get the associated continuum orbital
        continuum_orbital = this % csfs(config_idx_a) % orbital_sequence_number

        !This is the L2 matrix position
        j = l1 + this % L2_mat_offset

        !do n1=1, num_states

        !        call master_prototype(n1)%synchronize_symbols
        !enddo

        !For each state, evaluate the integrals without expansion and store
        !do n1=1,num_states
        !    mat_coeffs(n1) = this%evaluate_integrals(master_prototype(n1),NORMAL_PHASE)
        !    call matrix_elements%insert_matrix_element(j,matrix_index(n1),mat_coeffs(n1),4)
        !    call matrix_elements%update_continuum(.false.)
        !enddo

        !Start of continuum orbital expansion
        do ido = 1, num_continuum_orbitals
            !ido = loop_idx + my_idx
            !call matrix_elements%update_pure_L2(.false.,num_states)
            !if(ido > num_continuum_orbitals) cycle

            call this % contr_expand_continuum_L2_eval(continuum_orbital, master_prototype, ido, mat_coeffs)

            do n1 = 1, num_states
                !Expand and evaluate
                !Store coefficient into matrix object
                call matrix_elements % insert_matrix_element(j, matrix_index(n1) + ido - 1, mat_coeffs(n1), 8, 0.0_wp)
                !Perform update if needed
                call matrix_elements % update_pure_L2(.false.)
            end do
        end do

        !Clean up prototype symbols
        call master_prototype % destroy
        call temp_prototype % destroy

        this % non_zero_prototype = this % non_zero_prototype + 1

        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8, size(mat_coeffs))

        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

    end subroutine evaluate_class_4


    !> \brief Builds the class 5 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to building the 'diagonal' off-diagonal elements of the hamiltonian of differing target symmetries
    !> but same continua symmetry.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               Target symmetry 1
    !> \param[in]  num_target_1     The number of target states within target symmetry 1
    !> \param[in]  i2               Target symmetry 2
    !> \param[in]  num_target_2     The number of target states within target symmetry 2
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements
    !>
    subroutine evaluate_class_5 (this, i1, num_target_1, i2, num_target_2, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        integer, intent(in)              :: i1, i2, num_target_1, num_target_2

        !>Contracted prototypes for each target state combination
        !type(SymbolicElementVector),target        ::    master_prototype(num_target_1,num_target_2)
        !type(SymbolicElementVector),pointer        ::    mst_pntr(:)
        type(ContractedSymbolicElementVector) :: master_prototype
        !>Starting index for each
        integer,  allocatable :: matrix_index(:,:,:)
        real(wp), allocatable :: alpha_coeffs(:,:)
        !>Temporary storage of Slater rule calculations
        type(SymbolicElementVector) :: temp_prototype
        !>Configuration and target state variables
        integer                    ::    config_idx_a,config_idx_b,n1,n2
        !> Number of continuum orbitals
        integer                    ::    num_continuum_orbitals
        !>Counts for number of configurations of target symmetry 1,2 and in total TS1*TS2
        integer                    ::    num_configurations_a,num_configurations_b,total_configurations
        !------------------------------MPI LOOP VARIABLES----------------------------------!
        integer                    ::    loop_idx,my_idx,loop_skip,loop_ido
        !---------------------------------LOOP VARIABLES--------------------------------!
        integer                    ::    ido,jdo,starting_index_a,starting_index_b,i,j,err
        !> Current continuum orbital
        integer                    ::    continuum_orbital
        !

        !> Target coefficient and matrix coefficient
        real(wp)               :: alpha
        real(wp), allocatable  :: mat_coeffs(:,:)

        allocate(matrix_index(2, num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS5::MATRIX_INDEX')

        allocate(alpha_coeffs(num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS5::alpha_coeffs')

        allocate(mat_coeffs(num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS5::mat_coeffs')

        !Get number of continuum orbitals, which is a minimum of the two
        num_continuum_orbitals = min(this % options % num_continuum_orbitals_target(i1), &
                                     this % options % num_continuum_orbitals_target(i2))

        !get the starting index of the configuration of target symmetry 1
        starting_index_a = this % get_starting_index(i1)

        !get the starting index of the configurations of target symmetry 2
        starting_index_b = this % get_starting_index(i2)

        !total number of configurations of target symmetry 1
        num_configurations_a = this % options % num_ci_target_sym(i1)

        !total number of configurations of target symmetry 2
        num_configurations_b = this % options % num_ci_target_sym(i2)

        !Compute starting matrix index and construct contracted symbolic prototypes for each combination of target states
        do n1 = 1, num_target_1
            do n2 = 1, num_target_2
                call this % compute_matrix_ij(i1, n1, 1, i2, n2, 1, matrix_index(1,n1,n2), matrix_index(2,n1,n2))
            end do
        end do
        call master_prototype % construct(num_target_1, num_target_2)

        !Construct
        call temp_prototype % construct
        call temp_prototype % clear

        !mst_pntr(1:num_target_1*num_target_2) => master_prototype(:,:)
        !mat_ptr(1:num_target_1*num_target_2) => mat_coeffs(:,:)
        !Compute the total number of configurations to loop through (num_config_a*num_config_b)
        total_configurations = compute_total_box(num_configurations_a, num_configurations_b)

        !Set up our mpi variables
        loop_skip = max(1, grid % lprocs)
        my_idx = max(grid % lrank, 0)

        call master_timer % start_timer("Class 5 prototype")
        !This loop differs from classes 1,3 and 4. Instead the 2 dimensional loop is collapsed into one dimension
        !And the expected CSFs are computed from the single index. This was done to improve load balancing significantly
        !and ensure better OpenMP looping when it is added eventually. If this wasn;t done then the last MPI process will have the largest
        !chunk of configurations to deal compared to the others.
        !Again this loop is strided by MPI processes
        do loop_ido = 1, total_configurations, loop_skip
            alpha_coeffs = 0.0_wp

            !Calculate the real loop index
            loop_idx = loop_ido + my_idx
            !Skip if beyond total
            if (loop_idx > total_configurations) cycle
            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx, num_configurations_a, ido, jdo)

            alpha_coeffs = 0.0_wp

            !Use the two dimensional id to calculate the CSF index
            config_idx_a = starting_index_a + (ido - 1) * this % csf_skip
            config_idx_b = starting_index_b + (jdo - 1) * this % csf_skip

            !Slater for prototype symbols
            call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 1, .false.)

            if (temp_prototype % is_empty()) cycle
            !Perform contraction on each combination of target state of the form cimn*ci'm'n'*Himj,i'm'j' seen in Eq. 8
            !do n1=1,num_target_1
            !    do n2=1,num_target_2
            !        alpha = this%get_target_coeff(i1,ido,n1)*this%get_target_coeff(i2,jdo,n2)
            !        call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
            !    enddo
            !enddo
            !call this%fast_contract_class_567(i1,i2,num_target_1,num_target_2,ido,jdo,temp_prototype,master_prototype)
            do n1 = 1, num_target_1
                do n2 = 1, num_target_2
                    alpha_coeffs(n1,n2) = this % get_target_coeff(i1,ido,n1) * this % get_target_coeff(i2,jdo,n2)
                    !call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
                end do
            end do
            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)

            !clear the temporary symbols
            call temp_prototype % clear
        end do

        call master_timer % stop_timer("Class 5 prototype")
        call temp_prototype % destroy

        !call master_memory%print_memory_report()
        !write(stdout,"('Class 5 num symbols ',4i16)") master_prototype%get_size(),num_target_1,num_target_2,master_prototype%get_size()*num_target_1*num_target_2
        !Get the last CSf and find the associated orbital
        config_idx_a = starting_index_a + (num_configurations_a - 1) * this % csf_skip
        continuum_orbital = this % csfs(config_idx_a) % orbital_sequence_number

        !Evaluate integrals for the first element without expansion. store and possible update the matrix.
        !do n1=1,num_target_1
        !    do n2=1,num_target_2
        !        mat_coeffs(n1,n2) = this%evaluate_integrals(master_prototype(n1,n2),NORMAL_PHASE)
        !        call matrix_elements%insert_matrix_element(matrix_index(2,n1,n2),matrix_index(1,n1,n2),mat_coeff,5)
        !        call matrix_elements%update_continuum(.false.)
        !    enddo
        !enddo

        call master_timer % start_timer("Class 5 Synchonize")
        !write(stdout,"('Synch-Start ')")
                !call master_prototype(n1,n2)%print
        call master_prototype % synchronize_symbols
        !write(stdout,"('Synch-End ')")
        call master_timer % stop_timer("Class 5 Synchonize")

        !Set up our mpi variables
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        !write(stdout,"('Expand ')")
        do loop_ido = 1, num_continuum_orbitals, loop_skip
            !Calculate the real loop index
            ido = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false., num_target_1 * num_target_2)
            !Skip if beyond total
            if (ido > num_continuum_orbitals) cycle
            mat_coeffs = 0
            call this % contr_expand_diagonal_eval(continuum_orbital, master_prototype, ido, mat_coeffs)
            !write(stdout,"('Orb ',i8)") ido
            do n1 = 1, num_target_1
                do n2 = 1, num_target_2
                    !Expand and evaluate
                    !Insert into the matrix
                    call matrix_elements % insert_matrix_element(matrix_index(2,n1,n2) + ido - 1, &
                                            matrix_index(1,n1,n2) + ido - 1, mat_coeffs(n1,n2), 8)
                    !Update matrix if needed.
                end do
            end do
        end do

        !---Cleanup----!

        call master_prototype % destroy

        this % non_zero_prototype = this % non_zero_prototype + 1

        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8,   size(mat_coeffs))

        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

    end subroutine evaluate_class_5


    !> \brief   Builds the class 6 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to building the off-diagonal elements of the hamiltonian of differing target symmetries
    !> but same continua symmetry.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               Target symmetry 1
    !> \param[in]  num_target_1     The number of target states within target symmetry 1
    !> \param[in]  i2               Target symmetry 2
    !> \param[in]  num_target_2     The number of target states within target symmetry 2
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine evaluate_class_6 (this, i1, num_target_1, i2, num_target_2, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        integer,           intent(in)    :: i1, i2, num_target_1, num_target_2
        class(BaseMatrix), intent(inout) ::  matrix_elements

        !>Contracted prototypes for each target state combination
       !type(SymbolicElementVector),target     :: master_prototype(num_target_1,num_target_2)
       !type(SymbolicElementVector),pointer    :: mst_pntr(:)
        type(ContractedSymbolicElementVector)  :: master_prototype

        !>Starting index for each
        integer,  allocatable :: matrix_index(:,:,:)
        real(wp), allocatable :: alpha_coeffs(:,:)

        !>Temporary storage of Slater rule calculations
        type(SymbolicElementVector) :: temp_prototype

        !>Configuration and target state variables
        integer :: config_idx_a, config_idx_b, n1, n2

        !> Number of continuum orbitals
        integer :: num_continuum_orbitals_a, num_continuum_orbitals_b

        !>Counts for number of configurations of target symmetry 1,2 and in total TS1*TS2
        integer :: num_configurations_a, num_configurations_b, total_configurations

        !------------------------------MPI LOOP VARIABLES----------------------------------!
        integer :: loop_idx, my_idx, loop_skip, loop_ido

        !---------------------------------LOOP VARIABLES--------------------------------!
        integer :: ido, jdo, starting_index_a, starting_index_b, i, j, err

        !> Current continuum orbital for target symmetry 1 and 2
        integer :: continuum_orbital_a, continuum_orbital_b, total_orbs

        !> Target coefficient and matrix coefficient
        real(wp)              :: alpha
        real(wp), allocatable :: mat_coeffs(:,:)

        allocate(matrix_index(2, num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS6::MATRIX_INDEX')

        allocate(alpha_coeffs(num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS6::alpha_coeffs')

        allocate(mat_coeffs(num_target_1, num_target_2), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS6::mat_coeffs')

        !Get number of continuum orbitals for target symmetry 1
        num_continuum_orbitals_a = this % options % num_continuum_orbitals_target(i1)

        !Get number of continuum orbitals for target symmetry 2
        num_continuum_orbitals_b = this % options % num_continuum_orbitals_target(i2)

        !get the starting index of the configurations for target symmetry 1
        starting_index_a = this % get_starting_index(i1)

        !get the starting index of the configurations for target symmetry 2
        starting_index_b = this % get_starting_index(i2)

        !get the number of CSFS for target symmetry 1
        num_configurations_a = this % options % num_ci_target_sym(i1)

        !get the number of CSFS for target symmetry 2
        num_configurations_b = this % options % num_ci_target_sym(i2)

        !Compute starting matrix index and construct contracted symbolic prototypes for each combination of target states
        do n1 = 1, num_target_1
            do n2 = 1, num_target_2
                call this % compute_matrix_ij(i1, n1, 1, i2, n2, 1, matrix_index(1,n1,n2), matrix_index(2,n1,n2))
            end do
        end do
        call master_prototype % construct(num_target_1, num_target_2)

        !Construct
        call temp_prototype % construct
        call temp_prototype % clear
        !mst_pntr(1:num_target_1*num_target_2) => master_prototype(:,:)
        !mat_ptr(1:num_target_1*num_target_2) => mat_coeffs(:,:)
        !Compute the total number of configurations to loop through (num_config_a*num_config_b)
        total_configurations = compute_total_box(num_configurations_a, num_configurations_b)

        !Set up our mpi variables
        loop_skip = max(1, grid % lprocs)
        my_idx = max(grid % lrank, 0)

        config_idx_a = 0
        config_idx_b = 0

        call master_timer % start_timer("Class 6 prototype")
        !This loop differs from classes 1,3 and 4. Instead the 2 dimensional loop is collapsed into one dimension
        !And the expected CSFs are computed from the single index. This was done to improve load balancing significantly
        !and ensure better OpenMP looping when it is added eventually. If this wasn't done then the last MPI process will have the largest
        !chunk of configurations to deal compared to the others.
        !Again this loop is strided by MPI processes
        do loop_ido = 1, total_configurations, loop_skip
            alpha_coeffs = 0.0_wp

            !Calculate the real loop index
            loop_idx = loop_ido + my_idx
            !Skip if beyond total
            if (loop_idx > total_configurations) cycle
            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx, num_configurations_a, ido, jdo)

            !Use the two dimensional id to calculate the CSF index
            config_idx_a = starting_index_a + (ido - 1) * this % csf_skip
            config_idx_b = starting_index_b + 1 + (jdo - 1) * this % csf_skip

            !Slater for prototype symbols
            call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 1, .false.)

            if (temp_prototype % is_empty()) cycle
            !Perform contraction on each combination of target state of the form cimn*ci'm'n'*Himj,i'm'j' seen in Eq. 8
            !do n1=1,num_target_1
            !    do n2=1,num_target_2
            !        alpha = this%get_target_coeff(i1,ido,n1)*this%get_target_coeff(i2,jdo,n2)
            !        call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
            !    enddo
            !enddo
            !call this%fast_contract_class_567(i1,i2,num_target_1,num_target_2,ido,jdo,temp_prototype,master_prototype)
            do n1 = 1, num_target_1
                do n2 = 1, num_target_2
                    alpha_coeffs(n1,n2) = this % get_target_coeff(i1,ido,n1) * this % get_target_coeff(i2,jdo,n2)
                    !call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
                end do
            end do
            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)

            !clear the temporary symbols
            call temp_prototype % clear
        end do

        call master_timer % stop_timer("Class 6 prototype")
        call temp_prototype % clear

        ! When there are more ranks in the shared memory communicator than pairs of target CSFs,
        ! some ranks exit the previous DO cycle without setting config_idx_a and config_idx_b. These are needed
        ! just to get number of
        n1 = config_idx_a; call mpi_reduceall_max(n1, config_idx_a, grid % lcomm)
        n2 = config_idx_b; call mpi_reduceall_max(n2, config_idx_b, grid % lcomm)

        !Get the continuum orbital numbers for both of the last CSFS
        continuum_orbital_a = this % csfs(config_idx_a) % orbital_sequence_number
        continuum_orbital_b = this % csfs(config_idx_b) % orbital_sequence_number

        call master_timer % start_timer("Class 6 Synchonize")

        !call master_prototype(n1,n2)%print
        call master_prototype % synchronize_symbols
        call master_timer % stop_timer("Class 6 Synchonize")

        total_orbs = compute_total_box(num_continuum_orbitals_a, num_continuum_orbitals_b)

        !Set up our mpi variables
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        !Loop through each orbital in target symmetry 1
        do loop_ido = 1, total_orbs, loop_skip
            !Loop through each orbital in target symmetry 2
            !do jdo=1,num_continuum_orbitals_b
            loop_idx = loop_ido + my_idx
            call matrix_elements%update_pure_L2(.false.,num_target_1*num_target_2)

            !Skip if beyond total
            if(loop_idx > total_orbs) cycle

            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx,num_continuum_orbitals_a,ido,jdo)
            if (ido == jdo) cycle !Class 5 already handled the diagonal so ignore

            !Loop through each state in target symmetry 1
            !Expand and evaluate immediately
            call this % contr_expand_off_diagonal_eval (continuum_orbital_a, continuum_orbital_b, &
                                                        master_prototype, ido, jdo, mat_coeffs)
            do n1 = 1, num_target_1
                !Loop through each state in target symmetry 2
                do n2 = 1, num_target_2
                    !Insert into the matrix
                    call matrix_elements % insert_matrix_element (matrix_index(2,n1,n2) + jdo - 1, &
                                                                  matrix_index(1,n1,n2) + ido - 1, &
                                                                  mat_coeffs(n1,n2), 8)
                    !Update if neccessary
                end do
            end do
        end do

        !Cleanup

        call master_prototype % destroy
        call temp_prototype % destroy

        this % non_zero_prototype = this % non_zero_prototype + 1

        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8,   size(mat_coeffs))

        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

    end subroutine evaluate_class_6


    !> \brief   Builds the class 7 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This relates to building the off-diagonal elements of the hamiltonian of differing target symmetries and continua symmetry
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[in]  i1               Target symmetry 1
    !> \param[in]  num_target_1     The number of target states within target symmetry 1
    !> \param[in]  i2               Target symmetry 2
    !> \param[in]  num_target_2     The number of target states within target symmetry 2
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements
    !>
    subroutine evaluate_class_7 (this, i1, num_target_1, i2, num_target_2, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        integer,           intent(in)    :: i1, i2, num_target_1, num_target_2
        class(BaseMatrix), intent(inout) :: matrix_elements

        ! Contracted prototypes for each target state combination
       !type(SymbolicElementVector), target   :: master_prototype(num_target_1, num_target_2)
       !type(SymbolicElementVector), pointer  :: mst_pntr(:)
        type(ContractedSymbolicElementVector) :: master_prototype

        ! Starting index for each
        integer,  allocatable :: matrix_index(:,:,:)
        real(wp), allocatable :: alpha_coeffs(:,:)

        ! Temporary storage of Slater rule calculations
        type(SymbolicElementVector) :: temp_prototype

        ! Configuration and target state variables
        integer :: config_idx_a, config_idx_b, n1, n2

        ! Number of continuum orbitals
        integer :: num_continuum_orbitals_a, num_continuum_orbitals_b

        ! Counts for number of configurations of target symmetry 1,2 and in total TS1*TS2
        integer :: num_configurations_a, num_configurations_b, total_configurations

        !------------------------------MPI LOOP VARIABLES----------------------------------!
        integer :: loop_idx, my_idx, loop_skip, loop_ido

        !---------------------------------LOOP VARIABLES--------------------------------!
        integer :: ido, jdo, starting_index_a, starting_index_b, i, j, err

        ! Current continuum orbital for target symmetry 1 and 2
        integer :: continuum_orbital_a, continuum_orbital_b, total_orbs

        ! Target coefficient and matrix coefficient
        real(wp)              :: alpha
        real(wp), allocatable :: mat_coeffs(:,:)

        allocate(matrix_index(2,num_target_1,num_target_2), stat = err)
        call master_memory % track_memory(storage_size(matrix_index)/8, size(matrix_index), err, 'CLASS7::MATRIX_INDEX')

        allocate(alpha_coeffs(num_target_1,num_target_2), stat = err)
        call master_memory % track_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs), err, 'CLASS7::alpha_coeffs')

        allocate(mat_coeffs(num_target_1,num_target_2), stat = err)
        call master_memory % track_memory(storage_size(mat_coeffs)/8, size(mat_coeffs), err, 'CLASS7::mat_coeffs')

        config_idx_a = 0
        config_idx_b = 0

        !Get number of continuum orbitals for target symmetry 1
        num_continuum_orbitals_a = this % options % num_continuum_orbitals_target(i1)
        !Get number of continuum orbitals for target symmetry 2
        num_continuum_orbitals_b = this % options % num_continuum_orbitals_target(i2)

        !get the starting index of the configurations for target symmetry 1
        starting_index_a = this % get_starting_index(i1)
        !get the starting index of the configurations for target symmetry 2
        starting_index_b = this % get_starting_index(i2)

        !get the number of CSFS for target symmetry 1
        num_configurations_a = this % options % num_ci_target_sym(i1)
        !get the number of CSFS for target symmetry 2
        num_configurations_b = this % options % num_ci_target_sym(i2)

        !Compute starting matrix index and construct contracted symbolic prototypes for each combination of target states
        do n1 = 1, num_target_1
            do n2= 1, num_target_2
                call this % compute_matrix_ij(i1, n1, 1, i2, n2, 1, matrix_index(1,n1,n2), matrix_index(2,n1,n2))
            end do
        end do

        call master_prototype % construct(num_target_1, num_target_2)

        !Construct
        call temp_prototype % construct
        call temp_prototype % clear
        !mst_pntr(1:num_target_1*num_target_2) => master_prototype(:,:)
        !mat_ptr(1:num_target_1*num_target_2) => mat_coeffs(:,:)
        !Compute the total number of configurations to loop through (num_config_a*num_config_b)
        total_configurations = compute_total_box(num_configurations_a, num_configurations_b)

        !Set up our mpi variables
        loop_skip = max(1, grid % lprocs)
        my_idx = max(grid % lrank, 0)

        call master_timer % start_timer("Class 7 prototype")
        !This loop differs from classes 1,3 and 4. Instead the 2 dimensional loop is collapsed into one dimension
        !And the expected CSFs are computed from the single index. This was done to improve load balancing significantly
        !and ensure better OpenMP looping when it is added eventually. If this wasn't done then the last MPI process will have the largest
        !chunk of configurations to deal compared to the others.
        !Again this loop is strided by MPI processes
        do loop_ido = 1, total_configurations, loop_skip
            alpha_coeffs = 0.0_wp

            !Calculate the real loop index
            loop_idx = loop_ido + my_idx
            !Skip if beyond total
            if (loop_idx > total_configurations) cycle
            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx, num_configurations_a, ido, jdo)

            !Use the two dimensional id to calculate the CSF index
            config_idx_a = starting_index_a + (ido - 1) * this % csf_skip
            config_idx_b = starting_index_b + (jdo - 1) * this % csf_skip

            !Slater for prototype symbols
            call this % slater_rules(this % csfs(config_idx_b), this % csfs(config_idx_a), temp_prototype, 1, .false.)

            !Cycle if empty
            if (temp_prototype % is_empty()) cycle
            !Perform contraction on each combination of target state of the form cimn*ci'm'n'*Himj,i'm'j' seen in Eq. 8
            !do n1=1,num_target_1
            !    do n2=1,num_target_2
            !        alpha = this%get_target_coeff(i1,ido,n1)*this%get_target_coeff(i2,jdo,n2)

            !        call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
            !    enddo
            !enddo
            !call this%fast_contract_class_567(i1,i2,num_target_1,num_target_2,ido,jdo,temp_prototype,master_prototype)
            do n1 = 1, num_target_1
                do n2 = 1, num_target_2
                    alpha_coeffs(n1,n2) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i2, jdo, n2)
                    !call master_prototype(n1,n2)%add_symbols(temp_prototype,alpha)
                end do
            end do
            call master_prototype % add_symbols(temp_prototype, alpha_coeffs)

            !clear the temporary symbols
            call temp_prototype % clear
        end do

        call master_timer % stop_timer("Class 7 prototype")
        call temp_prototype % clear

        ! When there are more ranks in the shared memory communicator than pairs of target CSFs,
        ! some ranks exit the previous DO cycle without setting config_idx_a and config_idx_b. These are needed
        ! just to get number of
        n1 = config_idx_a; call mpi_reduceall_max(n1, config_idx_a, grid % lcomm)
        n2 = config_idx_b; call mpi_reduceall_max(n2, config_idx_b, grid % lcomm)

        !Get the continuum orbital numbers for both of the last CSFS
        continuum_orbital_a = this % csfs(config_idx_a) % orbital_sequence_number
        continuum_orbital_b = this % csfs(config_idx_b) % orbital_sequence_number

        call master_timer % start_timer("Class 7 Synchonize")
        call master_prototype % synchronize_symbols
        call master_timer % stop_timer("Class 7 Synchonize")

        total_orbs = compute_total_box(num_continuum_orbitals_a, num_continuum_orbitals_b)

        !Set up our mpi variables
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        call master_timer % start_timer("Class 7 expansion")

        !Loop through each orbital in target symmetry 1
        do loop_ido = 1, total_orbs, loop_skip
            !Loop through each orbital in target symmetry 2
            !do jdo=1,num_continuum_orbitals_b
            loop_idx = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false., num_target_1 * num_target_2)
            !Skip if beyond total
            if (loop_idx > total_orbs) cycle
            !Convert the 1-D index into 2-D
            call box_index_to_ij(loop_idx, num_continuum_orbitals_a, ido, jdo)
            !Loop through each state in target symmetry 1
            !Expand and evaluate immediately
            call this % expand_off_diagonal_gen_eval(continuum_orbital_a, continuum_orbital_b, master_prototype, &
                                                     ido, jdo, mat_coeffs)

            do n1 = 1, num_target_1
                !Loop through each state in target symmetry 2
                do n2 = 1, num_target_2
                    !Insert into the matrix
                    call matrix_elements % insert_matrix_element(matrix_index(2,n1,n2) + jdo - 1, &
                                                                 matrix_index(1,n1,n2) + ido - 1, &
                                                                 mat_coeffs(n1,n2), 8)
                    !Update if neccessary
                end do
            end do
        end do

        call master_timer % stop_timer("Class 7 expansion")

        !-------------Cleanup---------------!
        call temp_prototype % destroy
        call master_prototype % destroy

        this % non_zero_prototype = this % non_zero_prototype + 1

        call master_memory % free_memory(storage_size(matrix_index)/8, size(matrix_index))
        call master_memory % free_memory(storage_size(alpha_coeffs)/8, size(alpha_coeffs))
        call master_memory % free_memory(storage_size(mat_coeffs)/8, size(mat_coeffs))

        deallocate(matrix_index, alpha_coeffs, mat_coeffs)

    end subroutine evaluate_class_7


    !> \brief   Builds the class 8 part of the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \deprecated Now replaced with the fast class_2_and_8 combined subroutine.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements
    !>
    subroutine evaluate_class_8 (this, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements

        type(SymbolicElementVector) :: symbolic_elements  ! The smbolic elements
        integer                     :: l1, l2             ! L2 CSf indicies
        real(wp)                    :: mat_coeff          ! The calculated matrix coefficient
        integer                     :: loop_ido, my_idx, loop_skip

        !Construct the symbolic element object
        call symbolic_elements % construct

        my_idx = max(1, grid % grank)
        loop_skip = max(1, grid % gprocs)

        !Loop through each L2 function
        do loop_ido = this % start_L2, this % total_num_csf - 1, loop_skip
            l1 = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false.)

            !Loop through the lower triangular of L2
            do l2 = this % start_L2, l1 - 1

                !Check for an update here to ensure synchonization with all processes
                if (l1 > this % total_num_csf - 1) exit

                !Slater the L2 functions
                call this % slater_rules(this % csfs(l1), this % csfs(l2), symbolic_elements, 1, .false.)
                !If empty then cycle
                if (symbolic_elements % is_empty()) cycle
                !Evaluate the integral
                mat_coeff = this % evaluate_integrals(symbolic_elements, NORMAL_PHASE)
                !Insert the coffieicnt into the correct matrix position
                call matrix_elements % insert_matrix_element(l2 + this % L2_mat_offset, l1 + this % L2_mat_offset, mat_coeff, 8)
                !Clear for the next step
                call symbolic_elements % clear

                this % non_zero_prototype = this % non_zero_prototype + 1
            end do
        end do

        !Cleanup the symbols
        call symbolic_elements % destroy

    end subroutine evaluate_class_8


    !> \brief   Builds both the class 2 and class 8 (Pure L2) parts of the matrix elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This suborutine builds both the diagonal and off-diagonal elements of the L2 part of the matrix.
    !> It collapses the 2-dimensional loop into a one dimensional loop and then moves through the entire lower
    !> triangular portion o the matrix to calculate very quicky and in parallel (if mpi enabled) the matrix elements.
    !>
    !> \param[inout] this           Hamiltonian object to query.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine evaluate_class_2_and_8 (this, matrix_elements)
        class(Contracted_Hamiltonian)    :: this
        class(BaseMatrix), intent(inout) :: matrix_elements

        integer                     :: l1, l2               ! L2 indicies to be calculated
        type(SymbolicElementVector) :: symbolic_elements    ! The symbolic elements evaluated from Slaters rules
        real(wp)                    :: mat_coeff            ! The calculated matrix coefficients
        integer                     :: loop_idx, ido        ! Looping indicies
        integer                     :: total_vals           ! Total number of CSFS to loop through
        integer                     :: loop_skip, my_idx    ! MPI variables
        integer                     :: trig_n               ! The total number of L2 CSFS
        integer                     :: diag                 ! Wheter we are diagonal or not

        call symbolic_elements % construct              ! Construct our symbol object
        trig_n = this % num_L2_functions                ! Calculate the number of L2 functions
        total_vals = compute_total_box(trig_n, trig_n)  ! Compute the total number of CSFS to loop through

        this % diagonal_flag = 0

        !Our MPI looping variasbles
        loop_skip = max(1, grid % gprocs)
        my_idx= max(grid % grank, 0)

        write (stdout, "('start,total: ',2i12)") this % start_L2, total_vals

        ! Whilst similar in concept to classes 5,6,7. Instead of collapsing an N^2 to one dimension it collapses an
        ! N(N+1)/2 loop into one dimension. The index calculation is a little mroe complex but offers the same performance
        ! benefits as classes 5,6 and 7
        do ido = 1, total_vals, loop_skip

                !Update if needed
                call matrix_elements % update_pure_L2(.false.)
                call symbolic_elements % clear

                !Calculate the real index
                loop_idx = ido + my_idx
                !If we've gone over then cycle
                if (loop_idx > total_vals) cycle
                !Convert the one dimensional index to two dimensions
                call box_index_to_ij(loop_idx, trig_n, l1, l2)
               !write (stdout, "('idx = ',i12,' l1 = ',i12,' l2 ',i12)") loop_idx, l1, l2
                !Get the correct L2 offsets
                l1 = l1 + this % start_L2 - 1
                l2 = l2 + this % start_L2 - 1
               !write (stdout, "('idx = ',i12,' d l1 = ',i12,'d l2 ',i12,'mat ',2i8)") loop_idx, l1, l2, l2 + this % L2_mat_offset, l1 + this % L2_mat_offset

                if ((l2 + this % L2_mat_offset) < (l1 + this % L2_mat_offset)) cycle

                !As a precaution, if we've gone past either, skip
                if (l1 > this % total_num_csf .or. l2 > this % total_num_csf) cycle
                if (l1 < this % start_L2 .or. l2 < this % start_L2 ) cycle
                !write(stdout,"('l1,2 = ',2i8)") l1,l2
                !stop "Error in L2 indexing"
                !endif

                diag = 1
                if (l1 == l2) diag = 0

                !Slater the CSFS
                call this % slater_rules(this % csfs(l1), this % csfs(l2), symbolic_elements, diag, .false.)
                !If empty, skip
                if (symbolic_elements % is_empty()) cycle
                !Evaluate normally
                mat_coeff = this % evaluate_integrals(symbolic_elements, NORMAL_PHASE)

                !if(l2/=l1 .and. abs(mat_coeff) >= DEFAULT_MATELEM_THRESHOLD) write(1235,"(3i8,es16.5)") diag,l2+this%L2_mat_offset,l1+this%L2_mat_offset,mat_coeff
                !if(l1==l2) call this%rmat_ci(1)%modify_L2_diagonal(mat_coeff)

                !Insert the matrix element into the correct place
                call matrix_elements % insert_matrix_element(l2 + this % L2_mat_offset, l1 + this % L2_mat_offset, mat_coeff, 8)

                !Clear for the next cycle
                this % non_zero_prototype = this % non_zero_prototype + 1

        end do

        !Cleanup
        call symbolic_elements % destroy

    end subroutine evaluate_class_2_and_8


    !> \brief   Finds the starting matrix index for paticular target symmetries and target states
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this  Hamiltonian object to query.
    !> \param[in] i1 Target symmetry 1
    !> \param[in] n1 Target state of target symmetry 1
    !> \param[in] j1 Continuum orbital of Target state of target symmetry 1
    !> \param[in] i2 Target symmetry 2
    !> \param[in] n2 Target state of target symmetry 2
    !> \param[in] j2 Continuum orbital of Target state of target symmetry 2
    !> \param[out] i Resultant matrix co-ordinate i
    !> \param[out] j Resultant matrix co-ordinate j
    !>
    subroutine compute_matrix_ij (this, i1, n1, j1, i2, n2, j2, i, j)
        class(Contracted_Hamiltonian) :: this
        integer, intent(in)           :: i1, n1, j1, i2, n2, j2
        integer, intent(out)          :: i, j
        integer                       :: ido

        i = 0

        do ido = 2, i1
            i = i + this % options % num_target_state_sym(ido - 1) &
                  * this % options % num_continuum_orbitals_target(ido - 1)
        end do

        do ido = 2, n1
            i = i + this % options % num_continuum_orbitals_target(i1)
        end do

        i = i + j1
        j = 0

        do ido = 2, i2
            j = j + this % options % num_target_state_sym(ido - 1) &
                  * this % options % num_continuum_orbitals_target(ido - 1)
        enddo

        do ido = 2, n2
            j = j + this % options % num_continuum_orbitals_target(i2)
        end do

        j = j + j2

    end subroutine compute_matrix_ij

!----------------------------------EXPANSION-----------------------------------------------------!

    !> Deprecated
    subroutine expand_diagonal (this, continuum_orbital, master_prototype, j1, expanded_prototype)
        class(Contracted_Hamiltonian)            :: this
        class(SymbolicElementVector), intent(in) :: master_prototype, expanded_prototype
        integer,                      intent(in) :: continuum_orbital, j1

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lw1, lw2, ido, i, num_integrals

        num_integrals = master_prototype % get_size()

        do ido = 1, num_integrals

            !Get the integral
            call master_prototype % get_coeff_and_integral(ido, integral_coeff, integral)
            call unpack_ints(integral, lwd)

            lw1 = 0
            lw2 = 0

            do i = 1, 4
                if (lwd(i) == continuum_orbital) then
                    if(lw1 == 0) then
                        lw1 = i
                    else
                        lw2 = i
                    end if
                end if
            end do

            if (lw1 == 0) then
                call expanded_prototype % insert_symbol(integral, integral_coeff, .false.)
                cycle
            else if (lw2 == 0) then
                 lwd(lw1) = continuum_orbital + j1 - 1
            else
                 lwd(lw1) = continuum_orbital + j1 - 1
                 lwd(lw2) = continuum_orbital + j1 - 1
            end if

            call expanded_prototype % insert_ijklm_symbol(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), integral_coeff, .false.)

        end do

    end subroutine expand_diagonal

    !> Deprecated
    subroutine expand_off_diagonal (this, continuum_orbital_a, continuum_orbital_b, master_prototype, ja, jb, expanded_prototype)
        class(Contracted_Hamiltonian)            :: this
        class(SymbolicElementVector), intent(in) :: master_prototype, expanded_prototype
        integer,                      intent(in) :: continuum_orbital_a, continuum_orbital_b, ja, jb

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lwa, lwb, ido, i, num_integrals

        num_integrals = master_prototype % get_size()

        do ido = 1, num_integrals

            !Get the integral
            call master_prototype % get_coeff_and_integral(ido, integral_coeff, integral)
            call unpack_ints(integral, lwd)

            lwa = 0
            lwb = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital_a) lwa = i
               if (lwd(i) == continuum_orbital_b) lwb = i
            end do

            !No occurances
            if (max(lwa, lwb) == 0) then
                call expanded_prototype % insert_symbol(integral, integral_coeff, .false.)
                cycle
            else if (lwb == 0) then
                lwd(lwa) = continuum_orbital_a + ja - 1
            else if (lwa == 0) then
                lwd(lwb) = continuum_orbital_b + jb - 2
            else
                lwd(lwa) = continuum_orbital_a + ja - 1
                lwd(lwb) = continuum_orbital_b + jb - 2
            endif

            !One occurance

            call expanded_prototype % insert_ijklm_symbol(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), integral_coeff, .false.)

        end do

    end subroutine expand_off_diagonal


    !> Deprecated
    subroutine expand_continuum_L2 (this, continuum_orbital, master_prototype, j1, expanded_prototype)
        class(Contracted_Hamiltonian)            :: this
        class(SymbolicElementVector), intent(in) :: master_prototype, expanded_prototype
        integer,                      intent(in) :: continuum_orbital, j1

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lw, ido, i, num_integrals

        num_integrals = master_prototype % get_size()

        do ido = 1, num_integrals

            !Get the integral
            call  master_prototype % get_coeff_and_integral(ido, integral_coeff, integral)
            call  unpack_ints(integral, lwd)

            lw = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital) lw=i
            end do

            if (lw == 0) then
                call expanded_prototype % insert_symbol(integral, integral_coeff, .false.)
                cycle
            else
                 lwd(lw) = continuum_orbital + j1 - 1
            end if

            call expanded_prototype % insert_ijklm_symbol(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), integral_coeff, .false.)

        end do

    end subroutine expand_continuum_L2


    !> Deprecated
    subroutine expand_off_diagonal_general (this, continuum_orbital_a, continuum_orbital_b, master_prototype, &
                                            ja, jb, expanded_prototype)
        class(Contracted_Hamiltonian)            :: this
        class(SymbolicElementVector), intent(in) :: master_prototype, expanded_prototype
        integer,                      intent(in) :: continuum_orbital_a, continuum_orbital_b, ja, jb

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lwa, lwb, ido, i, num_integrals

        num_integrals = master_prototype % get_size()

        do ido = 1, num_integrals

            !Get the integral
            call master_prototype % get_coeff_and_integral(ido, integral_coeff, integral)
            call unpack_ints(integral, lwd)

            lwa = 0
            lwb = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital_a) lwa = i
               if (lwd(i) == continuum_orbital_b) lwb = i
            end do

            !No occurances
            if (max(lwa, lwb) == 0) then
                call expanded_prototype % insert_symbol(integral, integral_coeff, .false.)
                cycle
            else if(lwb == 0)then
                lwd(lwa) = continuum_orbital_a + ja - 1
            else if (lwa == 0) then
                lwd(lwb) = continuum_orbital_b + jb - 1
            else
                lwd(lwa) = continuum_orbital_a + ja - 1
                lwd(lwb) = continuum_orbital_b + jb - 1
            end if

            !IF(ipair(lwd(1))+lwd(2) > ipair(lwd(3))+lwd(4)) then
                call expanded_prototype % insert_ijklm_symbol(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), integral_coeff, .false.)
            !else
                !call expanded_prototype%insert_ijklm_symbol(lwd(3),lwd(4),lwd(1),lwd(2),lwd(5),integral_coeff)
            !endif
        end do

    end subroutine expand_off_diagonal_general


    !> \brief   Expands the diagonal prototype symbols and evaluates for a single  matrix element for co-ordiate (j1)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this            Hamiltonian object to query.
    !> \param[in] continuum_orbital  Prototype continuum orbital
    !> \param[in] master_prototype   Prototype symbolic elements to be expanded
    !> \param[in] j1                 Desired continuum orbital
    !> \param[out] mat_coeffs        Output array for the matrix element.
    !>
    subroutine expand_diagonal_eval (this, continuum_orbital, master_prototype, j1, mat_coeffs)
        class(Contracted_Hamiltonian)            :: this
        integer,                     intent(in)  :: continuum_orbital, j1
        type(SymbolicElementVector), intent(in)  :: master_prototype(:)
        real(wp),                    intent(out) :: mat_coeffs(:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lw1, lw2, ido, i, prototypes, num_prototypes, num_integrals

        num_prototypes = size(master_prototype)
        num_integrals = master_prototype(1) % get_size()
        mat_coeffs = 0

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            call master_prototype(1) % get_coeff_and_integral(ido, integral_coeff, integral)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp

            lw1 = 0
            lw2 = 0

            do i = 1, 4
                if (lwd(i) == continuum_orbital) then
                    if (lw1 == 0) then
                        lw1 = i
                    else
                        lw2 = i
                    end if
                end if
            end do

            if (lw1 == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)
                do prototypes = 1, num_prototypes
                    mat_coeffs(prototypes) = mat_coeffs(prototypes) &
                                           + master_prototype(prototypes) % get_coefficient(ido) * integral_coeff
                end do
                cycle
            else if (lw2 == 0) then
                lwd(lw1) = continuum_orbital + j1 - 1
            else
                lwd(lw1) = continuum_orbital + j1 - 1
                lwd(lw2) = continuum_orbital + j1 - 1
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)
            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do prototypes = 1, num_prototypes
                mat_coeffs(prototypes) = mat_coeffs(prototypes) + master_prototype(prototypes)%get_coefficient(ido) * integral_coeff
            end do

        end do

    end subroutine expand_diagonal_eval


    !> \brief   Expands the diagonal prototype symbols and evaluates for a single  matrix element for co-ordiate (j1)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this            Hamiltonian object to query.
    !> \param[in] continuum_orbital  Prototype continuum orbital
    !> \param[in] master_prototype   Prototype symbolic elements to be expanded
    !> \param[in] j1                 Desired continuum orbital
    !> \param[out] mat_coeffs        Output array for the matrix element.
    !>
    subroutine contr_expand_diagonal_eval (this, continuum_orbital, master_prototype, j1, mat_coeffs)
        class(Contracted_Hamiltonian)                      :: this
        type(ContractedSymbolicElementVector), intent(in)  :: master_prototype
        integer,                               intent(in)  :: continuum_orbital, j1
        real(wp),                              intent(out) :: mat_coeffs(:,:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff

        integer :: lwd(8), lw1, lw2, ido, i, prototypes, num_prototypes, num_integrals, num_targets_1, num_targets_2, n1, n2

        num_integrals = master_prototype % get_size()
        num_targets_1 = master_prototype % get_num_targets_sym1()
        num_targets_2 = master_prototype % get_num_targets_sym2()
        mat_coeffs = 0

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            integral = master_prototype % get_integral_label(ido)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp
            lw1 = 0
            lw2 = 0

            do i = 1, 4
                if (lwd(i) == continuum_orbital) then
                    if (lw1 == 0) then
                        lw1 = i
                    else
                        lw2 = i
                    end if
                end if
            end do

            if (lw1 == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)
                do n1 = 1,num_targets_1
                    do n2 = 1,num_targets_2
                        mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                    end do
                end do
                cycle
            else if (lw2 == 0) then
                lwd(lw1) = continuum_orbital + j1 - 1
            else
                lwd(lw1) = continuum_orbital + j1 - 1
                lwd(lw2) = continuum_orbital + j1 - 1
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)
            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do n1 = 1,num_targets_1
                do n2 = 1,num_targets_2
                    mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                end do
            end do

        end do

    end subroutine contr_expand_diagonal_eval


    !> \brief   Expands the continuum-L2 prototype symbols and evaluates for a single  matrix element for co-ordinate (j1),L2
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this            Hamiltonian object to query.
    !> \param[in] continuum_orbital  Prototype continuum orbital
    !> \param[in] master_prototype   Prototype symbolic elements to be expanded
    !> \param[in] j1                 Desired continuum orbital
    !> \param[out] mat_coeffs        Output array for the matrix element.
    !>
    subroutine expand_continuum_L2_eval (this, continuum_orbital, master_prototype, j1, mat_coeffs)
        class(Contracted_Hamiltonian)             :: this
        class(SymbolicElementVector), intent(in)  :: master_prototype(:)
        integer,                      intent(in)  :: continuum_orbital, j1
        real(wp),                     intent(out) :: mat_coeffs(:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lw, ido, i, prototypes, num_prototypes, num_integrals

        num_prototypes = size(master_prototype)
        num_integrals = master_prototype(1) % get_size()

        mat_coeffs = 0.0_wp

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            call master_prototype(1) % get_coeff_and_integral(ido, integral_coeff, integral)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp
            lw = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital) lw = i
            end do

            if (lw == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)
                do prototypes = 1, num_prototypes
                    mat_coeffs(prototypes) = mat_coeffs(prototypes) &
                                           + master_prototype(prototypes) % get_coefficient(ido) * integral_coeff
                end do
            else
                 lwd(lw) = continuum_orbital + j1 - 1
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)

            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do prototypes = 1, num_prototypes
                mat_coeffs(prototypes) = mat_coeffs(prototypes) + master_prototype(prototypes)%get_coefficient(ido) * integral_coeff
            end do

        end do

    end subroutine expand_continuum_L2_eval


    !> \brief   Expands the continuum-L2 prototype symbols and evaluates for a single matrix element for co-ordinate (j1),L2
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this            Hamiltonian object to query.
    !> \param[in] continuum_orbital  Prototype continuum orbital
    !> \param[in] master_prototype   Prototype symbolic elements to be expanded
    !> \param[in] j1                 Desired continuum orbital
    !> \param[out] mat_coeffs        Output array for the matrix element.
    !>
    subroutine contr_expand_continuum_L2_eval (this, continuum_orbital, master_prototype, j1, mat_coeffs)
        class(Contracted_Hamiltonian)                       :: this
        class(ContractedSymbolicElementVector), intent(in)  :: master_prototype
        integer,                                intent(in)  :: continuum_orbital, j1
        real(wp),                               intent(out) :: mat_coeffs(:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lw, ido, i, num_targets_1, num_targets_2, n1, n2, num_integrals

        num_targets_1 = master_prototype % get_num_targets_sym1()
        num_integrals = master_prototype % get_size()
        mat_coeffs = 0.0_wp

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            integral =  master_prototype % get_integral_label(ido)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp

            lw = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital) lw = i
            end do

            if (lw == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)
                do n1 = 1, num_targets_1
                    mat_coeffs(n1) = mat_coeffs(n1) + master_prototype % get_coefficient(ido, n1, 1) * integral_coeff
                end do
                cycle
            else
                 lwd(lw) = continuum_orbital + j1 - 1
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)

            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do n1 = 1, num_targets_1
                mat_coeffs(n1) = mat_coeffs(n1) + master_prototype % get_coefficient(ido, n1, 1) * integral_coeff
            end do

        end do

    end subroutine contr_expand_continuum_L2_eval


    !> \brief   Expands the off-diagonal same/similar symmetry prototype symbols and evaluates for a single matrix element for co-ordinate offset (ja,jb)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this              Hamiltonian object to query.
    !> \param[in] continuum_orbital_a  Prototype continuum orbital for target symmetry a
    !> \param[in] continuum_orbital_b  Prototype continuum orbital for target symmetry b
    !> \param[in] master_prototype     Prototype symbolic elements to be expanded
    !> \param[in] ja                   Desired continuum orbital of target symmetry a
    !> \param[in] jb                   Desired continuum orbital of target symmetry b
    !> \param[in] mat_coeffs           Output array for the matrix element.
    !>
    subroutine expand_off_diagonal_eval (this, continuum_orbital_a, continuum_orbital_b, master_prototype, ja, jb, mat_coeffs)
        class(Contracted_Hamiltonian)             :: this
        class(SymbolicElementVector), intent(in)  :: master_prototype(:)
        integer,                      intent(in)  :: continuum_orbital_a, continuum_orbital_b, ja, jb
        real(wp),                     intent(out) :: mat_coeffs(:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lwa, lwb, ido, i, prototypes, num_prototypes, num_integrals

        num_prototypes = size(master_prototype)
        num_integrals = master_prototype(1) % get_size()
        mat_coeffs = 0

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            call  master_prototype(1) % get_coeff_and_integral(ido, integral_coeff, integral)
            call  unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp

            lwa = 0
            lwb = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital_a) lwa = i
               if (lwd(i) == continuum_orbital_b) lwb = i
            end do

            !No occurances
            if (max(lwa, lwb) == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

                do prototypes = 1,num_prototypes
                    mat_coeffs(prototypes) = mat_coeffs(prototypes) &
                                           + master_prototype(prototypes) % get_coefficient(ido) * integral_coeff
                end do
                !expand_off_diagonal_eval = expand_off_diagonal_eval+this%evaluate_integrals_singular(integral,integral_coeff,NORMAL_PHASE)
                cycle
            else if (lwb == 0) then
                lwd(lwa) = continuum_orbital_a + ja - 1
            else if (lwa == 0) then
                lwd(lwb) = continuum_orbital_b + jb - 2
            else
                lwd(lwa) = continuum_orbital_a + ja - 1
                lwd(lwb) = continuum_orbital_b + jb - 2
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)

            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do prototypes = 1,num_prototypes
                mat_coeffs(prototypes) = mat_coeffs(prototypes) + master_prototype(prototypes)%get_coefficient(ido) * integral_coeff
            end do

        end do

    end subroutine expand_off_diagonal_eval


    !> \brief   Expands the off-diagonal same/similar symmetry prototype symbols and evaluates for a single matrix element for co-ordinate offset (ja,jb)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this              Hamiltonian object to query.
    !> \param[in] continuum_orbital_a  Prototype continuum orbital for target symmetry a
    !> \param[in] continuum_orbital_b  Prototype continuum orbital for target symmetry b
    !> \param[in] master_prototype     Prototype symbolic elements to be expanded
    !> \param[in] ja                   Desired continuum orbital of target symmetry a
    !> \param[in] jb                   Desired continuum orbital of target symmetry b
    !> \param[out] mat_coeffs          Output array for the matrix element.
    !>
    subroutine contr_expand_off_diagonal_eval (this, continuum_orbital_a, continuum_orbital_b, master_prototype, ja, jb, mat_coeffs)
        class(Contracted_Hamiltonian)                       :: this
        class(ContractedSymbolicElementVector), intent(in)  :: master_prototype
        integer,                                intent(in)  :: continuum_orbital_a, continuum_orbital_b, ja, jb
        real(wp),                               intent(out) :: mat_coeffs(:,:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff

        integer :: lwd(8), lwa, lwb, ido, i, prototypes, num_prototypes, num_integrals, num_targets_1, num_targets_2, n1, n2


        num_integrals = master_prototype % get_size()
        num_targets_1 = master_prototype % get_num_targets_sym1()
        num_targets_2 = master_prototype % get_num_targets_sym2()
        mat_coeffs = 0

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            integral = master_prototype % get_integral_label(ido)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp

            lwa = 0
            lwb = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital_a) lwa = i
               if (lwd(i) == continuum_orbital_b) lwb = i
            end do

            !No occurances
            if (max(lwa, lwb) == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

                do n1 = 1, num_targets_1
                    do n2 = 1, num_targets_2
                        mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                    end do
                end do
                !expand_off_diagonal_eval = expand_off_diagonal_eval+this%evaluate_integrals_singular(integral,integral_coeff,NORMAL_PHASE)
                cycle
            else if (lwb == 0) then
                lwd(lwa) = continuum_orbital_a + ja - 1
            else if (lwa == 0) then
                lwd(lwb) = continuum_orbital_b + jb - 2
            else
                lwd(lwa) = continuum_orbital_a + ja - 1
                lwd(lwb) = continuum_orbital_b + jb - 2
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)

            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do n1 = 1, num_targets_1
                do n2 = 1, num_targets_2
                    mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                end do
            end do

        end do

    end subroutine contr_expand_off_diagonal_eval


    !> \brief   Expands the off-diagonal generic prototype symbols and evaluates for a single matrix element for co-ordinate (ja,jb)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this              Hamiltonian object to query.
    !> \param[in] continuum_orbital_a  Prototype continuum orbital for target symmetry a
    !> \param[in] continuum_orbital_b  Prototype continuum orbital for target symmetry b
    !> \param[in] master_prototype     Prototype symbolic elements to be expanded
    !> \param[in] ja                   Desired continuum orbital of target symmetry a
    !> \param[in] jb                   Desired continuum orbital of target symmetry b
    !> \param[out] mat_coeffs          Output array for the matrix element.
    !>
    subroutine expand_off_diagonal_gen_eval (this, continuum_orbital_a, continuum_orbital_b, master_prototype, ja, jb, mat_coeffs)
        class(Contracted_Hamiltonian)                       :: this
        class(ContractedSymbolicElementVector), intent(in)  :: master_prototype
        integer,                                intent(in)  :: continuum_orbital_a, continuum_orbital_b, ja, jb
        real(wp),                               intent(out) :: mat_coeffs(:,:)

        integer(longint) :: integral(NIDX)
        real(wp)         :: integral_coeff
        integer          :: lwd(8), lwa, lwb, ido, i, n1, n2, num_targets_1, num_targets_2, ii, num_integrals

        num_targets_1 = master_prototype % get_num_targets_sym1()
        num_targets_2 = master_prototype % get_num_targets_sym2()
        num_integrals = master_prototype % get_size()
        mat_coeffs = 0

        if (num_integrals == 0) return

        do ido = 1, num_integrals

            !Get the integral
            integral = master_prototype % get_integral_label(ido)
            call unpack_ints(integral, lwd)

            integral_coeff = 1.0_wp

            lwa = 0
            lwb = 0

            do i = 1, 4
               if (lwd(i) == continuum_orbital_a) lwa = i
               if (lwd(i) == continuum_orbital_b) lwb = i
            end do

            !No occurances
            if (max(lwa, lwb) == 0) then
                integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)
                do n1 = 1,num_targets_1
                    do n2 = 1,num_targets_2
                        mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                    end do
                end do
                cycle
            else if(lwb == 0)then
                lwd(lwa) = continuum_orbital_a + ja - 1
            else if (lwa == 0) then
                lwd(lwb) = continuum_orbital_b + jb - 1
            else
                lwd(lwa) = continuum_orbital_a + ja - 1
                lwd(lwb) = continuum_orbital_b + jb - 1
            end if

            call pack_ints(lwd(1), lwd(2), lwd(3), lwd(4), lwd(5), 0, 0, 0, integral)

            integral_coeff = this % evaluate_integrals_singular(integral, 1.0_wp, NORMAL_PHASE)

            do n1 = 1, num_targets_1
                do n2 = 1, num_targets_2
                    mat_coeffs(n1,n2) = mat_coeffs(n1,n2) + master_prototype % get_coefficient(ido, n1, n2) * integral_coeff
                end do
            end do

        end do

    end subroutine expand_off_diagonal_gen_eval


    !> \brief   (Not used) Compresses several prototypes into a single prototype at position 1
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> To be used for the OpenMP implementation by combining multiple symbols.
    !>
    !> \param[inout] this     Hamiltonian object to query.
    !> \param[in] prototypes  A symbolic prototype array to be compressed.
    !>
    subroutine reduce_prototypes (this, prototypes)
        class(Contracted_Hamiltonian), intent(in) :: this
        class(SymbolicElementVector),  intent(in) :: prototypes(:)

        integer :: num_prototypes, ido

        num_prototypes = size(prototypes)

        if (num_prototypes <= 1) return

        do ido = 2, num_prototypes
            call prototypes(1) % add_symbols(prototypes(ido))
        end do

    end subroutine reduce_prototypes


    !> \brief   Gets the starting index configuration state functions at a specific symmetry
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this   Hamiltonian object to query.
    !> \param[in] symmetry  Target symmetry.
    !>
    !> \result The starting CSF index for a given target symmetry
    !>
    integer function get_starting_index (this, symmetry)
        class(Contracted_Hamiltonian), intent(in) :: this
        integer,                       intent(in) :: symmetry

        integer :: ido

        get_starting_index = 1
        do ido = 2, symmetry
            get_starting_index = get_starting_index + this % options % num_ci_target_sym(ido - 1) * this % csf_skip
        end do

    end function get_starting_index


    !> \brief   Simply a wrapper to aid in clarity
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this  Hamiltonian object to query.
    !> \param[in]    i     Target symmetry.
    !> \param[in]    m     Target state.
    !> \param[in]    n     Target configuration.
    !>
    !> \result target coefficient
    !>
    real(wp) function get_target_coeff (this, i, m, n)
         class(Contracted_Hamiltonian), intent(in) :: this

         integer :: i, m, n

         get_target_coeff = this % rmat_ci(i) % eigenvectors(m, n)

    end function get_target_coeff


    subroutine fast_contract_class_1_3_diag (this, i1, num_targets, ido, temp_prototype, master_prototypes)
        class(Contracted_Hamiltonian)            :: this
        integer,                      intent(in) :: i1, ido, num_targets
        class(SymbolicElementVector), intent(in) :: temp_prototype, master_prototypes(num_targets * (num_targets + 1) / 2)

        integer(longint) :: label(NIDX)

        integer  :: found_idx, ii, n1, n2, num_symbols, jdo
        real(wp) :: alpha(num_targets * (num_targets + 1) / 2), coeff

        ii = 1

        !get our coeffs
        do n1 = 1, num_targets
            do n2 = n1, num_targets
                ! Apply our coefficeints for each combination to the master, Left hand of Eq. 6
                alpha(ii) = this % get_target_coeff(i1, ido, n1) * this % get_target_coeff(i1, ido, n2)
                !call master_prototype(ii)%add_symbols(temp_prototype,alpha)
                ii = ii + 1
            end do
        end do

        num_symbols = temp_prototype % get_size()

        do jdo = 1, num_symbols
            call temp_prototype % get_coeff_and_integral(jdo, coeff, label)
            !We only need to perform the expensive check once per insertion
            if (master_prototypes(1) % check_same_integral(label, found_idx)) then
                ii = 1
                !get our coeffs
                do n1 = 1, num_targets
                    do n2 = n1, num_targets
                        call master_prototypes(ii) % modify_coeff(found_idx, alpha(ii) * coeff)
                        ii = ii + 1
                    end do
                end do
            else
                ii = 1
                !get our coeffs
                do n1 = 1, num_targets
                    do n2 = n1, num_targets
                        !Insert without a check
                        call master_prototypes(ii) % insert_symbol(label, alpha(ii) * coeff, .false.)
                        ii = ii + 1
                    end do
                end do
            end if
        end do

    end subroutine fast_contract_class_1_3_diag


    subroutine fast_contract_class_1_3_offdiag (this, i1, num_targets, s1, s2, temp_prototype, master_prototypes)
        class(Contracted_Hamiltonian)            :: this
        integer,                      intent(in) :: i1, s1, s2, num_targets
        class(SymbolicElementVector), intent(in) :: temp_prototype, master_prototypes(num_targets * (num_targets + 1) / 2)

        integer(longint) :: label(NIDX)

        integer  :: found_idx, ii, n1, n2, num_symbols, jdo
        real(wp) :: alpha(num_targets * (num_targets + 1) / 2), coeff

        ii = 1

        !get our coeffs
        do n1 = 1, num_targets
            do n2 = n1, num_targets
                ! Apply our coefficeints for each combination to the master, Left hand of Eq. 6
                alpha(ii) = this % get_target_coeff(i1, s1, n1) * this % get_target_coeff(i1, s2, n2) &
                            + this % get_target_coeff(i1, s2, n1) * this % get_target_coeff(i1, s1, n2)
                !call master_prototype(ii)%add_symbols(temp_prototype,alpha)
                ii = ii + 1
            end do
        end do

        num_symbols = temp_prototype % get_size()

        do jdo = 1, num_symbols
            call temp_prototype % get_coeff_and_integral(jdo, coeff, label)
            !We only need to perform the expensive check once per insertion
            if (master_prototypes(1) % check_same_integral(label, found_idx)) then
                ii = 1
                !get our coeffs
                do n1 = 1, num_targets
                    do n2 = n1, num_targets
                        call master_prototypes(ii) % modify_coeff(found_idx, alpha(ii) * coeff)
                        ii = ii + 1
                    end do
                end do
            else
                ii = 1
                !get our coeffs
                do n1 = 1, num_targets
                    do n2 = n1, num_targets
                        !Insert without a check
                        call master_prototypes(ii) % insert_symbol(label, alpha(ii) * coeff, .false.)
                        ii = ii + 1
                    end do
                end do
            end if
        end do

    end subroutine fast_contract_class_1_3_offdiag


    subroutine fast_contract_class_567 (this, i1, i2, num_targets1, num_targets2, s1, s2, temp_prototype, master_prototypes)
        class(Contracted_Hamiltonian)            :: this
        integer,                      intent(in) :: i1, i2, s1, s2, num_targets1, num_targets2
        class(SymbolicElementVector), intent(in) :: temp_prototype, master_prototypes(num_targets1, num_targets2)

        integer(longint) :: label(NIDX)

        integer  :: found_idx, ii, n1, n2, num_symbols, jdo
        real(wp) :: alpha(num_targets1, num_targets2), coeff

        alpha = 0
        !get our coeffs
        do n1 = 1, num_targets1
            do n2 = 1, num_targets2
                ! Apply our coefficeints for each combination to the master, Left hand of Eq. 6
                alpha(n1,n2) = this % get_target_coeff(i1, s1, n1) * this % get_target_coeff(i2, s2, n2)
                !call master_prototype(ii)%add_symbols(temp_prototype,alpha)
            end do
        end do

        num_symbols = temp_prototype % get_size()

        do jdo = 1, num_symbols
            call temp_prototype % get_coeff_and_integral(jdo, coeff, label)
            !We only need to perform the expensive check once per insertion
            if (master_prototypes(1,1) % check_same_integral(label, found_idx)) then
                !get our coeffs
                do n1 = 1, num_targets1
                    do n2 = 1, num_targets2
                        call master_prototypes(n1,n2) % modify_coeff(found_idx, alpha(n1,n2) * coeff)
                    end do
                end do
            else
                !get our coeffs
                do n1 = 1, num_targets1
                    do n2 = 1, num_targets2
                        !Insert without a check
                        call master_prototypes(n1,n2) % insert_symbol(label, alpha(n1,n2) * coeff, .false.)
                    end do
                end do
            end if
        end do

    end subroutine fast_contract_class_567

end module Contracted_Hamiltonian_Module
