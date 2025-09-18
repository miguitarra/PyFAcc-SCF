module basis
    use iso_c_binding
    use basis_struct
    use utils
    implicit none

contains

    subroutine build_basis(atom_symbols, z_num, coords, natoms, basis_file, molecule, shells)
        character(len=*), dimension(natoms), intent(in) :: atom_symbols
        integer(c_int), intent(in) :: z_num(:)
        real(c_double), dimension(3, natoms), intent(in) :: coords
        integer(c_int), intent(in) :: natoms
        character(len=*), intent(in) :: basis_file
        type(Molecule_type), intent(inout) :: molecule
        type(Shell_type), allocatable, intent(out) :: shells(:)

        integer :: i, j, n_shells, ios, unit
        character(len=200) :: line, atom_label, shell_label
        integer :: num_prim
        real(wp) :: scale
        integer, allocatable :: atom_indices(:)
        integer :: n_letter

        type(Shell_type), allocatable :: tmp_shells(:)

        ! initialize molecule atoms
        allocate(molecule%atoms(natoms))
        molecule%n_atoms = natoms
        do i = 1, natoms
            molecule%atoms(i)%x = coords(1, i)
            molecule%atoms(i)%y = coords(2, i)
            molecule%atoms(i)%z = coords(3, i)
            molecule%atoms(i)%z_num = z_num(i)
            molecule%atoms(i)%n_shells = 0
        end do

        ! read basis file
        open(newunit=unit, file=trim(basis_file), status='old', action='read')

        n_shells = 0
        allocate(tmp_shells(0))

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (trim(line) == '' .or. line(1:1) == '!') cycle

            ! start of atom section
            read(line, *, iostat=ios) atom_label
            if (ios /= 0) cycle
            if (trim(atom_label) == '****') cycle

            ! find all indices of this atom_label in the molecule
            n_letter = 0
            allocate(atom_indices(0))
            do i = 1, natoms
                if (trim(atom_symbols(i)) == trim(atom_label)) then
                    n_letter = n_letter + 1
                    call extend_int_array(atom_indices, n_letter)
                    atom_indices(n_letter) = i
                end if
            end do
            if (n_letter == 0) then
                deallocate(atom_indices)
                cycle  ! skip atoms not in molecule
            end if

            ! loop shells for this atom until ****
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, '****') > 0) exit

                read(line, *) shell_label, num_prim, scale
                if (trim(shell_label) == 'S') then
                    call read_shell(unit, atom_indices, 0, num_prim, n_shells, tmp_shells, molecule)
                else if (trim(shell_label) == 'P') then
                    call read_shell(unit, atom_indices, 1, num_prim, n_shells, tmp_shells, molecule)
                else if (trim(shell_label) == 'SP') then
                    call read_sp_shell(unit, atom_indices, num_prim, n_shells, tmp_shells, molecule)
                end if
            end do
            deallocate(atom_indices)
        end do

        close(unit)

        ! return final shells
        call move_alloc(tmp_shells, shells)

    end subroutine build_basis


    subroutine read_shell(unit, atom_indices, lval, num_prim, n_shells, shells, molecule)
        integer, intent(in) :: unit, lval, num_prim
        integer, intent(inout) :: n_shells
        integer, intent(in), dimension(:) :: atom_indices
        type(Shell_type), allocatable, intent(inout) :: shells(:)
        type(Molecule_type), intent(inout) :: molecule
        real(wp) :: expv, coeff


        integer :: n_shells_old, i, idx, n_letter

        n_letter = size(atom_indices)
        n_shells_old = n_shells
        n_shells = n_shells + n_letter
        call extend_shells(shells, n_shells)

        do idx = 1, n_letter
            shells(n_shells_old + idx)%atom_index = atom_indices(idx)
            shells(n_shells_old + idx)%l_num = lval
            shells(n_shells_old + idx)%num_prim = num_prim
            shells(n_shells_old + idx)%x = molecule%atoms(atom_indices(idx))%x
            shells(n_shells_old + idx)%y = molecule%atoms(atom_indices(idx))%y
            shells(n_shells_old + idx)%z = molecule%atoms(atom_indices(idx))%z
            allocate(shells(n_shells_old + idx)%exponents(num_prim))
            allocate(shells(n_shells_old + idx)%coefficients(num_prim))
            molecule%atoms(atom_indices(idx))%n_shells = molecule%atoms(atom_indices(idx))%n_shells + 1
        end do

        do i = 1, num_prim
            read(unit, *) expv, coeff
            do idx = 1, n_letter
                shells(n_shells_old + idx)%exponents(i) = expv
                shells(n_shells_old + idx)%coefficients(i) = coeff
            end do
        end do
    end subroutine read_shell


    subroutine read_sp_shell(unit, atom_indices, num_prim, n_shells, shells, molecule)
        integer, intent(in) :: unit, num_prim
        integer, intent(inout) :: n_shells
        integer, intent(in), dimension(:) :: atom_indices
        type(Shell_type), allocatable, intent(inout) :: shells(:)
        type(Molecule_type), intent(inout) :: molecule

        integer :: n_shells_old, i, idx, n_letter
        real(wp) :: expv, coeffS, coeffP

        n_letter = size(atom_indices)
        ! S part
        n_shells_old = n_shells
        n_shells = n_shells + n_letter
        call extend_shells(shells, n_shells)

        do idx = 1, n_letter
            shells(n_shells_old + idx)%atom_index = atom_indices(idx)
            shells(n_shells_old + idx)%l_num = 0
            shells(n_shells_old + idx)%num_prim = num_prim
            shells(n_shells_old + idx)%x = molecule%atoms(atom_indices(idx))%x
            shells(n_shells_old + idx)%y = molecule%atoms(atom_indices(idx))%y
            shells(n_shells_old + idx)%z = molecule%atoms(atom_indices(idx))%z
            allocate(shells(n_shells_old + idx)%exponents(num_prim))
            allocate(shells(n_shells_old + idx)%coefficients(num_prim))
            molecule%atoms(atom_indices(idx))%n_shells = molecule%atoms(atom_indices(idx))%n_shells + 1
        end do

        ! P part
        n_shells_old = n_shells
        n_shells = n_shells + n_letter
        call extend_shells(shells, n_shells)
        do idx = 1, n_letter
            shells(n_shells_old + idx)%atom_index = atom_indices(idx)
            shells(n_shells_old + idx)%l_num = 1
            shells(n_shells_old + idx)%num_prim = num_prim
            shells(n_shells_old + idx)%x = molecule%atoms(atom_indices(idx))%x
            shells(n_shells_old + idx)%y = molecule%atoms(atom_indices(idx))%y
            shells(n_shells_old + idx)%z = molecule%atoms(atom_indices(idx))%z
            allocate(shells(n_shells_old + idx)%exponents(num_prim))
            allocate(shells(n_shells_old + idx)%coefficients(num_prim))
            molecule%atoms(atom_indices(idx))%n_shells = molecule%atoms(atom_indices(idx))%n_shells + 1
        end do

        do i = 1, num_prim
            read(unit, *) expv, coeffS, coeffP
            ! S shells
            do idx = 1, n_letter
                shells(n_shells - 2 * n_letter + idx)%exponents(i) = expv
                shells(n_shells - 2 * n_letter + idx)%coefficients(i) = coeffS
            end do
            ! P shells
            do idx = 1, n_letter
                shells(n_shells - n_letter + idx)%exponents(i) = expv
                shells(n_shells - n_letter + idx)%coefficients(i) = coeffP
            end do
        end do
    end subroutine read_sp_shell

    subroutine extend_shells(shells, n_new)
        type(Shell_type), allocatable, intent(inout) :: shells(:)
        integer, intent(in) :: n_new
        type(Shell_type), allocatable :: tmp(:)
        if (.not. allocated(shells)) then
            allocate(shells(n_new))
        else if (size(shells) < n_new) then
            call move_alloc(shells, tmp)
            allocate(shells(n_new))
            shells(1:size(tmp)) = tmp
        end if
    end subroutine extend_shells

    subroutine extend_int_array(arr, n_new)
        integer, allocatable, intent(inout) :: arr(:)
        integer, intent(in) :: n_new
        integer, allocatable :: tmp(:)
        if (.not. allocated(arr)) then
            allocate(arr(n_new))
        else if (size(arr) < n_new) then
            call move_alloc(arr, tmp)
            allocate(arr(n_new))
            arr(1:size(tmp)) = tmp
        end if
    end subroutine extend_int_array

end module basis