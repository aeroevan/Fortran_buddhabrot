program buddhabrot
    implicit none

    character (len=*), parameter :: filename = 'buddhabrot.ppm'

    integer, parameter :: grid_resolution = 512, n_max = 100
    real, parameter :: intensity = 2048.

    integer, dimension(grid_resolution, grid_resolution) :: &
        exposure_map_R, exposure_map_G, exposure_map_B

    integer :: max_exposure_R, max_exposure_G, max_exposure_B

    real :: x, y
    integer :: x_tmp, y_tmp, y_tmp_m
    real, parameter :: x_min = -1.0, x_max = 2.0, y_min = -1.3, y_max = 1.3
    real, parameter :: x_delta = x_max - x_min, x_delta_inv = 1.0/x_delta
    real, parameter :: y_delta = y_max - y_min, y_delta_inv = 1.0/y_delta
    complex :: z, c

    integer, parameter :: batch_size = 10000000
    integer :: i, j_low, j_high, j, iter

    exposure_map_R = 0
    exposure_map_G = 0
    exposure_map_B = 0

    !$omp parallel do default(private), shared(exposure_map_R, exposure_map_G, exposure_map_B)
    do i = 1, batch_size
        call random_number(x)
        call random_number(y)

        z = cmplx(0, 0)
        c = cmplx(2.5*x - 2.0, 1.3*y)

        if (not_in_M_set(c, n_max)) then
            do iter = 1, n_max
                z = z*z + c
                x_tmp = int(grid_resolution * (real(z) + x_max) * x_delta_inv)
                y_tmp = int(grid_resolution * (aimag(z) + y_max) * y_delta_inv)
                y_tmp_m = int(grid_resolution - y_tmp)
                if ((x_tmp > 0) .and. (x_tmp < grid_resolution) .and. &
                    (y_tmp > 0) .and. (y_tmp < grid_resolution)) then
                    if ((iter > 2) .and. (iter < 50)) then
                        exposure_map_B(x_tmp, y_tmp)    = exposure_map_B(x_tmp, y_tmp) + 1
                        exposure_map_B(x_tmp, y_tmp_m)  = exposure_map_B(x_tmp, y_tmp_m) + 1
                    end if
                    if ((iter > 25) .and. (iter < 75)) then
                        exposure_map_G(x_tmp, y_tmp)    = exposure_map_G(x_tmp, y_tmp) + 1
                        exposure_map_G(x_tmp, y_tmp_m)  = exposure_map_G(x_tmp, y_tmp_m) + 1
                    end if
                    if ((iter > 50) .and. (iter < 100)) then
                        exposure_map_R(x_tmp, y_tmp)    = exposure_map_R(x_tmp, y_tmp) + 1
                        exposure_map_R(x_tmp, y_tmp_m)  = exposure_map_R(x_tmp, y_tmp_m) + 1
                    end if
                end if
            end do
        end if
    end do
    !$omp end parallel do

    max_exposure_R = maxval(exposure_map_R)
    max_exposure_G = maxval(exposure_map_G)
    max_exposure_B = maxval(exposure_map_B)

    exposure_map_R = int((exposure_map_R / real(max_exposure_R))*intensity)
    exposure_map_G = int((exposure_map_G / real(max_exposure_G))*intensity)
    exposure_map_B = int((exposure_map_B / real(max_exposure_B))*intensity)

    open(unit=33, file=filename, status='replace', form='formatted', &
        access='sequential')
    write(33, '(a2)')       'P3'
    write(33, '(i5,2x,i5)') grid_resolution, grid_resolution
    write(33, '(i5)')       int(intensity)

    do i = 1, grid_resolution
        do j_low = 1, grid_resolution, 4
            j_high = min(j_low + 3, grid_resolution)
            write(33, '(12i5)') &
                (exposure_map_R(i, j), exposure_map_G(i, j), exposure_map_B(i, j), &
                    j = j_low, j_high)
        end do
    end do

    close(33)
contains
    pure function not_in_M_set(c, n_max)
        implicit none
        complex, intent(in) :: c
        integer, intent(in) :: n_max

        real, parameter :: escape_orbit = 2.0

        integer :: n
        complex :: z
        logical :: not_in_M_set

        z = cmplx(0,0)
        n = 0
        if ((abs(c - cmplx(-1,0)) < 0.25) .or. (abs(1.0 - sqrt(1.0 - 4*c)) < 1.0)) then
            not_in_M_set = .false.
        else
            do while ((abs(z) < escape_orbit) .and. (n < n_max))
                z = z*z + c
                n = n + 1
            end do
            if (n >= n_max) then
                not_in_M_set = .false.
            else
                not_in_M_set = .true.
            end if
        end if
    end function not_in_M_set
end program buddhabrot
