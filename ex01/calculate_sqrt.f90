program sqrt_newton
    use iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none

    real(8) :: a, x_old, x_new, epsilon
    integer :: max_iter, i

    epsilon = 1.0e-10_8
    max_iter = 100

    write(output_unit, '(A)', advance='no') "a: "
    read(input_unit, *) a

    if (a < 0.0_8) then
        write(error_unit, '(G0)') "Error: Input must be a non-negative number."
        stop 1, quiet=.true.
    else if (a == 0.0_8) then
        write(output_unit, '(A, G0)') "result: ", 0.0_8
        stop
    end if

    x_new = a

    do i = 1, max_iter
        x_old = x_new
        x_new = 0.5_8 * (x_old + a / x_old)

        if (abs(x_new - x_old) < epsilon) exit
    end do

    write(output_unit, '(A, G0)') "result: ", x_new

end program sqrt_newton
