program gauss_quad
    use iso_fortran_env, only: output_unit
    implicit none

    real(8) :: exact_val, approx_val
    integer :: q

    ! 1. 手計算による厳密解の算出
    exact_val = 46.0_8 / 15.0_8

    write(output_unit, '(A, G0)') "Exact value (Hand calculation) : ", exact_val
    write(output_unit, *) "------------------------------------------------"

    ! 2. 1次〜4次までのGauss求積を実行し、結果を出力
    do q = 1, 4
        call calc_gauss(q, approx_val)

        write(output_unit, '(A, I1, A, G0)') "Gauss q = ", q, " : ", approx_val
        write(output_unit, '(A, G0)')        "Error       : ", abs(approx_val - exact_val)
        write(output_unit, *) "------------------------------------------------"
    end do

contains

    ! ========================================================
    ! 積分対象の関数 f(x)
    ! ========================================================
    real(8) function f(x)
        real(8), intent(in) :: x
        f = 1.0_8 + x + x**2_8 + x**3_8 + x**4_8
    end function f

    ! ========================================================
    ! Gauss求積を実行するサブルーチン
    ! ========================================================
    subroutine calc_gauss(q, result)
        integer, intent(in) :: q
        real(8), intent(out) :: result

        ! 積分点(xi)と重み(w)を格納する動的配列
        real(8), allocatable :: xi(:), w(:)
        integer :: i

        result = 0.0_8

        ! 次数(q)に応じて配列のサイズを確保し、表の値を代入
        select case(q)
            case(1)
                allocate(xi(1), w(1))
                xi = [0.0_8]
                w  = [2.0_8]

            case(2)
                allocate(xi(2), w(2))
                xi = [-0.5773502692_8, 0.5773502692_8]
                w  = [1.0_8, 1.0_8]

            case(3)
                allocate(xi(3), w(3))
                xi = [-0.774596697_8, 0.0_8, 0.774596697_8]
                w  = [0.5555555556_8, 0.8888888889_8, 0.5555555556_8]

            case(4)
                allocate(xi(4), w(4))
                xi = [-0.8611363116_8, -0.3399810436_8, 0.3399810436_8, 0.8611363116_8]
                w  = [0.3478548451_8, 0.6521451549_8, 0.6521451549_8, 0.3478548451_8]

            case default
                write(output_unit, '(A)') "Error: Unsupported degree for Gauss Quadrature."
                stop 1, quiet=.true.
        end select

        ! 求積計算の実行： I_q = Σ f(xi) * Wi
        do i = 1, q
            result = result + f(xi(i)) * w(i)
        end do

        ! メモリの解放
        deallocate(xi, w)
    end subroutine calc_gauss

end program gauss_quad
