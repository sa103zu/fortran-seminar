program pi_montecarlo
    use iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none

    real(8) :: x, y, pi_approx, pi_true
    integer(8) :: n, count, i

    write(output_unit, '(A)', advance='no') "n: "
    read(input_unit, *) n

    if (n <= 0) then
        write(error_unit, *) "Error: Input must be a positive integer."
        stop 1, quiet=.true.
    end if

    count = 0
    call random_seed()

    ! モンテカルロシミュレーション
    do i = 1, n
        ! 0.0 <= x < 1.0 の一様乱数を生成
        call random_number(x)
        call random_number(y)

        ! 原点からの距離の二乗が1以下か判定（四分の一円の内部か）
        if (x**2 + y**2 <= 1.0_8) then
            count = count + 1
        end if
    end do

    ! 結果の計算 (整数同士の割り算を避けるため、real(..., 8) で実数に変換)
    pi_approx = 4.0_8 * real(count, 8) / real(n, 8)

    ! 比較用の真の円周率を計算 (arccos(-1) = π)
    pi_true = acos(-1.0_8)

    ! 結果の出力
    write(output_unit, '(A, G0)') "result (Monte Carlo) : ", pi_approx
    write(output_unit, '(A, G0)') "true value (acos(-1)): ", pi_true
    write(output_unit, '(A, G0)') "Error                : ", abs(pi_approx - pi_true)

end program pi_montecarlo
