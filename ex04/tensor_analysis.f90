program tensor_calc
    use iso_fortran_env, only: output_unit
    implicit none

    ! 変数宣言 (すべて整数型を使用)
    integer :: u(3), v(3)
    integer :: S_mat(3,3), T_mat(3,3)

    ! 作業用・結果格納用
    integer :: ans_scalar
    integer :: ans_vec(3)
    integer :: ans_mat(3,3)
    integer :: i, j, k, l  ! 指標表示のループ計算用変数

    ! ==========================================
    ! 1. 初期値の設定
    ! ==========================================
    u = [1, 2, 3]
    v = [4, 5, 6]

    ! Fortranは列優先(Column-major)なので、行ごとに明示的に代入する
    S_mat(1,:) = [1, 4, 6]
    S_mat(2,:) = [4, 2, 5]
    S_mat(3,:) = [6, 5, 3]

    T_mat(1,:) = [ 7, 10, 12]
    T_mat(2,:) = [10,  8, 11]
    T_mat(3,:) = [12, 11,  9]

    ! ==========================================
    ! 2. 課題1〜7の計算
    ! ==========================================
    write(output_unit, *) "========== 基本課題 =========="
    write(output_unit, *) ""

    ! 1) u・v (内積)
    ! 指標表示: u_i v_i
    ans_scalar = 0
    do i = 1, 3
        ans_scalar = ans_scalar + u(i) * v(i)
    end do
    call print_scalar("1) u . v", "u_i v_i", ans_scalar)

    write(output_unit, *) ""

    ! 2) u × v (外積)
    ! 指標表示: e_ijk u_j v_k
    ans_vec = 0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                ans_vec(i) = ans_vec(i) + eps(i,j,k) * u(j) * v(k)
            end do
        end do
    end do
    call print_vec("2) u x v", "e_ijk u_j v_k", ans_vec)

    write(output_unit, *) ""

    ! 3) u (x) v (テンソル積 / 直積)
    ! 指標表示: u_i v_j
    ans_mat = 0
    do i = 1, 3
        do j = 1, 3
            ans_mat(i,j) = u(i) * v(j)
        end do
    end do
    call print_mat("3) u (x) v", "u_i v_j", ans_mat)

    write(output_unit, *) ""

    ! 4) Su (テンソルとベクトルの積)
    ! 指標表示: S_ij u_j
    ans_vec = 0
    do i = 1, 3
        do j = 1, 3
            ans_vec(i) = ans_vec(i) + S_mat(i,j) * u(j)
        end do
    end do
    call print_vec("4) S u", "S_ij u_j", ans_vec)

    write(output_unit, *) ""

    ! 5) S × u (テンソルとベクトルの外積)
    ! 指標表示: S_ik e_jkl u_l
    ans_mat = 0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    ans_mat(i,j) = ans_mat(i,j) + S_mat(i,k) * eps(j,k,l) * u(l)
                end do
            end do
        end do
    end do
    call print_mat("5) S x u", "S_ik e_jkl u_l", ans_mat)

    write(output_unit, *) ""

    ! 6) ST (テンソル同士の積)
    ! 指標表示: S_ik T_kj
    ans_mat = 0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                ans_mat(i,j) = ans_mat(i,j) + S_mat(i,k) * T_mat(k,j)
            end do
        end do
    end do
    call print_mat("6) S T", "S_ik T_kj", ans_mat)

    write(output_unit, *) ""

    ! 7) S・T (テンソルの内積)
    ! 指標表示: S_ij T_ij
    ans_scalar = 0
    do i = 1, 3
        do j = 1, 3
            ans_scalar = ans_scalar + S_mat(i,j) * T_mat(i,j)
        end do
    end do
    call print_scalar("7) S . T", "S_ij T_ij", ans_scalar)


    ! ==========================================
    ! 3. 課題8: 演習問題 1.2の計算 (1〜3のみ)
    ! ==========================================
    write(output_unit, *) ""
    write(output_unit, *) "========== 演習問題 1.2 =========="
    write(output_unit, *) ""

    ! 1) u・T v
    ! 指標表示: u_i v_j T_ij
    ans_scalar = 0
    do i = 1, 3
        do j = 1, 3
            ans_scalar = ans_scalar + u(i) * v(j) * T_mat(i,j)
        end do
    end do
    call print_scalar("1.2-1) u . T v", "u_i v_j T_ij", ans_scalar)

    write(output_unit, *) ""

    ! 2) S^T T^T u
    ! 指標表示: S_ji u_k T_kj
    ans_vec = 0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                ans_vec(i) = ans_vec(i) + S_mat(j,i) * u(k) * T_mat(k,j)
            end do
        end do
    end do
    call print_vec("1.2-2) S^T T^T u", "S_ji u_k T_kj", ans_vec)

    write(output_unit, *) ""

    ! 3) v・T S u
    ! 指標表示: S_ij T_ki u_j v_k
    ans_scalar = 0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                ans_scalar = ans_scalar + S_mat(i,j) * T_mat(k,i) * u(j) * v(k)
            end do
        end do
    end do
    call print_scalar("1.2-3) v . T S u", "S_ij T_ki u_j v_k", ans_scalar)

contains

    ! ==========================================
    ! サブルーチン・関数群
    ! ==========================================

    ! 交代記号 eps_ijk を返す関数
    integer function eps(i, j, k)
        integer, intent(in) :: i, j, k
        eps = (i - j) * (j - k) * (k - i) / 2
    end function eps

    ! スカラー出力
    subroutine print_scalar(title, idx_str, s)
        character(len=*), intent(in) :: title, idx_str
        integer, intent(in) :: s
        write(output_unit, '(A)') title
        write(output_unit, '(A)') "  [Index notation]: " // idx_str
        write(output_unit, '("    ", I6)') s
    end subroutine print_scalar

    ! ベクトル出力
    subroutine print_vec(title, idx_str, v)
        character(len=*), intent(in) :: title, idx_str
        integer, intent(in) :: v(3)
        write(output_unit, '(A)') title
        write(output_unit, '(A)') "  [Index notation]: " // idx_str
        write(output_unit, '("    ", 3(I6, 2X))') v
    end subroutine print_vec

    ! テンソル出力
    subroutine print_mat(title, idx_str, A)
        character(len=*), intent(in) :: title, idx_str
        integer, intent(in) :: A(3,3)
        integer :: r
        write(output_unit, '(A)') title
        write(output_unit, '(A)') "  [Index notation]: " // idx_str
        do r = 1, 3
            write(output_unit, '("    ", 3(I6, 2X))') A(r,:)
        end do
    end subroutine print_mat

end program tensor_calc
