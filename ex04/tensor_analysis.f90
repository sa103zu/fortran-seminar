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
    ans_scalar = dot_product(u, v)
    write(output_unit, *) "1) u . v"
    write(output_unit, '(G0)') ans_scalar

    write(output_unit, *) ""

    ! 2) u × v (外積)
    ans_vec = cross_product(u, v)
    call print_vec("2) u x v", ans_vec)

    write(output_unit, *) ""

    ! 3) u (x) v (テンソル積 / 直積)
    ans_mat = dyadic_product(u, v)
    call print_mat("3) u (x) v", ans_mat)

    write(output_unit, *) ""

    ! 4) Su (テンソルとベクトルの積)
    ans_vec = matmul(S_mat, u)
    call print_vec("4) S u", ans_vec)

    write(output_unit, *) ""

    ! 5) S × u (テンソルとベクトルの外積)
    ans_mat = tensor_cross_vec(S_mat, u)
    call print_mat("5) S x u", ans_mat)

    write(output_unit, *) ""

    ! 6) ST (テンソル同士の積)
    ans_mat = matmul(S_mat, T_mat)
    call print_mat("6) S T", ans_mat)

    write(output_unit, *) ""

    ! 7) S・T (テンソルの内積)
    ! Σ(S_ij * T_ij) なので、配列の要素ごとの積の和をとる
    ans_scalar = sum(S_mat * T_mat)
    write(output_unit, *) "7) S . T"
    write(output_unit, '(G0)') ans_scalar

    ! ==========================================
    ! 3. 課題8: 演習問題 1.2の計算 (1〜3のみ)
    ! ==========================================
    write(output_unit, *) ""
    write(output_unit, *) "========== 演習問題 1.2 =========="
    write(output_unit, *) ""

    ! 1) u・T v
    ans_scalar = dot_product(u, matmul(T_mat, v))
    write(output_unit, *) "1.2-1) u . T v"
    write(output_unit, '(G0)') ans_scalar

    write(output_unit, *) ""

    ! 2) S^T T^T u
    ans_vec = matmul(transpose(S_mat), matmul(transpose(T_mat), u))
    call print_vec("1.2-2) S^T T^T u", ans_vec)

    write(output_unit, *) ""

    ! 3) v・T S u
    ans_scalar = dot_product(v, matmul(T_mat, matmul(S_mat, u)))
    write(output_unit, *) "1.2-3) v . T S u"
    write(output_unit, '(G0)') ans_scalar

contains

    ! ==========================================
    ! サブルーチン・関数群 (すべて整数型で再定義)
    ! ==========================================

    ! 交代記号（エディントンのイプシロン） eps_ijk を返す関数
    integer function eps(i, j, k)
        integer, intent(in) :: i, j, k
        if (i == j .or. j == k .or. k == i) then
            eps = 0
        else if ((i==1 .and. j==2 .and. k==3) .or. &
                 (i==2 .and. j==3 .and. k==1) .or. &
                 (i==3 .and. j==1 .and. k==2)) then
            eps = 1
        else
            eps = -1
        end if
    end function eps

    ! ベクトル同士の外積: c = a x b
    function cross_product(a, b) result(c)
        integer, intent(in) :: a(3), b(3)
        integer :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    ! ベクトル同士のテンソル積（直積）: C_ij = a_i * b_j
    function dyadic_product(a, b) result(C)
        integer, intent(in) :: a(3), b(3)
        integer :: C(3,3)
        integer :: i, j
        do i = 1, 3
            do j = 1, 3
                C(i,j) = a(i) * b(j)
            end do
        end do
    end function dyadic_product

    ! テンソルとベクトルの外積: C_ij = A_ik * eps_jkl * b_l
    function tensor_cross_vec(A, b) result(C)
        integer, intent(in) :: A(3,3), b(3)
        integer :: C(3,3)
        integer :: i, j, k, l
        C = 0
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do l = 1, 3
                        if (eps(j,k,l) /= 0) then
                            C(i,j) = C(i,j) + A(i,k) * eps(j,k,l) * b(l)
                        end if
                    end do
                end do
            end do
        end do
    end function tensor_cross_vec

    ! ベクトルをきれいに表示するサブルーチン
    subroutine print_vec(title, v)
        character(len=*), intent(in) :: title
        integer, intent(in) :: v(3)
        write(output_unit, '(A)') title
        write(output_unit, '(3(G0, 2X))') v
    end subroutine print_vec

    ! テンソル（行列）をきれいに表示するサブルーチン
    subroutine print_mat(title, A)
        character(len=*), intent(in) :: title
        integer, intent(in) :: A(3,3)
        integer :: r
        write(output_unit, '(A)') title
        do r = 1, 3
            write(output_unit, '(3(G0, 2X))') A(r,:)
        end do
    end subroutine print_mat

end program tensor_calc
