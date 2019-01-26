!  compile with F90 or newer, gfortran >= 4.8 (checked)

subroutine arr_pow_d(arr_in,  size, alpha, beta, pows)
!		double version calculate (alpha*arr+ beta)^pows
		implicit none
		integer, intent(in) :: size
		real(kind = 8)  :: arr_in(size)
		real(kind = 8), intent(in) :: alpha, beta, pows
!		real(kind = 8), intent(out)  :: arr_out(size)		

		arr_in = (alpha*arr_in+beta)**pows

end subroutine arr_pow_d
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_pow_f(arr_in, arr_out, size, alpha, beta, pows)
!		float version calculate (alpha*arr+ beta)^pows
		implicit none
		integer, intent(in)::size
		real(kind = 4), intent(in)  :: arr_in(size)
		real(kind = 4), intent(in) :: alpha, beta, pows
		real(kind = 4), intent(out)  :: arr_out(size)		

		arr_out = (alpha*arr_in+ beta)**pows

end subroutine arr_pow_f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_base_d(arr_in, arr_out, size, alpha, beta, base)
!		double version calculate base^(alpha*arr + beta)
		implicit none
		integer, intent(in) :: size
		real(kind = 8), intent(in)  :: arr_in(size)
		real(kind = 8), intent(in) :: alpha, beta, base
		real(kind = 8), intent(out)  :: arr_out(size)

		arr_out = base**(alpha*arr_in + beta)

end subroutine arr_base_d
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_base_f(arr_in, arr_out, size, alpha, beta, base)
!		float version calculate base^(alpha*arr + beta)
		implicit none
		integer, intent(in):: size
		real(kind = 4), intent(in)  :: arr_in(size)
		real(kind = 4), intent(in) :: alpha, beta, base
		real(kind = 4), intent(out)  :: arr_out(size)

		arr_out = base**(alpha*arr_in + beta)

end subroutine arr_base_f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_exp_d(arr_in, arr_out, size, alpha, beta)
!		double version calculate exp^(alpha*arr+ beta)
		implicit none
		integer, intent(in) :: size
		real(kind = 8), intent(in) :: alpha, beta
		real(kind = 8), intent(in)  :: arr_in(size)		
		real(kind = 8), intent(out)  :: arr_out(size)		

		arr_out = EXP(alpha*arr_in+beta)

end subroutine arr_exp_d
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_exp_f(arr_in, arr_out, size, alpha, beta)
!		float version calculate exp^(alpha*arr)
		implicit none
		integer, intent(in) :: size
		real(kind = 4), intent(in) :: alpha, beta
		real(kind = 4), intent(in)  :: arr_in(size)		
		real(kind = 4), intent(out)  :: arr_out(size)

		arr_out = EXP(alpha*arr_in+ beta)

end subroutine arr_exp_f
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_sum_d(arr_in, total, size, alpha, beta)
!		double version calculate SUM(arr*alpha + beta)
		implicit none
		integer, intent(in) :: size
		real(kind = 8), intent(in) :: alpha, beta
		real(kind = 8), intent(out) :: total
		real(kind = 8), intent(in)  :: arr_in(size)

		total = SUM(arr_in*alpha+beta)

end subroutine arr_sum_d
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine arr_sum_f(arr_in, total, size, alpha, beta)
!		float version calculate SUM(arr*alpha + beta)
		implicit none
		integer, intent(in) :: size
		real(kind = 4), intent(in) :: alpha, beta
		real(kind = 4), intent(out) :: total
		real(kind = 4), intent(in)  :: arr_in(size)

		total = SUM(arr_in*alpha + beta)

end subroutine arr_sum_f
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc