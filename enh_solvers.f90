module lib

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sqr_it(x, switch)

    real*8 :: sqr_it, x
    integer*4 :: switch

    if (switch == 1) then
      if (x < 0) then
        sqr_it = -1
      else
        sqr_it = 1
      end if
    else
      sqr_it = x
    end if

  end function sqr_it

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function func_apod(x, xlen, ascale)

    real*8 :: func_apod, x, xlen, ascale
    func_apod = exp(-ascale*(x-xlen/2)**2/xlen**2)

  end function func_apod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function func_apod_filter(x, apod, filter_df)

    real*8 :: func_apod_filter, x, apod, filter_df

    func_apod_filter = 0
    if (x < filter_df) then
      func_apod_filter = 1 - apod
    end if

  end function func_apod_filter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk4_cmt_moire_bwd(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x0, dx_calc, y0, xv_nopo, xv_nopos, env)

    implicit none

    !input variables
    integer*4, intent(in) :: xv_nopo, xv_nopos(xv_nopo)
    real*8, intent(in) :: alpha, detuning, phi, ascale, xlen, x0, dx_calc
    complex*16, intent(in) :: kappa(2)
    integer*4, intent(in) :: cs, pig

    !output variables
    complex*16, intent(inout) :: y0(2), env(2,xv_nopo)

    !internal
    integer*4 :: ii, jj
    real*8 :: x, mag
    complex*16 :: i, u(2), fout(2), k1(2), k2(2), k3(2), k4(2)

    i = (0.0, 1.0)

    u(:) = y0(:)
    env(:,xv_nopo) = u(:)
    x = x0

    do ii = 2,xv_nopo

      do jj = xv_nopos(ii-1)+2,xv_nopos(ii)+1

        call func_moire(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x, u, fout)
        k1(:) = -dx_calc * fout(:)

        call func_moire(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x-dx_calc/2, u(:)+k1(:)/2, fout)
        k2(:) = -dx_calc * fout(:)

        call func_moire(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x-dx_calc/2, u(:)+k2(:)/2, fout)
        k3(:) = -dx_calc * fout(:)

        call func_moire(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x-dx_calc, u(:)+k3(:), fout)
        k4(:) = -dx_calc * fout(:)

        u(:) = u(:) + (k1(:) + 2*k2(:) + 2*k3(:) + k4(:)) / 6

        x = x - dx_calc

      end do

      env(:,xv_nopo-ii+1) = u(:)

      !!!!!! normalise
      mag = env(1,xv_nopo-ii+1)
      env(:,:) = env(:,:) / mag
      u(:) = u(:) / mag

    end do

  end subroutine rk4_cmt_moire_bwd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine func_moire(kappa, alpha, detuning, phi, cs, pig, ascale, xlen, x, u, fout)

    implicit none

    !input variables
    real*8, intent(in) :: alpha, detuning, phi, ascale, xlen, x
    complex*16, intent(in) :: u(2), kappa(2)
    integer*4, intent(in) :: cs, pig

    !output variables
    complex*16, intent(inout) :: fout(2)

    !internal
    real*8 :: apod, pi
    complex*16 :: i, moire

    i = (0.0, 1.0)
    pi = 3.14159265359

    apod = func_apod(x, xlen, ascale)
    moire = i*apod*sqr_it(cos(2*x*alpha+phi+cs*(pi/2-alpha*xlen)),pig)

    fout(1) = kappa(1)*moire*u(2)*exp(-i*detuning*x)
    fout(2) = -kappa(2)*moire*u(1)*exp(i*detuning*x)

  end subroutine func_moire

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk4_cmt_x2_moire_cw_fwd(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x0, dx_calc, y0, xv_nopo, xv_nopos, env)

    implicit none

    !input variables
    integer*4, intent(in) :: xv_nopo, xv_nopos(xv_nopo)
    real*8, intent(in) :: n1, n2, pm_scale, filter, alpha, phi, ascale, xlen, x0, dx_calc
    complex*16, intent(in) :: kappa(4)
    integer*4, intent(in) :: cs, pig

    !output variables
    complex*16, intent(inout) :: y0(4), env(4,xv_nopo)

    !internal
    integer*4 :: ii, jj
    real*8 :: x
    complex*16 :: i, u(4), fout(4), k1(4), k2(4), k3(4), k4(4)

    i = (0.0, 1.0)

    u(:) = y0(:)
    env(:,1) = u(:)
    x = x0

    do ii = 2,xv_nopo

      do jj = xv_nopos(ii-1)+2,xv_nopos(ii)+1

        call func_x2_moire_cw(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x, u, fout)
        k1(:) = dx_calc * fout(:)

        call func_x2_moire_cw(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x+dx_calc/2, u(:)+k1(:)/2, fout)
        k2(:) = dx_calc * fout(:)

        call func_x2_moire_cw(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x+dx_calc/2, u(:)+k2(:)/2, fout)
        k3(:) = dx_calc * fout(:)

        call func_x2_moire_cw(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x+dx_calc, u(:)+k3(:), fout)
        k4(:) = dx_calc * fout(:)

        u(:) = u(:) + (k1(:) + 2*k2(:) + 2*k3(:) + k4(:)) / 6

        x = x + dx_calc

      end do

      env(:,ii) = u(:)

    end do

  end subroutine rk4_cmt_x2_moire_cw_fwd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine func_x2_moire_cw(n1, n2, pm_scale, filter, kappa, alpha, phi, cs, pig, ascale, xlen, x, u, fout)

    implicit none

    !input variables
    real*8, intent(in) :: n1, n2, pm_scale, filter, alpha, phi, ascale, xlen, x
    complex*16, intent(in) :: u(4), kappa(4)
    integer*4, intent(in) :: cs, pig

    !output variables
    complex*16, intent(inout) :: fout(4)

    !internal
    real*8 :: pi, apod1, apod2
    complex*16 :: i, moire, pm

    i = (0.0, 1.0)
    pi = 3.14159265359

    apod1 = func_apod(x, xlen, ascale)
    apod2 = func_apod_filter(x, apod1, filter)
    moire = i*apod1*sqr_it(cos(2*x*alpha+phi+cs*(pi/2-alpha*xlen)),pig)

    fout(1) = kappa(1)*moire*u(2) + pm_scale*conjg(u(1))*u(3)/n1
    fout(2) = -kappa(2)*moire*u(1) + pm_scale*conjg(u(2))*u(4)/ n1
    fout(3) = i*kappa(3)*apod2*u(4) - pm_scale*u(1)*u(1)/n2
    fout(4) = -i*kappa(4)*apod2*u(3) - pm_scale*u(2)*u(2)/n2

  end subroutine func_x2_moire_cw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lib
