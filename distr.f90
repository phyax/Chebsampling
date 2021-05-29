module distr
  ! use util
  implicit none
  ! add any new type of distribution data, 
  ! or modify old distribution data in this module.
  ! Each subroutine corresponds to one type of distribution function.

  contains


!----------------------------------------------------------------------
  subroutine force_free(fv, nvxf, nvyf, fgk, np)
    implicit none
    integer, intent(in) :: nvxf, nvyf, np
    real, dimension(nvxf, nvyf), intent(inout) :: fv
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: vt = 1.0, a = 1.0, b = 2.0
    real, parameter :: beta = vt**(-2), ux = 2.0**0.5 * vt, uy = 2.0**0.5 * vt
    real, parameter :: Ax = 0.0, Ay = 0.0
    real, parameter :: vxmax = 5.0, vymax = 5.0, betah = 0.5 * beta 
    real, dimension(nvxf) :: vx
    real, dimension(nvyf) :: vy
    real, dimension(nvxf, nvyf) :: wgt
    integer :: j

    ! distribution data
    vx = linspace(-vxmax, vxmax, nvxf)
    vy = linspace(-vymax, vymax, nvyf)
    do j = 1, nvyf
      fv(:, j) = exp(-betah * (vx**2 + vy(j)**2))&
        & * (exp(beta * uy * (vy(j) + Ay)) + a * cos(beta * ux * (vx + Ax)) + b)
    enddo

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nvxf, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nvyf) = 0.5
    wgt(1, 1) = 0.25
    wgt(nvxf, 1) = 0.25
    wgt(1, nvyf) = 0.25
    wgt(nvxf, nvyf) = 0.25
    fgk = sum(wgt * fv) / real(np)
    return
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine power_maxwellian(fv, nvperf, nvparf, fgk, np)
    implicit none
    integer, intent(in) :: nvperf, nvparf, np
    real, dimension(nvperf, nvparf), intent(inout) :: fv
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: kappa = 0.2
    real, parameter :: vtmpar = 1.0, vtmper = vtmpar * 2.0**0.5
    real, parameter :: vtpper = vtmper * (3.0 / 800.0)**0.5 
    real, parameter :: vtppar = vtpper * 2.0**0.5
    real, parameter :: vpermax = 4.0, vparmax = 4.0
    real, dimension(nvperf) :: vper
    real, dimension(nvparf) :: vpar
    real, dimension(nvperf, nvparf) :: wgt
    integer :: j
    real :: rkappah

    ! distribution data
    vper = linspace(0.0, vpermax, nvperf)
    vpar = linspace(-vparmax, vparmax, nvparf)
    rkappah = 0.5 / kappa
    do j = 1, nvparf
      ! factor (2 * pi * vper) is due to the cylindrical coordinate.
      fv(:, j) = (1.0 + rkappah * ((vper/vtpper)**2&
        & + (vpar(j)/vtppar)**2))**(-kappa - 1.0)&
        & * exp(-0.5 * ((vper/vtmper)**2 + (vpar(j)/vtmpar)**2))&
        & * (6.28318530717959 * vper)
    enddo

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nvperf, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nvparf) = 0.5
    wgt(1, 1) = 0.25
    wgt(nvperf, 1) = 0.25
    wgt(1, nvparf) = 0.25
    wgt(nvperf, nvparf) = 0.25
    fgk = sum(wgt * fv) / real(np)
    return
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine halo(fv, nvperf, nvparf, fgk, np)
    implicit none
    integer, intent(in) :: nvperf, nvparf, np
    real, dimension(nvperf, nvparf), intent(inout) :: fv
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: kappa = 3.0
    real, parameter :: vthpar = 1.0, vthper = 1.0 / 2.0**0.5
    real, parameter :: vtcpar = 0.3, vtcper = 0.3
    real, parameter :: delta = 3.0 * vtcpar, p = 10.0, q = 1.0
    real, parameter :: vpermax = 5.0, vparmax = 5.0
    real, dimension(nvperf) :: vper
    real, dimension(nvparf) :: vpar
    real, dimension(nvperf, nvparf) :: fvflat, wgt
    integer :: j

    ! distribution data
    vper = linspace(0.0, vpermax, nvperf)
    vpar = linspace(-vparmax, vparmax, nvparf)
    do j = 1, nvparf
      ! factor (2 * pi * vper) is due to the cylindrical coordinate.
      fv(:, j) = (1.0 + ((vper/vthper)**2 + (vpar(j)/vthpar)**2)&
        & / (2.0 * kappa - 3.0))**(-kappa - 1.0) * (6.28318530717959 * vper)

      fvflat(:, j) = (1.0 + (((vper/vtcper)**2 + (vpar(j)/vtcpar)**2)&
        & / (2.0 * delta))**p)**(-q)
    enddo
    fv = (1.0 - fvflat) * fv

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nvperf, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nvparf) = 0.5
    wgt(1, 1) = 0.25
    wgt(nvperf, 1) = 0.25
    wgt(1, nvparf) = 0.25
    wgt(nvperf, nvparf) = 0.25
    fgk = sum(wgt * fv) / real(np)
    return
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine lembege_pellat(fxy, nxv, nyv, fgk, np)
    implicit none
    integer, intent(in) :: nxv, nyv, np
    real, dimension(nxv, nyv), intent(inout) :: fxy
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: epsln = 0.04, ascl = 2.0
    real, parameter :: Lx = 16.0, Ly = 32.0
    real, dimension(nxv) :: x
    real, dimension(nyv) :: y
    real, dimension(nxv, nyv) :: wgt
    integer :: j
    real :: Fy

    ! distribution data
    x = linspace(-Lx/2.0, Lx/2.0, nxv)
    y = linspace(-Ly, 0.0, nyv)
    do j = 1, nyv
      Fy = exp(-epsln * y(j) / ascl)
      fxy(:, j) = 1.0 / (cosh(x / (ascl * Fy)) * Fy)**2
    enddo

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nxv, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nyv) = 0.5
    wgt(1, 1) = 0.25
    wgt(nxv, 1) = 0.25
    wgt(1, nyv) = 0.25
    wgt(nxv, nyv) = 0.25
    fgk = sum(wgt * fxy) / real(np)
    return
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine butterfly(fxy, nxv, nyv, fgk, np)
    implicit none
    integer, intent(in) :: nxv, nyv, np
    real, dimension(nxv, nyv), intent(inout) :: fxy
    real, intent(inout) :: fgk

    ! local data
    real, dimension(nxv) :: x
    real, dimension(nyv) :: y
    real, dimension(nxv, nyv) :: wgt
    integer :: j

    ! distribution data
    x = linspace(-3.0, 3.0, nxv)
    y = linspace(-3.0, 3.0, nyv)
    do j = 1, nyv
      fxy(:, j) = exp(-x**2 - 2.0*y(j)**2) / cosh(5.0d0 * x * y(j)) &
        & * (x - y(j))**2
    enddo

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nxv, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nyv) = 0.5
    wgt(1, 1) = 0.25
    wgt(nxv, 1) = 0.25
    wgt(1, nyv) = 0.25
    wgt(nxv, nyv) = 0.25
    fgk = sum(wgt * fxy) / real(np)
    return
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine charged_LP_ics(fxy, nxv, nyv, fgk, np)
    implicit none
    integer, intent(in) :: nxv, nyv, np
    real, dimension(nxv, nyv), intent(inout) :: fxy
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: vDi = 1.0/3.0, vDe = -4.0/3.0
    real, parameter :: Te0 = 0.5/6.0, Ti0 = 0.5 - Te0
    real, parameter :: Tib = Ti0, Teb = Te0, nb = 0.2
    real, parameter :: epsln = 0.04
    real, parameter :: Lx = 16.0, Ly = 32.0
    real, dimension(nxv) :: x
    real, dimension(nyv) :: y
    real, dimension(nxv, nyv) :: pot, Az, wgt

    ! grid
    x = linspace(-Lx/2.0, Lx/2.0, nxv)
    y = linspace(-Ly, 0.0, nyv)
    ! get potential and vector potential
    call sol_ampere_poisson(pot, Az, x, y, nxv, nyv, &
      &vDi, Ti0, Tib, vDe, Te0, Teb, nb, epsln)

    ! distribution data
    fxy = exp((-pot + vDi * Az) / Ti0)

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nxv, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nyv) = 0.5
    wgt(1, 1) = 0.25
    wgt(nxv, 1) = 0.25
    wgt(1, nyv) = 0.25
    wgt(nxv, nyv) = 0.25
    fgk = sum(wgt * fxy) / real(np)
    return
  end subroutine


!----------------------------------------------------------------------
  subroutine charged_LP_ecs(fxy, nxv, nyv, fgk, np)
    implicit none
    integer, intent(in) :: nxv, nyv, np
    real, dimension(nxv, nyv), intent(inout) :: fxy
    real, intent(inout) :: fgk

    ! local data
    real, parameter :: vDi = 1.0/3.0, vDe = -4.0/3.0
    real, parameter :: Te0 = 0.5/6.0, Ti0 = 0.5 - Te0
    real, parameter :: Tib = Ti0, Teb = Te0, nb = 0.2
    real, parameter :: epsln = 0.04
    real, parameter :: Lx = 16.0, Ly = 32.0
    real, dimension(nxv) :: x
    real, dimension(nyv) :: y
    real, dimension(nxv, nyv) :: pot, Az, wgt

    ! grid
    x = linspace(-Lx/2.0, Lx/2.0, nxv)
    y = linspace(-Ly, 0.0, nyv)
    ! get potential and vector potential
    call sol_ampere_poisson(pot, Az, x, y, nxv, nyv, &
      &vDi, Ti0, Tib, vDe, Te0, Teb, nb, epsln)

    ! distribution data
    fxy = exp((pot - vDe * Az) / Te0)

    ! density normalization factor
    wgt = 1.0
    wgt(1, :) = 0.5
    wgt(nxv, :) = 0.5
    wgt(:, 1) = 0.5
    wgt(:, nyv) = 0.5
    wgt(1, 1) = 0.25
    wgt(nxv, 1) = 0.25
    wgt(1, nyv) = 0.25
    wgt(nxv, nyv) = 0.25
    fgk = sum(wgt * fxy) / real(np)
    return
  end subroutine


!----------------------------------------------------------------------
  subroutine sol_ampere_poisson(pot, Az, x, y, nxv, nyv, &
      &vDi, Ti0, Tib, vDe, Te0, Teb, nb, epsln)
    implicit none
    integer, intent(in) :: nxv, nyv
    real, dimension(nxv, nyv), intent(inout) :: pot, Az
    real, dimension(nxv), intent(in) :: x
    real, dimension(nyv), intent(in) :: y
    real, intent(in) :: vDi, Ti0, Tib, vDe, Te0, Teb, nb, epsln

    ! local data
    real, dimension(nyv) :: pot0, Az0
    real, dimension(3, (nxv-1)/2+1) :: u
    real, dimension(3) :: u0
    integer :: j, nxhd
    integer :: ierr

    ! Az at x = 0
    Az0 = epsln * y
    ! pot at x = 0
    call get_bnd_pot(pot0, Az0, nyv, ierr, vDi, Ti0, Tib, vDe, Te0, Teb, nb)
    if (ierr == 1) write(*, *) 'rho_deriv is too small'
    if (ierr == 2) write(*, *) 'Newton method did not converge'
    ! write(*, *) (pot0(j), j = 1, nyv)

    nxhd = size(u, 2)
    do j = 1, nyv
      ! initial data
      u0 = (/Az0(j), 0.0, pot0(j)/)
      ! solve IVP
      call RK4(u, u0, x, nxhd, vDi, Ti0, Tib, vDe, Te0, Teb, nb)
      ! store data
      Az(nxhd:, j) = u(1, :)
      pot(nxhd:, j) = u(3, :)
      Az(:nxhd-1, j) = u(1, nxhd:2:-1)
      pot(:nxhd-1, j) = u(3, nxhd:2:-1)
    enddo
  end subroutine


!----------------------------------------------------------------------
  subroutine RK4(u, u0, xv, nxhd, vDi, Ti0, Tib, vDe, Te0, Teb, nb)
    implicit none
    integer, intent(in) :: nxhd
    real, dimension(3, nxhd), intent(inout) :: u
    real, dimension(3), intent(in) :: u0
    real, dimension(nxhd), intent(in) :: xv
    real, intent(in) :: vDi, Ti0, Tib, vDe, Te0, Teb, nb

    ! local data
    real :: dx
    real, dimension(3) :: F0, F1, F2, F3
    integer :: j

    ! step size
    dx = xv(2) - xv(1)
    ! initial data
    u(:, 1) = u0

    do j = 1, nxhd-1
      F0 = ampere_poisson_deriv(u(:, j),&
        &vDi, Ti0, Tib, vDe, Te0, Teb, nb)
      F1 = ampere_poisson_deriv(u(:, j) + 0.5 * dx * F0,&
        &vDi, Ti0, Tib, vDe, Te0, Teb, nb) 
      F2 = ampere_poisson_deriv(u(:, j) + 0.5 * dx * F1,&
        &vDi, Ti0, Tib, vDe, Te0, Teb, nb)
      F3 = ampere_poisson_deriv(u(:, j) + dx * F2,&
        &vDi, Ti0, Tib, vDe, Te0, Teb, nb)
      u(:, j+1) = u(:, j) + dx * (F0 + 2.0 * F1 + 2.0 * F2 + F3) / 6.0
    enddo
    return
  end subroutine


!----------------------------------------------------------------------
  function ampere_poisson_deriv(u, vDi, Ti0, Tib, vDe, Te0, Teb, nb) result(dudx)
    ! Az = u[1]; By = u[2]; pot = u[3]
    implicit none
    ! output
    real, dimension(3) :: dudx
    ! input
    real, dimension(3), intent(in) :: u
    real, intent(in) :: vDi, Ti0, Tib, vDe, Te0, Teb, nb

    ! local data
    real :: exp_i0, exp_e0, exp_ib, exp_eb

    exp_i0 = exp((-u(3) + vDi * u(1)) / Ti0)
    exp_e0 = exp((u(3) - vDe * u(1)) / Te0)
    exp_ib = exp(-u(3) / Tib)
    exp_eb = exp(u(3) / Teb)

    ! d/dx (Az)
    dudx(1) = -u(2)
    ! d/dx (By)
    dudx(2) = vDi * exp_i0 - vDe * exp_e0
    ! d/dx (pot)
    dudx(3) = -u(2) * (vDi / Ti0 * exp_i0 + vDe / Te0 * exp_e0) &
      &/ (exp_i0 / Ti0 + nb * exp_ib / Tib + exp_e0 / Te0 + nb * exp_eb / Teb)
    return
  end function


!----------------------------------------------------------------------
  subroutine get_bnd_pot(pot, Az, nyv, ierr, vDi, Ti0, Tib, vDe, Te0, Teb, nb)
    ! find potential at the boundary using Newton's method
    implicit none
    integer, intent(in) :: nyv
    real, dimension(nyv), intent(inout) :: pot
    real, dimension(nyv), intent(in) :: Az
    integer, intent(inout) :: ierr
    real, intent(in) :: vDi, Ti0, Tib, vDe, Te0, Teb, nb

    ! local data
    integer, parameter :: maxiter = 20
    integer :: j, k
    real :: potg, pot_new, rho, rho_deriv
    logical :: sol

    ierr = 0
    do j = 1, nyv
      ! initial guess
      potg = 0.0
      sol = .false.
      do k = 1, maxiter
        rho = charge_density(potg, Az(j), vDi, Ti0, Tib, vDe, Te0, Teb, nb)
        rho_deriv = charge_density_deriv(potg, Az(j), vDi, Ti0, Tib,&
          & vDe, Te0, Teb, nb)
        ! check if denominator is too small
        if (abs(rho_deriv) < 1e-15) then
          ierr = 1
          return
        endif
        ! update pot_new
        pot_new = potg - rho / rho_deriv
        ! break loop if solution reaches convergence
        if (abs(pot_new - potg) < 1e-3) then
          sol = .true.
          exit
        endif
        ! update potg to iterate
        potg = pot_new
      enddo
      if (sol) then
        ! store data
        pot(j) = pot_new
      else
        ! Newton method did not converge
        ierr = 2
        return
      endif
    enddo
    return
  end subroutine


!----------------------------------------------------------------------
  function charge_density(pot, Az, vDi, Ti0, Tib, vDe, Te0, Teb, nb) result(rho)
    ! total charge density in Lembege-Pellat current sheet
    implicit none
    ! output
    real :: rho
    ! input
    real, intent(in) :: pot, Az, vDi, Ti0, Tib, vDe, Te0, Teb, nb

    rho = exp((-pot + vDi * Az) / Ti0) + nb * exp(-pot / Tib)&
      & - exp((pot - vDe * Az) / Te0) - nb * exp(pot / Teb)
    return
  end function
  

!----------------------------------------------------------------------
  function charge_density_deriv(pot, Az, vDi, Ti0, Tib, vDe, Te0, Teb, nb) &
      &result (rho_deriv)
    ! derivative of total charge density with respect to potential
    implicit none
    ! output
    real :: rho_deriv
    ! input
    real, intent(in) :: pot, Az, vDi, Ti0, Tib, vDe, Te0, Teb, nb

    rho_deriv = -exp((-pot + vDi * Az) / Ti0) / Ti0&
      & - nb * exp(-pot / Tib) / Tib&
      & - exp((pot - vDe * Az) / Te0) / Te0&
      & - nb * exp(pot / Teb) / Teb
    return
  end function
!
!
!----------------------------------------------------------------------
  function linspace(start, end, num) result(samples)
    implicit none
    ! input
    real, intent(in) :: start, end
    integer, intent(in) :: num
    ! output
    real, dimension(num) :: samples
    ! local data
    real :: step
    integer :: j

    step = (end - start) / real(num - 1)

    samples = (/ (start + (j - 1) * step, j = 1, num) /)
  end function
!
!
end module
