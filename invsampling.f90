module invsampling
  use mfft1
  ! use util
  implicit none

  contains


!----------------------------------------------------------------------
! This subroutine distributes particles based on the edges array.
  subroutine pp_inv_sampling_2D(part, npp, fxy, npx, npy, nxv, nyv, edges,&
      &idimp, npmax, idps, ierr, tol, eps, tcheb, tbisc)
    implicit none
    integer, intent(in) :: npx, npy, nxv, nyv, idimp, npmax, idps
    integer, intent(inout) :: npp, ierr
    real, dimension(idimp, npmax), intent(inout) :: part
    real, dimension(nxv, nyv), intent(in) :: fxy
    real, dimension(idps), intent(in) :: edges
    real, intent(in) :: tol, eps
    real, intent(inout) :: tcheb, tbisc

    ! local data
    integer :: kl, kr, npt, k, ly, lyr, joff, nps
    real :: dy
    integer, dimension(1) :: ierr1, iwork1
    double precision, dimension(1) :: sum1, work1
    real, dimension(nyv) :: margy
    real, dimension(npy) :: party
    real, dimension(nxv) :: condx
    real, dimension(npx) :: partx
    double precision :: dnpxy

    nps = npp

    ! marginal distribution in y direction
    margy = (fxy(1, :) + fxy(nxv, :)) * 0.5 + sum(fxy(2:nxv-1, :), 1)
    ! sample margy
    call inv_sampling_1D(margy, party, nyv, npy, tol, eps, tcheb, tbisc, ierr)
    ! error handling
    if (ierr == 1) return

    ! find particle partition on this processor
    kl = minloc(abs(party - edges(1)), 1)
    if (party(kl) <= edges(1)) kl = kl + 1
    kr = minloc(abs(party - edges(2)), 1)
    if (party(kr) > edges(2)) kr = kr - 1
    ! check if particle partition is correct
    ierr = kr - kl
    if (kr < kl) return
    ierr = 0
    ! check if particle overflow will occur
    npt = npp + npx * (kr - kl + 1)
    ! This chunk of code will cause problem if the execution return to the main
    ! program on cerntain processors. Comment out for now.
    ! ierr1(1) = npt
    ! call PPIMAX(ierr1, iwork1, 1)
    ! ierr = ierr1(1)
    ! if (ierr > npmax) return
    ! ierr = 0

    do k = kl, kr
      ! conditional distribution in x direction for a given y
      ly = party(k) + 1
      dy = party(k) - ly + 1
      lyr = ly + 1
      condx = fxy(:, ly) * (1.0 - dy) + fxy(:, lyr) * dy
      ! sample condx
      call inv_sampling_1D(condx, partx, nxv, npx, tol, eps, tcheb, tbisc, ierr)
      ! error handling
      if (ierr == 1) return
      ! store co-ordinates
      joff = npx * (k - kl) + npp
      part(1, joff+1:joff+npx) = partx
      part(2, joff+1:joff+npx) = party(k)
    enddo

    ! update number of particles
    npp = npt
    ! This chunk of code will cause problem if the execution return to the main
    ! program on cerntain processors. Comment out for now.
    ! check if not all particles were distributed
    ! dnpxy = dble(npx) * dble(npy)
    ! sum1(1) = dble(npp - nps)
    ! call PPDSUM(sum1, work1, 1)
    ! dnpxy = sum1(1) - dnpxy
    ! if (dnpxy /= 0.0d0) ierr = -1
    return
  end subroutine


!----------------------------------------------------------------------
! This subroutine distributes particles uniformly across processors.
  ! subroutine pp_inv_sampling_2D(part, npp, fxy, npx, npy, nxv, nyv, kstrt, nvp,&
  !     &idimp, npmax, ierr, tol, eps, tcheb, tbisc)
  !   implicit none
  !   integer, intent(in) :: npx, npy, nxv, nyv, kstrt, nvp, idimp, npmax
  !   integer, intent(inout) :: npp, ierr
  !   real, dimension(idimp, npmax), intent(inout) :: part
  !   real, dimension(nxv, nyv), intent(in) :: fxy
  !   real, intent(in) :: tol, eps
  !   real, intent(inout) :: tcheb, tbisc

  !   ! local data
  !   integer :: mpy, mpys, ks, npt, moff, kk, k, ly, lyr, joff, nps
  !   real :: dy
  !   integer, dimension(1) :: ierr1, iwork1
  !   double precision, dimension(1) :: sum1, work1
  !   real, dimension(nyv) :: margy
  !   real, dimension(npy) :: party
  !   real, dimension(nxv) :: condx
  !   real, dimension(npx) :: partx
  !   double precision :: dnpxy

  !   nps = npp
  !   ! particle distribution constants
  !   ks = kstrt - 1
  !   ! mpy = number of particles per processor in y direction
  !   mpy = (npy - 1) / nvp + 1
  !   mpys = min(mpy, max(0, npy - mpy*ks))
  !   ! check if particle overflow will occur
  !   npt = npp + npx * mpys
  !   ierr1(1) = npt
  !   call PPIMAX(ierr1, iwork1, 1)
  !   ierr = ierr1(1)
  !   if (ierr > npmax) return
  !   ierr = 0

  !   ! marginal distribution in y direction
  !   margy = (fxy(1, :) + fxy(nxv, :)) * 0.5 + sum(fxy(2:nxv-1, :), 1)
  !   ! sample margy
  !   call inv_sampling_1D(margy, party, nyv, npy, tol, eps, tcheb, tbisc, ierr)
  !   ! error handling
  !   if (ierr == 1) return

  !   moff = mpy * ks
  !   do kk = 1, mpys
  !     ! conditional distribution in x direction for a given y
  !     k = moff + kk
  !     ly = party(k) + 1
  !     dy = party(k) - ly + 1
  !     lyr = ly + 1
  !     condx = fxy(:, ly) * (1.0 - dy) + fxy(:, lyr) * dy
  !     ! sample condx
  !     call inv_sampling_1D(condx, partx, nxv, npx, tol, eps, tcheb, tbisc, ierr)
  !     ! error handling
  !     if (ierr == 1) return
  !     ! store co-ordinates
  !     joff = npx * (kk-1) + npp
  !     part(1, joff+1:joff+npx) = partx
  !     part(2, joff+1:joff+npx) = party(k)
  !   enddo

  !   ! update number of particles
  !   npp = npt
  !   ! check if not all particles were distributed
  !   dnpxy = dble(npx) * dble(npy)
  !   sum1(1) = dble(npp - nps)
  !   call PPDSUM(sum1, work1, 1)
  !   dnpxy = sum1(1) - dnpxy
  !   if (dnpxy /= 0.0d0) ierr = -1
  !   return
  ! end subroutine


!----------------------------------------------------------------------
  subroutine inv_sampling_2D(fdata, xsamples, ysamples, nxgrid, nygrid,&
      &nxsample, nysample, tol, eps, tcheb, tbisc, ierr)
    implicit none
    ! fdata: f(x) data on two dimensional grid.
    ! xsamples: position of samples in the x direction; 2D array.
    ! ysamples: position of samples in the y direction; 1D array.
    ! (nxgrid, nygrid): dimension of fdata.
    ! (nxsample, nysample): dimension of xsamples.
    ! nysample: dimension of ysamples.
    ! tol: tolerance in constructing chebyshev coefficients.
    ! eps: error control in bisection method.
    ! tcheb: time spent on calculating chebyshev coefficients.
    ! tbisc: time spent on root-finding with bisection method.
    ! ierr: ierr = 1 success; ierr = 0 failure.
    integer, intent(in) :: nxgrid, nygrid, nxsample, nysample
    real, dimension(nxgrid, nygrid), intent(in) :: fdata
    real, dimension(nxsample, nysample), intent(inout) :: xsamples
    real, dimension(nysample), intent(inout) :: ysamples
    real, intent(in) :: tol, eps
    real, intent(inout) :: tcheb, tbisc
    integer, intent(inout) :: ierr

    ! local data
    real, dimension(nygrid) :: margy
    real, dimension(nxgrid) :: condx
    integer :: j, ly, lyr
    real :: dy

    ! marginal distribution in y
    margy = (fdata(1, :) + fdata(nxgrid, :)) * 0.5 + sum(fdata(2:nxgrid-1, :), 1)
    ! sample margy
    call inv_sampling_1D(margy, ysamples, nygrid, nysample,&
      &tol, eps, tcheb, tbisc, ierr)
    ! error handling
    if (ierr == 1) return

    do j = 1, nysample
      ! conditional distribution in x for a given y
      ly = ysamples(j) + 1
      dy = ysamples(j) - ly + 1
      lyr = ly + 1
      condx = fdata(:, ly) * (1.0 - dy) + fdata(:, lyr) * dy
      ! sample condx
      call inv_sampling_1D(condx, xsamples(:, j), nxgrid, nxsample,&
        &tol, eps, tcheb, tbisc, ierr)
      ! error handling
      if (ierr == 1) return
    enddo
  end subroutine

!----------------------------------------------------------------------
  subroutine inv_sampling_1D(fdata, xsamples, ngrid, nsamples,&
      &tol, eps, tcheb, tbisc, ierr)
    implicit none
    ! fdata: original grid data of f(x).
    ! xsamples: position of samples in the grid coordinate; array.
    ! ngrid: number of original grids.
    ! nsamples: number of samples.
    ! tol: tolerance in constructing chebyshev coefficients.
    ! eps: error control in bisection method.
    ! tcheb: time spent on calculating chebyshev coefficients.
    ! tbisc: time spent on root-finding with bisection method.
    ! ierr: ierr = 1 success; ierr = 0 failure.
    integer, intent(in) :: ngrid, nsamples
    real, dimension(ngrid), intent(in) :: fdata
    real, dimension(nsamples), intent(inout) :: xsamples
    real, intent(in) :: tol, eps
    real, intent(inout) :: tcheb, tbisc
    integer, intent(inout) :: ierr

    ! local data
    real, dimension(ngrid) :: cdfx
    real, dimension(262145) :: coeffx
    integer :: ncutoff
    integer :: j
    real, dimension(nsamples) :: y, xtmp
    integer, dimension(4) :: itime
    double precision :: dtime

    ! cumulative sum
    cdfx = my_cumsum(fdata, ngrid)
    ! normalization
    cdfx = cdfx / cdfx(ngrid)
    ! initialize timer
    call dtimer(dtime, itime, -1)
    ! chebyshev coefficients of cdfx
    call cheb_coeff(cdfx, coeffx, ngrid, ncutoff, tol, ierr)
    ! record time
    call dtimer(dtime, itime, 1)
    tcheb = tcheb + real(dtime)
    ! error handling
    if (ierr == 1) return

    ! initialize timer
    call dtimer(dtime, itime, -1)
    ! inverse transform sampling
    y = (/ ((j - 0.5d0) / real(nsamples), j = 1, nsamples) /)
    xsamples = inv_bisection(coeffx, y, ncutoff, nsamples, eps)
    ! record time
    call dtimer(dtime, itime, 1)
    tbisc = tbisc + real(dtime)

    ! map x back to grid coordinate
    call map_unit_to_grid(xsamples, xtmp, nsamples, ngrid)
    xsamples = xtmp
  end subroutine
!
!
!----------------------------------------------------------------------
  function my_cumsum(pdf, ngrid) result(cdf)
    implicit none
    ! input
    integer, intent(in) :: ngrid
    real, dimension(ngrid), intent(in) :: pdf
    ! output
    real, dimension(ngrid) :: cdf

    ! local data
    integer :: j

    cdf(1) = 0.0
    do j = 2, ngrid
      cdf(j) = cdf(j-1) + (pdf(j-1) + pdf(j)) / 2.0
    enddo
  end function
!
!
!----------------------------------------------------------------------
  function inv_bisection(cdf, y, nc, ny, eps) result(x)
    implicit none
    ! this function finds x = cdf^{-1}(y) using bisection method.
    ! cdf: chebyshev coefficients of cumulative distribution function.
    ! y: random numbers on [0, 1]; array.
    ! nc: effective length of chebyshev coefficients cdf.
    ! ny: number of random numbers y.
    ! eps: error control. 
    ! input:
    integer, intent(in) :: nc, ny
    real, dimension(nc), intent(in) :: cdf
    real, dimension(ny), intent(in) :: y
    real, intent(in) :: eps
    ! output:
    real, dimension(ny) :: x

    ! local data
    real, dimension(ny) :: a, b, c, vals, zeros
    logical, dimension(ny) :: I1, I2, I3

    a = -1.0
    b = 1.0
    c = 0.0
    zeros = 0.0

    do while (maxval(abs(b - a), 1) >= eps)
      vals = cheb_clenshaw(c, cdf, ny, nc) - y
      I1 = (vals <= -eps)
      I2 = (vals >= eps)
      I3 = ((.not. I1) .and. (.not. I2))
      a = merge(c, zeros, I1) + merge(a, zeros, I2) + merge(c, zeros, I3)
      b = merge(b, zeros, I1) + merge(c, zeros, I2) + merge(c, zeros, I3)
      c = (a + b) / 2.0
    enddo

    x = c
  end function
!
!
!----------------------------------------------------------------------
  function cheb_clenshaw(x, c, nx, nc) result(sum_cheb)
    implicit none
    ! clenshaw algorithm to perform the sum:
    ! output = c_0 * T_0(x) + c_1 * T_1(x) + ... + c_n * T_n(x).
    ! x: independent variable; array.
    ! c: coefficient of the chebyshev sum.
    ! nx: number of x points to be evaluated.
    ! nc: effective length of chebyshev coefficients.
    ! input
    integer, intent(in) :: nx, nc
    real, dimension(nx), intent(in) :: x
    real, dimension(nc), intent(in) :: c
    ! output
    real, dimension(nx) :: sum_cheb

    ! local data
    real, dimension(nx) :: bk1, bk2, xx, bk
    integer :: k

    bk1 = 0.0
    bk2 = 0.0
    xx = 2.0 * x

    do k = nc, 1, -1
      bk = c(k) + xx * bk1 - bk2
      bk2 = bk1
      bk1 = bk
    enddo

    sum_cheb = bk1 - x * bk2
  end function
!
!
!----------------------------------------------------------------------
  subroutine cheb_coeff(fdata, coeffs, nf, ncutoff, tol, ierr)
    implicit none
    ! fdata: original data of f(x).
    ! coeffs: chebyshev coefficients (c_0, c_1, ..., c_n) of f(x).
    ! f(x) = c_0 * T_0(x) + c_1 * T_1(x) + ... + c_n * T_n(x).
    ! c_k = inverse discrete chebyshev transform of f(x).
    ! ncutoff: the number of non-negligible chebyshev coefficients.
    ! tol: target relative accuracy used in standard chop.
    ! ierr: ierr = 1, construction failed; ierr = 0, construction succeeded.
    integer, intent(in) :: nf
    real, dimension(nf), intent(in) :: fdata
    real, dimension(262145), intent(inout) :: coeffs
    integer, intent(inout) :: ncutoff
    real, intent(in) :: tol
    integer, intent(inout) :: ierr
  
    ! local data
    integer :: j
    integer :: n
    real, dimension(262145) :: x, mapx
  
    ! nf = size(fdata, 1)
  
    do j = 4, 18
      ! n: number of chebyshev points.
      ! the degree of chebyshev polynomial is n-1.
      n = 2**j + 1
      ! chebyshev points.
      call chebpts(x, n)
      ! map chebyshev points to original grid.
      call map_unit_to_grid(x, mapx, n, nf)
      ! evaluate f(x) at chebyshev points.
      call interp1d(mapx, fdata, coeffs, n, nf)
      ! discrete chebyshev transform.
      call DCT(coeffs, n-1, j+1)
      ! truncate negligible chebyshev coefficients.
      call cheb_standardChop(coeffs, n, ncutoff, tol)
      ! ncutoff represents the number of non-negligible chebyshev coefficients.
      ! we are happy if ncutoff < n.
      if (ncutoff < n) exit
    enddo
  
    if (n == (2**18 + 1)) then
      ierr = 1
    else
      ierr = 0
    endif

    return
  end subroutine
!  
!  
!----------------------------------------------------------------------
  subroutine chebpts(x, n)
    implicit none
    ! x: chebyshev points to be generated.
    ! n: number of chebyshev points.
    integer, intent(in) :: n
    real, dimension(:), intent(inout) :: x
  
    ! local data
    integer :: j
  
    x(:n) = (/ (j - 1, j = 1, n) /)
    x(:n) = cos(x(:n) / real(n - 1) * 3.14159265358979)
  end subroutine
!  
!  
!----------------------------------------------------------------------
  subroutine map_unit_to_grid(xunit, xgrid, nx, ngrid)
    ! map xunit on [-1, 1] to grid coordinate [0, ngrid-1]
    implicit none
    ! xunit: chebyshev points on the interval [-1, 1].
    ! xgrid: points on the grid coordinate [0, ngrid-1].
    ! nx: number of data points.
    ! ngrid: number of original grids of f(x).
    integer, intent(in) :: nx, ngrid
    real, dimension(nx), intent(in) :: xunit
    real, dimension(ngrid), intent(inout) :: xgrid

    xgrid(:nx) = (xunit(:nx) + 1.0) / 2.0 * real(ngrid - 1)
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine map_grid_to_unit(xgrid, xunit, ngrid)
    ! map xgrid on grid coordinate [0, ngrid-1] to [-1, 1].
    implicit none
    ! xgrid: points on the grid coordinate [0, ngrid-1].
    ! xunit: points on the unit interval [-1, 1].
    ! ngrid: number of original grids.
    integer, intent(in) :: ngrid
    real, dimension(ngrid), intent(in) :: xgrid
    real, dimension(ngrid), intent(inout) :: xunit

    xunit = xgrid / real(ngrid - 1) * 2.0 - 1.0
  end subroutine
!  
!  
!----------------------------------------------------------------------
  subroutine interp1d(xa, fx, fxa, n, ngrid)
    implicit none
    ! xa: abscissa values for output.
    ! fx: function values at original grids.
    ! fxa: function values at xa.
    ! n: number of chebyshev points.
    ! ngrid: number of original grids of f(x).
    integer, intent(in) :: n, ngrid
    real, dimension(n), intent(in) :: xa
    real, dimension(ngrid), intent(in) :: fx
    real, dimension(n), intent(inout) :: fxa
  
    ! local data
    integer :: j
    integer :: lx, lxr
    real :: dx
    real, dimension(ngrid+1) :: fx1
  
    ! special case: xa = ngrid - 1. lxr exceeds ngrid.
    fx1(1:ngrid) = fx
    fx1(ngrid+1) = 0.0
  
    do j = 1, n 
      lx = xa(j) + 1
      dx = xa(j) - lx + 1
      lxr = lx + 1
      fxa(j) = fx1(lx) * (1.0 - dx) + fx1(lxr) * dx
    enddo
  end subroutine
!  
!  
!----------------------------------------------------------------------
  subroutine DCT(fx, n, indx)
    implicit none
    ! fx input: f(x) evaluated at the chebyshev points.
    ! fx output: chebyshev coefficients.
    ! n: degree of chebyshev polynomials.
    ! n + 1: number of chebyshev points.
    ! indx: 2**indx = 2*n ---> indx = log2(n) + 1
    integer, intent(in) :: n, indx
    real, dimension(n+1), intent(inout) :: fx
  
    ! local data
    real, dimension(2*n) :: fxx
    integer, dimension(n) :: mixup
    complex, dimension(n) :: sct
    integer :: isign
    ! Note: tfft will not be returned to the higher level subroutine
    ! because we use a separate timer in the higher level.
    ! keep tfft for compatibility with the mfft1 module.
    real :: tfft = 0.0
  
    ! extend fx from [0, 1, ..., n] to [0, 1, ..., 2n-1]
    fxx(:(n+1)) = fx(:(n+1))
    fxx((n+2):) = fx(n:2:-1)
  
    ! prepare fft tables: update mixup, sct.
    call mfft1_init(mixup, sct, indx)
    ! inverse FFT
    isign = -1
    call mfft1r(fxx, isign, mixup, sct, tfft, indx)
    ! take real part of the coefficients
    ! and scale the interior points: modes [1, 2, 3, ..., n-1].
    ! mode 0.
    fx(1) = fxx(1)
    ! mode n.
    fx(n+1) = fxx(2)
    ! modes [1, 2, 3, ..., n-1].
    fx(2:n) = 2.0 * fxx(3::2)
  end subroutine
!
!
!----------------------------------------------------------------------
  subroutine cheb_standardChop(coeffs, n, ncutoff, tol)
    implicit none
    ! This subroutine is a direct translation of standardChop in
    ! the Chebfun software system.
    ! coeffs: chebyshev coefficients.
    ! n: number of input chebyshev coefficients.
    ! ncutoff: number of output chebyshev coefficients.
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: coeffs
    integer, intent(inout) :: ncutoff
    real, intent(in) :: tol

    ! local data
    real, dimension(n) :: b, m, envelope
    integer :: j
    integer :: plateauPoint, j2
    real :: e1, e2, r
    logical :: plateau
    integer :: j3, d
    real, dimension(:), allocatable :: cc

    ncutoff = n

    ! step 1: convert coeffs to a nonincreasing envelope normalized
    ! to begin with value 1.
    b = abs(coeffs(:n))
    m = b(n)
    do j = n-1, 1, -1
      m(j) = max(b(j), m(j+1))
    enddo
    if (m(1) == 0) then
      ncutoff = 1
      return
    endif
    envelope = m / m(1)

    ! step 2: find a plateauPoint, if any.
    plateauPoint = n
    j2 = n
    do j = 2, n
      j2 = nint(1.25 * j + 5)
      if (j2 > n) return
      e1 = envelope(j)
      e2 = envelope(j2)
      r = 3.0 * (1.0 - log(e1) / log(tol))
      plateau = (e1 == 0) .or. (e2 / e1 > r)
      if (plateau) then
        plateauPoint = j - 1
        exit
      endif
    enddo

    ! step 3
    if (envelope(plateauPoint) == 0) then
      ncutoff = plateauPoint
    else
      j3 = count(envelope >= tol**(7.0/6.0))
      if (j3 < j2) then
        j2 = j3 + 1
        envelope(j2) = tol**(7.0/6.0)
      endif
      allocate(cc(j2))
      cc = log10(envelope(:j2))
      ! cc = cc + linspace(0.0, (-1.0/3.0) * log10(tol), j2)
      cc = cc + (/ (-1.0/3.0*log10(tol) * (j - 1) / real(j2 - 1), j = 1, j2) /)
      d = minloc(cc, 1)
      ncutoff = max(d - 1, 1)
    endif
  end subroutine
!
end module
