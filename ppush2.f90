module ppush2
  use invsampling
  implicit none

  contains


  subroutine pdcomp2(edges, nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps)
    implicit none
    integer, intent(in) :: ny, kstrt, nvp, idps
    real, dimension(idps), intent(inout) :: edges
    integer, intent(inout) :: nyp, noff, nypmx, nypmn

    ! local data
    integer :: kb, kr
    real :: at1
    integer, dimension(2) :: mypm, iwork2

    kb = kstrt - 1
    at1 = real(ny) / real(nvp)
    edges(1) = at1 * real(kb)
    noff = edges(1) + 0.5
    edges(2) = at1 * real(kb + 1)
    kr = edges(2) + 0.5
    nyp = kr - noff
    ! find maximum/minimum partition size
    mypm(1) = nyp
    mypm(2) = -nyp
    call PPIMAX(mypm, iwork2, 2)
    nypmx = mypm(1) + 1
    nypmn = -mypm(2)
    return
  end subroutine

  subroutine pfedges2(edges, nyp, noff, fmargy, nypmx, nypmn, ny,&
      &kstrt, nvp, idps)
    implicit none
    integer, intent(in) :: ny, kstrt, nvp, idps
    real, dimension(idps), intent(inout) :: edges
    integer, intent(inout) :: nyp, noff, nypmx, nypmn
    real, dimension(ny+1), intent(in) :: fmargy

    ! local data
    integer :: kb, ny1, ncutoff, ierr
    real, dimension(2) :: anplr, at2
    real, dimension(ny+1) :: cdfy
    real, dimension(262145) :: coeffy
    real, parameter :: tol = 1e-12, eps = 1e-15
    integer, dimension(2) :: mypm, iwork2

    ny1 = ny + 1

    ! particle distribution constants
    kb = kstrt - 1
    ! divide interval [0, 1] into nvp intervals.
    ! the interval on this processor:
    anplr(1) = real(kb) / real(nvp)
    anplr(2) = real(kb + 1) / real(nvp)
    ! cumulative sum
    cdfy = my_cumsum(fmargy, ny1)
    ! normalization
    cdfy = cdfy / cdfy(ny1)
    ! chebyshev coefficients of cdfy
    call cheb_coeff(cdfy, coeffy, ny1, ncutoff, tol, ierr)
    ! error handling
    if (ierr == 1) return
    ! inversion
    at2 = inv_bisection(coeffy, anplr, ncutoff, 2, eps)
    ! mapping
    call map_unit_to_grid(at2, edges, 2, ny1)
    ! set leftmost edge to zero
    if (kb == 0) edges(1) = 0.0
    ! set rightmost edge to ny
    if ((kb + 1) == nvp) edges(2) = real(ny)
    ! calculate number of grids and offsets in new partitions
    noff = edges(1) + 0.5
    kb = edges(2) + 0.5
    nyp = kb - noff
    edges(1) = real(noff)
    edges(2) = real(kb)

    ! find maximum/minimum partition size
    mypm(1) = nyp
    mypm(2) = -nyp
    call PPIMAX(mypm, iwork2, 2)
    nypmx = mypm(1) + 1
    nypmn = -mypm(2)
    return
  end subroutine


  subroutine ppgpost2l(part, q, npp, noff, qm, idimp, npmax, nxv, nypmx)
    implicit none
    integer, intent(in) :: npp, noff, idimp, npmax, nxv, nypmx
    real, dimension(idimp, npmax), intent(in) :: part
    real, dimension(nxv, nypmx), intent(inout) :: q
    real, intent(in) :: qm

    ! local data
    integer :: mnoff, j, nn, np, mm, mp
    real :: dxp, dyp, amx, amy
    mnoff = noff - 1

    do j = 1, npp
      ! find interpolation weights
      nn = part(1, j)
      mm = part(2, j)
      dxp = qm * (part(1, j) - real(nn))
      dyp = part(2, j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
      ! deposit charge
      q(np, mp) = q(np, mp) + dxp * dyp
      q(nn, mp) = q(nn, mp) + amx * dyp
      q(np, mm) = q(np, mm) + dxp * amy
      q(nn, mm) = q(nn, mm) + amx * amy
    enddo

    return
  end subroutine

end module
