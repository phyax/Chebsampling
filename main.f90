program main
  use invsampling
  use pplib2
  use ppush2
  use distr
  implicit none

  ! generate samples
  call gen_samples()

  contains


  subroutine gen_samples()
    implicit none

    ! icase = identiy of the example distribution data
    integer, dimension(1) :: icase

    ! number of processors (nvp = 1 is handled.)
    integer :: nvp = 64
    integer :: idproc, kstrt

    ! number of grids and auxillary variables
    integer :: ncx, ncy, ncx1, ncy1

    ! npx, npy = number of particle in x and y direction
    ! idimp = number of particle coordinates
    integer :: npx, npy
    integer, parameter :: idimp = 2

    ! fxy = pre-defined distribution data
    ! fgk = normalization factor
    ! real, dimension(ncx1, ncy1) :: fxy
    real, dimension(:, :), allocatable :: fxy
    real :: fgk

    ! density = particle density deposited on grids
    ! scrg = scratch array to collect density on different processors
    ! scr = scratch array to add density from guard cells
    real, dimension(:, :), allocatable :: density, scrg
    real, dimension(:), allocatable :: scr

    ! part = array for particle coordinates
    ! np = total number of particles
    ! npp = number of particles in each partition
    ! npmax = maximum number of particles in each partition
    real, dimension(:, :), allocatable :: part
    integer :: np, npp, npmax

    ! tol = relative tolerance in constructing chebyshev polynomials.
    ! eps = relative error in the bisection method.
    ! tcheb = time diagnostic for constructing chebyshev polynomials.
    ! tbisc = time diagnostic for root finding with bisection method.
    ! ierr = error flag.
    real, parameter :: tol = 1e-8, eps = 1e-14
    real :: tcheb = 0.0, tbisc = 0.0
    integer :: ierr

    ! edges(1:2) = lower:upper boundary of particle partition.
    ! nyp = number of primary (complete) gridpoints in particle partition.
    ! noff = lowermost global gridpoint in particle partition.
    ! nypmx = maximum size of particle partition, including guard cells.
    ! nypmn = minimum value of nyp.
    ! fmargy = marginal distribution in y direction.
    integer, parameter :: idps = 2
    real, dimension(idps) :: edges
    integer :: nyp, noff, nypmx, nypmn
    ! real, dimension(ncy1) :: fmargy
    real, dimension(:), allocatable :: fmargy
    
    integer :: j, k

    namelist /grids/ ncx, ncy
    namelist /particles/ npx, npy

    open(unit=1, file='input1.nml', status='old')
    read(1, grids)
    read(1, particles)
    close(unit=1)

    ! allocate model distribution array
    ncx1 = ncx + 1
    ncy1 = ncy + 1
    allocate(fxy(ncx1, ncy1), fmargy(ncy1))

    call PPINIT2(idproc, nvp)
    kstrt = idproc + 1

    ! check if too many processors
    if (nvp > ncy) then
      if (kstrt == 1) write(*,*) 'Too many nodes requested: ncy, nvp=', ncy, nvp
      call PPEXIT()
      stop
    endif

    ! allocate particle array
    np = dble(npx) * dble(npy)
    ! npmax = (np / nvp) * 1.5
    npmax = np
    allocate(part(idimp, npmax))

    if (kstrt == 1) then
      write(*, *) 'Input ID for the example distribution data:'
      write(*, *) '0 = Butterfly distribution;'
      write(*, *) '1 = Nonpolarized Lembege-Pellat current sheet;'
      write(*, *) '2 = Polarized Lembege-Pellat current sheet for ions;'
      write(*, *) '3 = Polarized Lembege-Pellat current sheet for electrons;'
      write(*, *) '4 = Halo electrons in the solar wind;'
      write(*, *) '5 = Electrons in the force-free current sheet;'
      write(*, *) "6 = Electrons in the injection region in the Earth's magnetotail;"
      read(*, *) icase(1)
      write(*, *) 'Distribution ID ', icase(1), ' selected!'
    endif
    call PPBICAST(icase, 1)
    if (icase(1) > 6) then
      if (kstrt == 1) then
        write(*, *) 'The desired distribution is not supported. Stop ...'
      endif
      stop
    endif

    ! get desired distribution data fxy and normalization factor fgk
    select case(icase(1))
      case(0)
        call butterfly(fxy, ncx1, ncy1, fgk, np)
      case(1)
        call lembege_pellat(fxy, ncx1, ncy1, fgk, np)
      case(2)
        call charged_LP_ics(fxy, ncx1, ncy1, fgk, np)
      case(3)
        call charged_LP_ecs(fxy, ncx1, ncy1, fgk, np)
      case(4)
        call halo(fxy, ncx1, ncy1, fgk, np)
      case(5)
        call force_free(fxy, ncx1, ncy1, fgk, np)
      case(6)
        call power_maxwellian(fxy, ncx1, ncy1, fgk, np)
    end select

    ! calculate partition variables: edges, nyp, noff, nypmx
    fmargy = (fxy(1, :) + fxy(ncx1, :)) * 0.5 + sum(fxy(2:ncx, :), 1)
    call pfedges2(edges, nyp, noff, fmargy, nypmx, nypmn, ncy, kstrt, nvp, idps)
    ! call pdcomp2(edges, nyp, noff, nypmx, nypmn, ncy, kstrt, nvp, idps)

    if (kstrt == 1) then
      write(*, *) 'nvp:'
      write(*, *) nvp
      write(*, *) 'nypmx:'
      write(*, *) nypmx
      write(*, *) 'nypmn:'
      write(*, *) nypmn
    endif

    ! check if too many partitions
    if (nypmn == 0) then
      if (kstrt == 1) write(*, *) 'too many partitions that leads to nypmn = 0'
      call PPEXIT()
      stop
    endif

    ! allocate data
    allocate(density(ncx1, nypmx), scr(ncx1), scrg(ncx1, nypmx))

    ! inverse transform sampling
    npp = 0
    call pp_inv_sampling_2D(part, npp, fxy, npx, npy, ncx1, ncy1, edges,&
      &idimp, npmax, idps, ierr, tol, eps, tcheb, tbisc)

    if (kstrt == nvp) then
      write(*, *) 'ierr:'
      write(*, *) ierr
      write(*, *) 'npp:'
      write(*, *) npp
      write(*, *) 'tcheb [seconds]:'
      write(*, *) tcheb
      write(*, *) 'tbisc [seconds]:'
      write(*, *) tbisc
    endif

    ! deposit particle density
    density = 0.0
    if (npp > 0) then
      call ppgpost2l(part, density, npp, noff, fgk, idimp, npmax, ncx1, nypmx)
    endif
    ! add guard cells
    call PPNLAGUARD2L(density, scr, nyp, ncx, kstrt, nvp, ncx1, nypmx)
    ! copy guard cells
    call PPNLCGUARD2L(density, nyp, kstrt, nvp, ncx1, nypmx) 

    ! open file
    if (kstrt == 1) then
      open(unit=1, file='data.bin', access='stream',&
        form='unformatted', status='new')
    endif
    ! collect and write data
    call PPWRVNDATA2(density, scrg, 1, ncx1, nyp, nypmx, 1)
    if (kstrt == 1) then
      ! write density function and close file
      write(unit=1) ((fxy(j, k), j = 1, ncx1), k = 1, ncy)
      close(unit=1)
    endif

    call PPEXIT()
  end subroutine
!
end program
