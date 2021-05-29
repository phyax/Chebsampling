!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications
! pplib2.f90 contains basic communications procedures for 1d partitions
! PPINIT2 initializes parallel processing for Fortran90, returns
!         number of processors and processor id.
! PPEXIT terminates parallel processing.
! PPABORT aborts parallel processing.
! PWTIMERA performs parallel local wall clock timing.
! PPSUM performs parallel sum of a real vector.
! PPDSUM performs parallel sum of a double precision vector.
! PPMAX performs parallel maximum of a real vector.
! PPIMAX performs parallel maximum of an integer vector.
! PPDMAX performs parallel maximum of a double precision vector.
! PPDSCAN performs parallel prefix reduction of a double precision
!         vector
! PPBICAST broadcasts integer data from node 0
! PPBDCAST broadcasts double precision data from node 0
!
! for non-periodic field boundary condtions
! PPNLCGUARD2L copies data to interior guard cells in y for scalar data,
!              linear interpolation, and distributed data with
!              non-uniform partition.
! PPNLAGUARD2L adds interior guard cells in y for scalar array, linear
!              interpolation, and distributed data with non-uniform
!              partition.
!
! PPWRVNDATA2 collects distributed real 2d vector non-uniform data f and
!             writes to a fortran unformatted file
!
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
      module pplib2
      use mpi
      implicit none
!
! common data for parallel processing
! lstat = length of status array
      integer, parameter :: lstat = MPI_STATUS_SIZE
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! msum = MPI_SUM
! mmax = MPI_MAX
      integer :: msum, mmax
      save
!
      private
      public :: lstat, nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      public :: PPINIT2, PPEXIT, PPABORT, PWTIMERA
      public :: PPSUM, PPDSUM, PPMAX, PPIMAX, PPDMAX, PPDSCAN
      public :: PPBICAST, PPBDCAST
      public :: PPNLCGUARD2L, PPNLAGUARD2L
      public :: PPWRVNDATA2
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
! this subroutine initializes parallel processing
! lgrp communicator = MPI_COMM_WORLD
! output: idproc, nvp
! idproc = processor id in lgrp communicator
! nvp = number of real or virtual processors obtained
      implicit none
      integer, intent(inout) :: idproc, nvp
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
! msum = MPI_SUM
! mmax = MPI_MAX
! local data
      integer :: ierror, ndprec, idprec
      integer :: iprec
      logical :: flag
      real :: prec
! ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (digits(prec) > 24) then
         ndprec = 1
      else
         ndprec = 0
      endif
! idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
      if (digits(iprec) > 31) then
         idprec = 1
      else
         idprec = 0
      endif
! this segment is used for mpi computers
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
! initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror /= 0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
! set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
! single precision real
      if (ndprec==0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
! double precision real
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
! single precision integer
!     if (idprec==0) then
!        mint = MPI_INTEGER
! double precision integer
!     else
!        mint = MPI_INTEGER8
!     endif
! operators
      msum = MPI_SUM
      mmax = MPI_MAX
      nvp = nproc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT()
! this subroutine terminates parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
! synchronize processes
         call MPI_BARRIER(lworld,ierror)
! terminate MPI execution environment
         call MPI_FINALIZE(ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT()
! this subroutine aborts parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: errorcode, ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
         errorcode = 1
! terminate MPI execution environment
         call MPI_ABORT(lworld,errorcode,ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
! this subroutine performs local wall clock timing
! input: icntrl, dtime
! icntrl = (-1,0,1) = (initialize,ignore,read) clock
! clock should be initialized before it is read!
! time = elapsed time in seconds
! dtime = current time
! written for mpi
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
! local data
      double precision :: jclock
! initialize clock
      if (icntrl==(-1)) then
         dtime = MPI_WTIME()
! read clock and write time difference from last clock initialization
      else if (icntrl==1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output real data
! g = scratch real array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! perform sum
      call MPI_ALLREDUCE(f,g,nxp,mreal,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! perform sum
       call MPI_ALLREDUCE(f,g,nxp,mdouble,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMAX(f,g,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, f(j,k) = maximum as a function of k of f(j,k)
! at the end, all processors contain the same maximum.
! f = input and output real data
! g = scratch real array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! mmax = MPI_MAX
! local data
      integer j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mreal,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, if(j,k) = maximum as a function of k of if(j,k)
! at the end, all processors contain the same maximum.
! if = input and output integer data
! ig = scratch integer array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mint = default datatype for integers
! mmax = MPI_MAX
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(if,ig,nxp,mint,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         if(j) = ig(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, f(j,k) = maximum as a function of k of f(j,k)
! at the end, all processors contain the same maximum.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! mmax = MPI_MAX
! local data
      integer j, ierr
! return if only one processor
      if (nproc==1) return
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSCAN(f,g,nxp)
! this subroutine performs a parallel prefix reduction of a vector,
! that is: f(j,k) = sum over k of f(j,k), where the sum is over k values
! less than idproc.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! msum = MPI_SUM
! local data
      integer :: j, ierr
! return if only one processor
      if (nproc==1) return
! performs a parallel prefixm sum
       call MPI_SCAN(f,g,nxp,mdouble,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBICAST(if,nxp)
! this subroutine broadcasts integer data from node 0
! if = input and output integer data
! nxp = number of data values
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mint = default datatype for integers
! local data
      integer :: ierr
! return if only one processor
      if (nproc==1) return
! broadcast integer
      call MPI_BCAST(if,nxp,mint,0,lgrp,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBDCAST(f,nxp)
! this subroutine broadcasts double precision data from node 0
! f = input and output double precision data
! nxp = number of data values
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mdouble = default double precision type
! local data
      integer :: ierr
! return if only one processor
      if (nproc==1) return
! broadcast integer
      call MPI_BCAST(f,nxp,mdouble,0,lgrp,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! this subroutine copies data to guard cells in non-uniform partitions
! guard cell on last processor is presumed already set.
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f
! nyp = number of primary gridpoints in field partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cell.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! special case for one processor
      if (nvp==1) return
      ks = kstrt - 1
      moff = nypmx*nvp + 2
! copy to guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kr < nvp) then
         call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      endif
      if (kl >= 0) then
         call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      endif
      if (kr < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNLAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! for scalar data.  no copying is done at the boundary edges.
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, nx1, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) return
      ks = kstrt - 1
      moff = nypmx*nvp + 1
! add guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kl >= 0) then
         call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      endif
      if (kr < nvp) then
         call MPI_SEND(f(1,ks),nxv,mreal,kr,moff,lgrp,ierr)
         do j = 1, nx1
            f(j,ks) = 0.0
         enddo
      endif
      if (kl >= 0) then
         call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
         do j = 1, nx1
            f(j,1) = f(j,1) + scr(j)
         enddo
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVNDATA2(f,g,ndim,nxv,nyp,nypmx,iunit)
! this subroutine collects distributed periodic real 2d vector
! non-uniform data f and writes to a fortran unformatted file
! f = input data to be written
! g = scratch data
! ndim = first dimension of data array f
! nyp = actual data written for third dimension
! nxv = second dimension of data array f
! nypmx = third dimension of data array f
! iunit = fortran unit number
! input: all
      implicit none
      integer, intent(in) :: ndim, nxv, nyp, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nyp) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: nvp, idproc, nps, nnxv, nnyp, nnypx, id, i, j, k, n
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nnxv = ndim*nxv
      nnyp = nnxv*nyp
      nnypx = nnxv*nypmx
! this segment is used for shared memory computers
!     write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nyp)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
! node 0 receives messages from other nodes
      if (idproc==0) then
! first write data for node 0
         write (unit=iunit) (((f(n,j,k),n=1,ndim),j=1,nxv),k=1,nyp)
! then write data from remaining nodes
         do i = 2, nvp
         id = i - 1
         call MPI_IRECV(g,nnypx,mreal,id,99,lgrp,msid,ierr)
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         nps = nps/nnxv
         write (unit=iunit) (((g(n,j,k),n=1,ndim),j=1,nxv),k=1,nps)
         enddo
! other nodes send data to node 0
      else
         call MPI_SEND(f,nnyp,mreal,0,99,lgrp,ierr)
      endif
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
      use pplib2, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use pplib2, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use pplib2, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use pplib2, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use pplib2, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use pplib2, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMAX(f,g,nxp)
      use pplib2, only: SUB => PPMAX
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use pplib2, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use pplib2, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSCAN(f,g,nxp)
      use pplib2, only: SUB => PPDSCAN
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBICAST(if,nxp)
      use pplib2, only: SUB => PPBICAST
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if
      call SUB(if,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPBDCAST(f,nxp)
      use pplib2, only: SUB => PPBDCAST
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f
      call SUB(f,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
      use pplib2, only: SUB => PPNLCGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      call SUB(f,nyp,kstrt,nvp,nxv,nypmx)
      end subroutine
!
      subroutine PPNLAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      use pplib2, only: SUB => PPNLAGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPWRVNDATA2(f,g,ndim,nxv,nyp,nypmx,iunit)
      use pplib2, only: SUB => PPWRVNDATA2
      implicit none
      integer, intent(in) :: ndim, nxv, nyp, nypmx, iunit
      real, intent(in), dimension(ndim,nxv,nyp) :: f
      real, intent(inout), dimension(ndim,nxv,nypmx) :: g
      call SUB(f,g,ndim,nxv,nyp,nypmx,iunit)
      end subroutine
!
