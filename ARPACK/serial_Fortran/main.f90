program arpack_diag
  implicit none

  ! Problem sizes / parameters
  integer, parameter :: N    = 1000
  integer, parameter :: nev  = 9
  integer, parameter :: ncv  = 2*nev + 1
  integer, parameter :: ldv  = N
  integer, parameter :: ldz  = N
  integer, parameter :: lworkl = ncv * (ncv + 8)

  ! ARPACK controls
  integer :: iparam(11), ipntr(11)
  integer :: info, ido
  character(len=1) :: bmat
  character(len=2) :: which
  character(len=1) :: howmny
  logical :: rvec
  logical :: select(ncv)

  ! Work / results
  real(8) :: tol, sigma
  real(8) :: resid(N)
  real(8) :: V(ldv, ncv)
  real(8) :: Z(ldz, nev)
  real(8) :: D(nev)
  real(8) :: workd(3*N)
  real(8) :: workl(lworkl)

  ! Local
  integer :: i
  real(8) :: val, ref, eps

  ! Settings
  tol    = 1.0d-6
  sigma  = 0.0d0
  bmat   = 'I'          ! standard eigenproblem
  which  = 'LA'         ! smallest algebraic
  howmny = 'A'          ! compute all requested Ritz vectors
  rvec   = .true.       ! return eigenvectors

  ! iparam setup
  iparam = 0
  iparam(1) = 1         ! ishift
  iparam(3) = 10*N      ! max iterations
  iparam(4) = 1         ! NB (block size) must be 1 for dsaupd
  iparam(7) = 1         ! mode = 1 (standard A*x)

  ! Reverse communication loop
  info = 0
  ido  = 0
  do
    call dsaupd(ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, &
                iparam, ipntr, workd, workl, lworkl, info)

    if (ido == -1 .or. ido == 1) then
      call dmatvec(workd(ipntr(1)), workd(ipntr(2)), N)
    else
      exit
    end if
  end do

  ! Check convergence: iparam(5) is number of converged Ritz values
  if (info < 0 .or. iparam(5) < nev) then
    write(*,*) 'Error in dsaupd: iparam(5)=', iparam(5), ' nev=', nev, ' info=', info
    stop 1
  end if

  ! Extract eigenvalues/eigenvectors
  call dseupd(rvec, howmny, select, D, Z, ldz, sigma, bmat, N, which, nev, tol, &
              resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info)

  if (info < 0) then
    write(*,*) 'Error in dseupd: info=', info
    stop 1
  end if

  ! Validate results against exact eigenvalues 1..nev (ascending)
  do i = 1, nev
    val = D(i)
    ref = dble(N-nev+i)
    eps = abs(val - ref)
    write(*,'(f12.6," - ",f12.6," = ",f12.6)') val, ref, eps
    if (eps > 1.0d-5) then
      write(*,*) 'Eigenvalue ', i, ' does not match: ', val, ' vs ', ref
      stop 1
    end if
  end do

  write(*,*) 'Done'

contains

  subroutine dmatvec(x, y, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in)  :: x(n)
    real(8), intent(out) :: y(n)
    integer :: k
    do k = 1, n
      y(k) = dble(k) * x(k)   ! A is diagonal with entries 1,2,...,n
    end do
  end subroutine dmatvec

end program
 
