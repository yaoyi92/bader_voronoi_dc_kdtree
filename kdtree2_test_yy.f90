
!
! MAIN PROGRAM HERE.  This is just an example so you know
! how to use it.   This file is in the public domain.
! 

module time_kdtree
  use kdtree2_precision_module
  use kdtree2_module
contains

  real function time_search(tree,nsearch,mode, nn, r2)
    !
    !  Return CPU time, in seconds, for searching 'nsearch' reference points
    !  using any specific search mode. 
    !

    type(kdtree2), pointer :: tree
    integer, intent(in)               :: nsearch ! how many reference points
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real(kdkind), intent(in)                  :: r2      ! radius^2
    !
    real(kdkind) :: qv(tree%dimen), rv               ! query vector, random variate
    integer :: i, random_loc, nf
    real    :: t0, t1
    type(kdtree2_result), allocatable :: results(:) 

    real(kind(0.0d0)) :: nftotal
    call cpu_time(t0)   ! get initial time.
    nftotal = 0.0

    allocate(results(nn)) 
    do i=1,nsearch

       select case (mode)
       case (1)
          !
          !  Fixed NN search around randomly chosen point
          !
          call random_number(qv)
          call kdtree2_n_nearest(tp=tree,qv=qv,nn=nn,results=results) 
       case (2)
          !
          ! Fixed NN seasrch around randomly chosen point in dataset
          ! with 100 correlation time
          !
          call random_number(rv)
          random_loc = floor(rv*tree%n) + 1
          call kdtree2_n_nearest_around_point(tp=tree,idxin=random_loc,&
           correltime=100,nn=nn,results=results) 
       case (3)
          !
          ! fixed r2 search
          !
          call random_number(qv)
          call kdtree2_r_nearest(tp=tree,qv=qv,r2=r2,nfound=nf,&
           nalloc=nn,results=results)
          nftotal = nftotal+nf
       case default
          write (*,*) 'Search type ', mode, ' not implemented.'
          time_search = -1.0 ! invalid
          return
       end select

    enddo
    call cpu_time(t1)

    time_search = t1-t0
    if (nftotal .gt. 0.0) then
!       write (*,*) 'Average number of neighbors found = ', nftotal / real(nsearch)
    endif
    deallocate(results)
    return
  end function time_search

  real function searches_per_second(tree,mode,nn,r2) result(res)
    !
    !
    ! return estimated number of searches per second.
    ! Will call "time_search" with increasing numbers of reference points 
    ! until CPU time taken is at least 1 second.
    !
    type(kdtree2), pointer :: tree
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real(kdkind), intent(in)                  :: r2      ! radius^2
    !
    integer :: nsearch
    real    :: time_taken

    nsearch = 50  ! start with 50 reference points
    do
       time_taken = time_search(tree,nsearch,mode,nn,r2)
       if (time_taken .lt. 1.0) then
!          write (*,*) 'Intermediate result : ', time_taken, &
!           'sec, nref=',nsearch
          nsearch = nsearch * 5
          cycle
       else
          res = real(nsearch) / time_taken
          return
       end if
    end do
    return
  end function searches_per_second

  real function average_number_within_ball(tree,navg,r2) result(res)
    ! return the arithmetical average number of points within ball of size
    ! 'r2'
    !
    !  log(avg) = 1/navg log(N_within_ball(i))
    type(kdtree2), pointer :: tree
    integer,intent(in) :: navg
    real(kdkind), intent(in)   :: r2
    !
    integer :: i, cnt
    real(kdkind) :: sum, qv(tree%dimen)
    sum = 0.0
    do i=1,navg
       call random_number(qv)
       cnt = kdtree2_r_count(tree,qv,r2)
       if (cnt .gt. 0)  sum = sum + real(cnt)
    end do

    res = sum/real(navg)
    return
  end function average_number_within_ball

end module time_kdtree

program kd_tree_test
  use kdtree2_module
  use time_kdtree

  integer :: n, d 
  real(kdkind) :: my_array(3,10)
  !real(kdkind), allocatable :: query_vec(:)

  type(kdtree2), pointer :: tree, tree2, tree3
  ! this is how you declare a tree in your main program

  integer :: k

  type(kdtree2_result),allocatable :: results(:), resultsb(:)
  integer   :: nnbrute, rind
  real      :: t0, t1, sps, avgnum, maxdeviation
  real(kdkind) :: rv 
  integer, parameter  :: nnn = 5
  integer, parameter  :: nr2 = 5
  real r2array(5)
  data r2array / 1.0e-4,1.0e-3,5.0e-3,1.0e-2,2.0e-2 / 
  integer   :: nnarray(5)
  data nnarray / 1, 5, 10, 25, 500/ 
  real(kdkind) :: qv(3)               ! query vector, random variate


  !allocate(my_array(d,n))
  !allocate(query_vec(d))
  n = 1
  allocate(results(n))
  call random_number(my_array)  !fills entire array with built-in-randoms

  call cpu_time(t0)
  tree => kdtree2_create(my_array,sort=.false.,rearrange=.false.)  ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for non-rearranged tree.'

  call random_number(qv)

  call kdtree2_n_nearest(tp=tree,qv=qv,nn=n,results=results) 
  write (*,*) my_array
  write (*,*) qv
  write (*,*) results


  


  call kdtree2_destroy(tree)  

  ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
  ! TO MAKE THE TREE.  

  !deallocate(my_array)
  ! deallocate the memory for the data.

end program kd_tree_test
