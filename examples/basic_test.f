
      include 'lplib3.ins'
	  integer*8 li
	  integer*4 i,ti,n
	  parameter (n=100000000)
      real*8 t,v1(n),v2(n),v3(n)
	  real*4 acc
	  external addvec

	  li = initparallel(0)
      if(li.eq.0) STOP ' cannot init LPLIB3'
	  print*, 'initparallel = ', li

	  ti = newtype(li,n)
      if(ti.eq.0) STOP ' cannot init new type'
	  print*, 'newtype = ', ti

      do i = 1,n
		  v1(i) = i
		  v2(i) = 2*i
      end do

	  t = getwallclock()
	  acc = launchparallel(li, ti, 0, addvec, 3, v1, v2, v3)
	  print*, 'wall time = ', getwallclock() - t, 'accelerate = ', acc

	  call stopparallel(li)

      end


      SUBROUTINE addvec(bi,ei,ti,v1,v2,v3)
	  INTEGER i,bi,ei,ti
	  REAL*8 v1(*),v2(*),v3(*)
	  print*, 'thread = ', ti, 'loop = ',bi,ei

	  do i = bi,ei
	  v3(i)=exp(log(acos(cos(v1(i)+1))))+exp(log(acos(cos(v2(i)+2))))
	  end do
	  return
      end
