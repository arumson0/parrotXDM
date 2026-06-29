subroutine c9_loops(tau,mtrx,c6,rc,c9,zinv,a1,a2,zdamp,damp,rmax2,n,l,e9_shell)
use omp_lib
implicit none
integer, intent(in) :: n, l
real(8), intent(in) :: tau(l,3), mtrx(3,3)
real(8), intent(in) :: c6(l,l), zinv(l,l), rc(l,l)
real(8), intent(in) :: c9(l,l,l)
real(8), intent(in) :: a1, a2, zdamp
integer, intent(in) :: damp
real(8), intent(in) :: rmax2
real(8), intent(out) :: e9_shell

integer :: nl1,nl2,nl3,ml1,ml2,ml3,i,j,k
real(8) :: rij2, rik2, rjk2, rij, rik, rjk
real(8) :: cosi, cosj, cosk, g, f
real(8) :: ri(3), rj(3), rk(3)
real(8) :: rijv(3), rikv(3), rjkv(3)
real(8) :: Tn(3), Tm(3)
real(8) :: rij3, rik3, rjk3
real(8) :: fij, fik, fjk

e9_shell = 0.0d0

! Parallelize the 6 outer loops with reduction
! Collapse tells OpenMP to treat nested loops as a single loop for better load balancing
! Reduction handles the e9_shell sum across threads
! Private lists all temporaries used inside loops
! No need for e_local or atomic
!$omp parallel do collapse(6) reduction(+:e9_shell) default(shared) &
!$omp private(nl1,nl2,nl3,ml1,ml2,ml3,i,j,k) &
!$omp private(ri,rj,rk,rijv,rikv,rjkv,Tn,Tm) &
!$omp private(rij2,rik2,rjk2,rij,rik,rjk) &
!$omp private(cosi,cosj,cosk,g,f) &
!$omp private(rij3,rik3,rjk3) &
!$omp private(fij,fik,fjk) &
!$omp schedule(dynamic, 100)

do nl1 = -n, n
  do nl2 = -n, n
    do nl3 = -n, n
      do ml1 = -n, n
        do ml2 = -n, n
          do ml3 = -n, n

            if (.not.(abs(nl1)==n .or. abs(nl2)==n .or. abs(nl3)==n)) cycle

            Tn = nl1*mtrx(:,1) + nl2*mtrx(:,2) + nl3*mtrx(:,3)
            Tm = ml1*mtrx(:,1) + ml2*mtrx(:,2) + ml3*mtrx(:,3)

            do i = 1, l
              ri = tau(i,:)
              do j = 1, l
                rj = tau(j,:) + Tn
                rijv = ri - rj
                rij2 = dot_product(rijv, rijv)
                if (rij2 == 0.0d0) cycle
                if (rij2 > rmax2 ) cycle
                rij = sqrt(rij2)

                do k = 1, l
                  rk = tau(k,:) + Tm
                  rikv = ri - rk
                  rjkv = rj - rk

                  rik2 = dot_product(rikv, rikv)
                  rjk2 = dot_product(rjkv, rjkv)
                  if (rik2 == 0.0d0 .or. rjk2 == 0.0d0) cycle
                  if (rik2 > rmax2  .or. rjk2 > rmax2 ) cycle

                  rik = sqrt(rik2)
                  rjk = sqrt(rjk2)

                  cosi = dot_product(rijv, rikv) / (rij*rik)
                  cosj =-dot_product(rijv, rjkv) / (rij*rjk)
                  cosk = dot_product(rikv, rjkv) / (rik*rjk)
                  g = 3.d0*cosi*cosj*cosk + 1.d0

                  rij3 = rij2 * rij
                  rik3 = rik2 * rik
                  rjk3 = rjk2 * rjk

                  if (damp == 0) then
                    f = (rij3 / (rij3 + (a1*rc(i,j)+a2)**3.d0)) &
                        * (rik3 / (rik3 + (a1*rc(i,k)+a2)**3.d0)) &
                        * (rjk3 / (rjk3 + (a1*rc(j,k)+a2)**3.d0))
                  else if (damp == 1) then
                    fij = (rij3) / (rij3 + (zdamp*c6(i,j)*zinv(i,j))**(1.d0/2.d0))
                    fik = (rik3) / (rik3 + (zdamp*c6(i,k)*zinv(i,k))**(1.d0/2.d0))
                    fjk = (rjk3) / (rjk3 + (zdamp*c6(j,k)*zinv(j,k))**(1.d0/2.d0))
                    f = fij*fik*fjk
                  else if (damp == 2) then
                    fij = (rij3*rij3 / (rij3*rij3 + (a1*rc(i,j) + a2)**6.d0 ))**(1.d0/2.d0)
                    fik = (rik3*rik3 / (rik3*rik3 + (a1*rc(i,j) + a2)**6.d0 ))**(1.d0/2.d0)
                    fjk = (rjk3*rjk3 / (rjk3*rjk3 + (a1*rc(i,j) + a2)**6.d0 ))**(1.d0/2.d0)
                    f = fij*fik*fjk
                  endif


                  e9_shell = e9_shell + c9(i,j,k)*g*f/(rij3*rik3*rjk3)

                end do
              end do
            end do

          end do
        end do
      end do
    end do
  end do
end do

!$omp end parallel do

end subroutine
