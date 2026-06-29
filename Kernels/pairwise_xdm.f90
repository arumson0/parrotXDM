subroutine pairwise_xdm(tau, mtrx, rc, c6, c8, c10, zinv, a1, a2, zdamp, damp, rmax2, n_vecs, l, exdm, e6, e8, e10)
    use omp_lib
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: n_vecs(3)
    real(8), intent(in) :: tau(l,3), mtrx(3,3)
    real(8), intent(in) :: rc(l,l), c6(l,l), c8(l,l), c10(l,l), zinv(l,l)
    real(8), intent(in) :: a1, a2, zdamp, rmax2
    integer, intent(in) :: damp  ! 0=bj, 1=z
    real(8), intent(out) :: exdm, e6, e8, e10

    integer :: nl1,nl2,nl3,i,j
    real(8) :: rij, rij2, f6,f8,f10
    real(8) :: coord_i(3), coord_j(3), vec_rij(3)
    real(8) :: e_local
    real(8) :: e6_local, e8_local, e10_local

    exdm = 0.0d0
    e6   = 0.0d0
    e8   = 0.0d0
    e10  = 0.0d0

    !$omp parallel default(shared) private(nl1,nl2,nl3,i,j,coord_i,coord_j,vec_rij,rij,rij2,f6,f8,f10,e_local) reduction(+:exdm)
    e_local = 0.0d0
    !$omp do collapse(3)
    do nl1 = -n_vecs(1), n_vecs(1)
      do nl2 = -n_vecs(2), n_vecs(2)
        do nl3 = -n_vecs(3), n_vecs(3)
          if (abs(nl1)==n_vecs(1) .or. abs(nl2)==n_vecs(2) .or. abs(nl3)==n_vecs(3)) then
            do i = 1, l
              coord_i = tau(i,:)
              do j = 1, l
                coord_j = tau(j,:) + nl1*mtrx(1,:) + nl2*mtrx(2,:) + nl3*mtrx(3,:)
                vec_rij = coord_i - coord_j
                rij2 = dot_product(vec_rij, vec_rij)
                if (rij2 == 0.0d0) cycle
                if (rij2 > rmax2) cycle
                rij = sqrt(rij2)

                if (damp == 0 .or. damp == 2) then
                  f6  = 1.d0 / (rij**6 + (a1*rc(i,j) + a2)**6)
                  f8  = 1.d0 / (rij**8 + (a1*rc(i,j) + a2)**8)
                  f10 = 1.d0 / (rij**10 + (a1*rc(i,j) + a2)**10)
                else if (damp == 1) then
                  f6  = 1.d0 / (rij**6  +  c6(i,j)*zdamp*zinv(i,j))
                  f8  = 1.d0 / (rij**8  +  c8(i,j)*zdamp*zinv(i,j))
                  f10 = 1.d0 / (rij**10 + c10(i,j)*zdamp*zinv(i,j))
                end if

                e6_local  = -c6(i,j)*f6
                e8_local  = -c8(i,j)*f8
                e10_local = -c10(i,j)*f10
                e_local = -(c6(i,j)*f6 + c8(i,j)*f8 + c10(i,j)*f10)

                exdm = exdm + e_local
                e6   = e6  + e6_local
                e8   = e8  + e8_local
                e10  = e10 + e10_local
              end do
            end do
          end if
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
end subroutine
