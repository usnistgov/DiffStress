c----------------------------------------------------------------------c
c     Rotate the coordinate of matrix a by rotation matrix 'rot' and
c     returns the rotated mat as 'b'
c
c     Transform a matrix by rotation mat rot and returns the result as b
c     b = rot  a  rot^T
c     b_ij = rot_ik a_kl rot_jl

c     rot{b_axes<-a_axes}
      subroutine matrot(a,rot,b)
      implicit none
      real*8 a(3,3), rot(3,3), b(3,3)
      integer i,j,k,l
cf2py intent(in) a, rot
cf2py intent(out) b
      do 100 i=1,3
      do 100 j=1,3
         b(i,j) = 0.d0
      do 100 k=1,3
      do 100 l=1,3
         b(i,j) = b(i,j) + rot(i,k) * a(k,l) * rot(j,l)
 100  continue
      return
      end subroutine matrot
c----------------------------------------------------------------------c
c     Rotate the coordinates of matrix a(6,6) by rotation matrix 'rot'
c     and returns the rotated mat as b(6,6) when expressed following
c     Voigt's notation.
c
c     Note that the matrices a and b are double-tensors in that each of
c     the tensors are as below:
c          a_(coordinate1,coordinate1) # ca
c          b_(coordinate2,coordinate2) # sa
c          rot_(coordinate2<-coordinate1) # sa<-ca
      subroutine matrot66(a,rot,b)
      implicit none
      real*8 a(6,6),rot(3,3),b(6,6),b3333(3,3,3,3),aux33(3,3),aux6(6),
     $     aux66(6,6),aux3333(3,3,3,3),dummy
      integer i,j,k,l,i1,j1,k1,l1
      call voigt(aux6, aux33, a, aux3333, 3)
      b(:,:) = 0.d0
      b3333(:,:,:,:) = 0.d0
      do 100 i=1,3
      do 100 j=1,3
      do 100 k=1,3
      do 100 l=1,3
         dummy = 0.d0
         do 20 i1=1,3
         do 20 j1=1,3
         do 20 k1=1,3
         do 20 l1=1,3
c     A[i,j,k,l] = r[i,i1] r[j,j1] A[i1,j1,k1,l1] r[k1,k]^t r[l1,l]^t
c     A[i,j,k,l] = r[i,i1] r[j,j1] A[i1,j1,k1,l1] r[k,k1]   r[l,l1]
            dummy = dummy + rot(i,i1) * rot(j,j1)
     $           * rot(k,k1) * rot(l,l1) * aux3333(i1,j1,k1,l1)
 20      continue
         b3333(i,j,k,l) = dummy
 100  continue
      call voigt(aux6,aux33,b,b3333,4)
      return
      end subroutine matrot66
c----------------------------------------------------------------------c
c     matrix multiplication: c_ij = a_ik * b_kj
      subroutine matply(a,b,c)
      implicit none
      real*8 a(3,3), b(3,3), c(3,3)
      integer i,j,k
cf2py intent(in) a, b
cf2py intent(out) c
c     c = a * b
      do 200 i=1,3
      do 200 j=1,3
         c(i,j) = 0.d0
      do 200 k=1,3
         c(i,j) = c(i,j) + a(i,k) * b(k,j)
 200  continue
      return
      end subroutine matply
c----------------------------------------------------------------------c
      subroutine vecply(a,b,c)
c     c[i] = b[i,j] * a[j]
      implicit none
      real*8 a(3), b(3,3), c(3)
      integer i,j
      do 200 i=1,3
         c(i) = 0.d0
      do 200 j=1,3
         c(i) = c(i) + b(i,j) * a(j)
 200  continue
      return
      end subroutine vecply
c----------------------------------------------------------------------c
c     indentity matrix returns.
      subroutine imat(a)
      implicit none
      integer i,j
      real*8 a(3,3)
cf2py intent(out) a
      do 100 i=1,3
      do 100 j=1,3
         a(i,j) = 0.d0
         if (i.eq.j) a(i,j) = 1.d0
 100  continue
      return
      end subroutine imat
c----------------------------------------------------------------------c
c     zeros matrix returns.
      subroutine zmat(a)
      implicit none
      integer i,j
      real*8 a(3,3)
cf2py intent(out) a
      do 100 i=1,3
      do 100 j=1,3
         a(i,j) = 0.d0
 100  continue
      return
      end subroutine zmat
c----------------------------------------------------------------------c
      subroutine zvec(a)
      implicit none
      integer i
      real*8 a(3)
      do i=1,3
         a(i) = 0.d0
      enddo
      return
      end subroutine zvec
c----------------------------------------------------------------------c
c     Cross operator can be used to obtain the cross product as a
c     matrix-vector product.
c     a x b = [a_x]_ij [b_j]
      subroutine crossop(u, m)
      implicit none
      real*8 u(3), m(3,3)
      integer i,j
cf2py intent(in) u
cf2py intent(out) m
c     Calculation of cross product operator
c     a x b = [a]_x b
c     Here [a]_x, is calculated and returned for vector u.
      do 100 i=1,3
      do 100 j=1,3
         m(i,j) = 0.d0
 100  continue
      m(1,2) = -u(3)
      m(1,3) =  u(2)
      m(2,1) =  u(3)
      m(2,3) = -u(1)
      m(3,1) = -u(2)
      m(3,2) =  u(1)
      end subroutine crossop
c----------------------------------------------------------------------c
c     Copy vector a to b.
      subroutine copyvect(a,b)
      implicit none
      integer i
      real*8 a(3), b(3)
      do i=1,3
         b(i) = a(i)
      enddo
      return
      end subroutine
c----------------------------------------------------------------------c
c     Copy matrix a(m,n) to b(m,n)
      subroutine copymatmn(a,b,m,n)
      implicit none
      integer i,j,m,n
      real*8 a(m,n), b(m,n)
      do 30 i=1,m
      do 30 j=1,n
         b(i,j) = a(i,j)
 30   continue
      return
      end subroutine
c----------------------------------------------------------------------c
c     Copy matrix a to b.
      subroutine copymat(a,b)
      implicit none
      integer i,j
      real*8 a(3,3), b(3,3)
      do 30 i=1,3
      do 30 j=1,3
         b(i,j) = a(i,j)
 30   continue
      return
      end subroutine
c----------------------------------------------------------------------c
c     Copy matrix a to b.
      subroutine copymatm(a,b,m)
      implicit none
      integer i,j,k,m
      real*8 a(3,3,m), b(3,3,m)
      do 30 k=1,m
      do 30 i=1,3
      do 30 j=1,3
         b(i,j,k) = a(i,j,k)
 30   continue
      return
      end subroutine
c----------------------------------------------------------------------c
c     Put a vector, a, into b's m-th 3 vector.
      subroutine addvect(a,n0,m,b)
      implicit none
      integer m, i, n0
      real*8 a(3), b(3,n0)
      do i=1,3
         b(i, m) = a(i)
      enddo

      return
      end subroutine
c----------------------------------------------------------------------c
c     Take a vector, a, from b's m-th 3 vector.
      subroutine takevect(a,n0,m,b)
      implicit none
      integer m, i, n0
      real*8 a(3), b(3,n0)
      do i=1,3
         a(i) = b(i, m)
      enddo

      return
      end subroutine
c----------------------------------------------------------------------c
c     put 33 matrix, a, into c's m-th 3x3 matrix.
      subroutine add33mat(a,n0,m,c)
      implicit none
      real*8 a(3,3), c(3,3,n0)
      integer i, j, m, n0
cf2py intent(in) a, n0, m, c
cf2py intent(out) c
      do 30 i=1,3
      do 30 j=1,3
         c(i,j,m) = a(i,j)
 30   continue
      return
      end subroutine add33mat
c----------------------------------------------------------------------c
c     Take out a 33-matrix, a, from c's m-th 3x3 matrix.
      subroutine take33mat(a,n0,m,c)
      implicit none
      real*8 a(3,3), c(3,3,n0)
      integer i,j,m,n0
cf2py intent(in) n0, m, c
cf2py intent(out) a
      do 30 i=1,3
      do 30 j=1,3
         a(i,j) = c(i,j,m)
 30   continue
      end subroutine take33mat
c----------------------------------------------------------------------c
      subroutine write33mat(m,a)
      implicit none
      real*8 a(3,3)
      integer i, j, m
      if (m.eq.-1) then
         do i=1,3
            write(*, '(3f7.3)') (a(i,j), j=1,3)
         enddo
         write(*, *)
      else
         do i=1,3
            write(m, '(3f7.3)') (a(i,j), j=1,3)
         enddo
         write(m, *)
      endif

      end subroutine write33mat
c----------------------------------------------------------------------c
      subroutine cross(a, b, c)
      implicit none
      real*8 a(3), b(3), c(3)
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine cross
c----------------------------------------------------------------------c
      real*8 function dot_prd(a,b,n)
      real*8 a(n), b(n)
      dot_prd = 0.
      do i=1, n
         dot_prd = dot_prd + a(i) * b(i)
      enddo
      end function dot_prd
c----------------------------------------------------------------------c
      real*8 function anorm(v,n)
      implicit none
      integer n, i
      real*8 v(n)
      anorm = 0.
      do i=1, n
         anorm = anorm + v(i)**2
      enddo
      anorm = sqrt(anorm)
      return
      end function anorm
c----------------------------------------------------------------------c
      subroutine mattrs(n,a,at)
      implicit none
      integer n, i, j
      real*8 a(n,n), at(n,n)
      do 100 i=1,n
      do 100 j=1,n
         at(i,j) = a(j,i)
 100  continue
      end subroutine mattrs
c----------------------------------------------------------------------c
      logical function isdecentrot(a, n, nt)
      implicit none
      integer i, j, k, n, nt
      logical isrotmat
      real*8 a(nt, n, n), dum(n, n)
      isdecentrot = .true.
      do 62 i=1, nt             ! (nt, n, n)
         do 61 j=1, n
         do 61 k=1, n
            dum(j, k) = a(i, j, k)
 61      continue

         if (.not.isrotmat(dum, n)) then
            isdecentrot = .false.
            goto 63
         endif
 62   enddo
 63   end function isdecentrot ! end of function isdecentrot

c     check if the given matrix is a rotation matrix!                  c
c     1). Orthogonality: A^T = A^-1                                    c
c     2). Det(A) = 1.                                                  c
      logical function isrotmat(a, n)
      implicit none
      integer i, j, k, l, m, n, irst
      real*8 a(n, n), inva(n, n), tiny, detmnt, dum(n,n),
     $     idx(3,3), det, dum02(n,n)
      logical ismateq
cf2py intent(in) a, n
cf2py intent(out) a
      isrotmat = .true.
      tiny = 0.000000001
      detmnt = det(a)
      if (abs(1-detmnt).ge.tiny) then
         write(*,*) 'Determinant of the matrix is not unity: ', detmnt
         isrotmat = .false.
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c     check if the a_ik inva_kj = identity
      do 50 i=1,n
      do 50 j=1,n
         inva(i,j) = a(i,j)
 50   continue
      call lu_inverse(inva, 3)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
      do 69 i=1,n
      do 69 j=1,n
         if (abs(a(j,i)-inva(i,j)).ge.tiny) then
            write(*,*) 'Inversed matrix is not equivalent',
     $           ' to transposed one.'
            write(*,*) 'R and R^T'
            do k=1,n
               write(*,'(3f8.3,2x, 3f8.3)')
     $              (a(k,m), m=1,n),
     $              (a(l,k), l=1,n)
            enddo
            write(*,*) 'Inverse matrix:'
            do k=1,n
               write(*,'(3f8.3)') (inva(l,k), l=1,n)
            enddo
            isrotmat = .false.
            goto 70
         endif
 69   continue
 70   end function isrotmat ! end of function isrotmat
c----------------------------------------------------------------------c
c     Normalize n-dimensional vector u and returns
      subroutine vnorm(u,n)
      implicit none
      real*8 u(n), norm
      integer i, n
cf2py intent(in) u, n
cf2py intetn(out) u
      norm = 0.d0
      do i=1,3
         norm = norm + u(i)**2
      enddo
      norm = sqrt(norm)
      do i=1,3
         u(i) = u(i) / norm
      enddo
      return
      end subroutine vnorm
c----------------------------------------------------------------------c
c     Reads a set of vectors (3, n0) and returns only the unique set of
c     them, b(3,n0), as well as the number of unique matrices, n1.
c     So the resulting b is meant to be b(3,1:n1).
      subroutine finduniquevect(a,b,n0,n1,icen,ipr)
      implicit none
      logical isin, isdupflag(n0)
      integer i,n0,n1,ip,in, inn, ind
      real*8 a(3,n0), b(3,n0), pvect(3), aa(3,n0), nan
      logical icen, ipr
cf2py intent(in) a, n0
cf2py intent(out) b, n1
      if (n0.le.0) then
         write(*,*) 'Wrong n0 is given!'
         stop
      endif
      nan = -17.
c      nan = nan/nan
      do i=1,n0
         isdupflag(i) = .false.
      enddo
c     Set a reference vector, ip-th, among a: probing vector.
c     Scan through for the rest of the vector (ip+1 to n0)
      ip = 1
      do i=1,3
         aa(i,ip) = a(i,ip)
      enddo
c----------------------------------------------------------------------c
c     Main do-while loop.
      do while(ip.ne.n0)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c        Probing vector (being a reference vector).
c        Find if pvect is in the given set in the range between
c        ip+1 and n0
         do 10 i=1,3
            pvect(i) = a(i,ip)
 10      continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
         inn = 0 !
         do 100 in=1, n0 ! loop over entire given set.
            inn = inn + 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
            if (in.gt.ip) then  ! if in is greater than the probing front
               if (isdupflag(in)) then ! if 'in' is in the flag
                  do 20 i=1,3
                     aa(i,inn) = nan ! Put intentionally nan.
 20               continue
               else ! if in is not duplicated with the current front.
                  do 30 i=1,3
                     aa(i,inn) = a(i,in)
 30               continue
               endif
            else if (in.eq.ip) then ! if in is at the probing front
               do 31 i=1,3
                  aa(i,in) = nan
 31            continue
            else if (in.lt.ip) then ! if in is already passed by the front
               do 40 i=1,3
                  aa(i,inn) = nan ! Put intentionally nan.
 40            continue
            endif
 100     continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c         write(*,'(3f7.3)') (aa(i,inn), i=1,3)
         call isvectin(aa, pvect, n0, ind, icen, isin) ! is pmat in aa?
         if (isin) then
            isdupflag(ind) = .true.
         else
            if (ipr) then
               write(*,'(3f7.3)') (aa(i,inn), i=1,3)
               write(*,*)
               write(*,*) '---------------------'
            endif
         endif
         ip = ip + 1
      enddo ! end of the main do while loop!
      inn = 0
      do 210 in=1, n0
         if (.not.(isdupflag(in))) then
            inn = inn + 1
            do i=1,3
               b(i,inn) = a(i,in)
            enddo
            if (ipr) write(*,'(3f7.3)') (b(i,inn), i=1,3)
         endif
 210  continue
      n1 = inn
      if (ipr) then
         write(*,*) 'n1: ', n1
c         read(*,*)
      endif
      return
      end subroutine finduniquevect
c----------------------------------------------------------------------c
c     Reads a set of 3x3 matrices a(3,3,n0) and returns only the unique
c     set of them, b(3,3,n0), as well as the number of unique matrices,
c     n1. So the resulting b is meant to be b(3,3,1:n1).
c     In addition, provides the array which arrays are unique.

c     iopt=0 Find nominally the unique matrices
c     iopt=1 Find unique 'rotation' matrices.
      subroutine finduniquemat(a,b,n0,n1,tag,ipr,iopt)
      implicit none
      logical isin, isdupflag(n0), tag(n0)
      integer i, j, n0, n1, ip, in, inn, ind, iopt
      real*8 a(3,3,n0), b(3,3,n0), pmat(3,3), aa(3,3,n0), nan
      logical ipr
cf2py intent(in) a, n0, ipr, iopt
cf2py intent(out) b, n1, tag
      nan = -17.
      n1 = 0
      do i=1,n0
         isdupflag(i) = .false.
      enddo
c      i = 1
c      write(*,*) (isdupflag(i), i=1,n0)
c     Set a reference matrix, ip-th, among a: probing maxt
c     scan through for the rest of the matrix (ip+1 to n0)
c     If once a duplicate is found, flag the index to the iskipind.
c     If not, increase the ip
c     if ip.eq.n0 then finishes.
      ip = 1
      do 1 i=1,3
      do 1 j=1,3
         aa(i,j,1) = a(i,j,1)
 1    continue
      do while(ip.ne.n0) ! Loop until ip meets n0
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
c        Probing matrix (being a reference matrix).
c        Find if pmat is in the given set in the range between
c        ip+1 and n0
         do 10 i=1,3
         do 9 j=1,3
            pmat(i,j) = a(i, j, ip)
 9       continue
         if (ipr) write(*, '(3f7.3)') (pmat(i,j), j=1,3)
 10      continue
c         if (ipr) read(*,*)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
         inn = 0 !
         if (ipr) write(*,*) 'in: loop over entire given set'
         do 100 in=1, n0 ! loop over entire given set.
            if (ipr) write(*,*) 'in: ', in
            inn = inn + 1
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
            if (in.gt.ip) then  ! if in is greater than the probing front
               if (isdupflag(in)) then ! if 'in' is in the flag
                  do 20 i=1,3
                  do 20 j=1,3
                     aa(i,j,inn) = nan ! Put intentionally nan.
 20               continue
               else ! if in is not duplicated with the current front.
                  do 30 i=1,3
                  do 30 j=1,3
                     aa(i,j,inn) = a(i,j,in)
 30               continue
               endif
            else if (in.eq.ip) then ! if in is at the probing front
               do 31 i=1,3
               do 31 j=1,3
                  aa(i,j,in) = nan
 31            continue
            else if (in.lt.ip) then ! if in is already passed by the front
               do 40 i=1,3
               do 40 j=1,3
                  aa(i,j,inn) = nan ! Put intentionally nan.
 40            continue
            endif
 100     continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - c
         call ismatin(aa, pmat, n0, ind, isin, iopt) ! is pmat in aa?
         if (isin) then
            isdupflag(ind) = .true.
         else
            if (ipr) then
               do i=1, 3
                  write(*,'(3f7.3)') (aa(i,j,inn),j=1,3)
                  write(*,*) '--------------------'
               enddo
            endif
         endif
         ip = ip + 1
      enddo ! end of the main do while loop!
      inn = 0
      do 210 in=1, n0
         if (.not.(isdupflag(in))) then
            tag(in) =.true.
            inn = inn + 1
            do 110 i=1,3
            do 109 j=1,3
               b(i,j,inn) = a(i,j,in)
 109        continue
c            write(*,'(3f7.3)') (b(i,j,inn), j=1,3)
 110        continue
         else
            tag(in) = .false.
         endif
 210  continue
      n1=inn
      if(ipr) write(*,*) 'n1:', n1
      return
      end subroutine finduniquemat
c----------------------------------------------------------------------c
c     if a vector b is in set of vectors a.
      subroutine isvectin(a,b,n,ind,icen,isin)
      implicit none
      logical isin, isvecteq
      real*8 a(3,n), b(3), dum(3), bm(3)
      integer i,n,in,ind
      logical icen
cf2py intent(in) a,b,n,icen
cf2py intent(out) ind,isin
      isin = .false.
      ind = -1
      do i=1,3
         bm(i) = -b(i)
      enddo
      do 100 in=1,n
         do 10 i=1,3
            dum(i) = a(i,in)
 10      continue
         if (isvecteq(dum,b)) then
            isin = .true.
            ind = in
            goto 110
         endif
         if (.not.icen) then
            if (isvecteq(dum,bm)) then
               isin = .true.
               ind = in
               goto 110
            endif
         endif
 100  continue
 110  return
      end subroutine
c----------------------------------------------------------------------c
c     if mat b is in a, returns ind and isin=.true.,
c     which otherwise is set to be .false.

c     iopt=0 Find nominally the unique matrices
c     iopt=1 Find unique 'rotation' matrices
      subroutine ismatin(a,b,n,ind,isin,iopt)
      implicit none
      logical isin, ismateq, iseqrot
      real*8 a(3,3,n), b(3,3), dum(3,3), nan
      integer i, j, n, in, ind, iopt

cf2py intent(in) a, b, n, iopt
cf2py intent(out) ind, isin
      isin = .false.
      nan = -17.
      ind = -1
      do 100 in=1,n
         do 10 i=1,3
         do 10 j=1,3
            dum(i,j) = a(i,j,in)
 10      continue
         if (iopt.eq.0) then
            isin = ismateq(dum, b, 3)
         else if (iopt.eq.1) then
            isin = iseqrot(dum,b)
         else
            write(*,*) 'Wrong option is given'
            stop
         endif
         if (isin) then
            ind = in
            goto 110
         endif
c$$$         if (iopt.eq.0) then
c$$$            if (ismateq(dum,b,3)) then
c$$$               isin = .true.
c$$$               ind = in
c$$$               goto 110
c$$$            endif
c$$$         else if (iopt.eq.1) then
c$$$            if (iseqrot(dum,b)) then
c$$$               isin = .true.
c$$$               ind = in
c$$$               goto 110
c$$$         endif
 100  continue
 110  return
      end subroutine
c----------------------------------------------------------------------c
c     If matrices a and b are equivalent?
      logical function ismateq(a,b,n)
      integer i,j,n
      real*8 a(n,n), b(n,n), tiny
      tiny = 1.d-6
      ismateq = .true.
      do 10 i=1,n
      do 10 j=1,n
         if (abs(a(i,j) - b(i,j)).gt.tiny) then
            ismateq = .false.
            goto 20
         endif
 10   continue
 20   return
      end function ismateq
c----------------------------------------------------------------------c
c     If vectors a and b are equal?
      logical function isvecteq(a,b)
      integer i
      real*8 a(3), b(3), tiny
c      tiny = 1.d-10
c      write(*,*) 'icen', icen
      tiny = 1.d-10
      isvecteq = .true.
      do 10 i=1,3
         If (abs(a(i)-b(i)).gt.tiny) then
            isvecteq = .false.
            goto 20
         endif
 10   continue
 20   return
      end function isvecteq
c     In below are some useful libraries                               c
c----------------------------------------------------------------------c
      logical function isperpen(a, b, n)
      implicit none
      real*8 trivial, sum, a(n), b(n)
      integer n, i
      trivial = 1.d-10
      sum = 0.D0
      do i=1, n
         sum = sum + a(i) * b(i)
      enddo
      if (abs(sum).lt.trivial) then
         isperpen = .true.
      else
         isperpen = .false.
      endif
      return
      end function isperpen ! end of function isperpen


      include "LIBRARY.SUB"
