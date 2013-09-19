!****
!-----------------------------------------------------------------------
!**** Prologue to LUPACK
!**** LUPACK is a set of 5 double precision subroutines for the
!**** solution and analysis of general systems of linear equations of 
!**** the form  Ax=r, where A is an N-by-N matrix and where r and x
!**** are rank N column vectors.
!**** Author= Micah Dembo of Boston University
!**** This version dated 11/22/04 was created for student use in
!**** the BME programming and modelling course: BE703.
!-----------------------------------------------------------------------
!****
!**** The routines of lupak are all based on the general approach of LU
!**** factorization with partial pivoting and are structured along the
!**** lines of similar routines that can be found in the classical 
!**** LINPACK collection that can be downloaded from netlib.  
!**** The LUPACK routines are not, in general, simply updated versions
!**** of certain linpack routines. The important distinguishing features 
!**** are as follows:
!**** 1)LUPACK routines are all entirely self contained. They do not 
!**** use subsiderary routines like the BLAS routines used by LINPACK.
!**** This is to insure that all oporations are transparent and can be
!**** easily modified for special situations (eg the default loop limits
!**** can be adjusted to increase speed if A is banded). 
!**** 2)LUPACK routines are written entirely in standard f90 rather
!**** than f77. There are no common blocks, line numbers,go-to's or
!**** other archaic features of F77.
!**** 3)LUPACK routines make maximal use of the whole array operations,
!**** dynamic memory allocation and intrinsic functions provided by the
!**** f90 standard.
!**** 4)LUPACK is a formal f90 "module". To obtain access to the
!**** contents a calling program or driver must include the statement 
!**** "use lupack" immediately after the program statement.
!**** 5)LUPACK routines are designed to be fast and to run well on both
!**** scalar and parellel machines. Hoever certain efficiencies have
!**** been sacrificed for the sake of pedigogical clarity. For example
!**** unlike LINPACK, the routines of LUPACK each have a single
!**** clear purpose and mode of operation. We avoid the use of flags
!**** which potentially allow the same routine to have several modes of
!**** operation. In the same vein, arguments to routines  are either
!**** used for input or output but never for both. Inother words we do
!**** not attempt to minimize storage requirements by using the same
!**** variable of the call sequence for both input and  output.
!**** 6)LUPACK routines include extensive inline commentaries which
!**** explain the definition of local variables and the purpose
!**** of various proceedures as they occur.
!-----------------------------------------------------------------------
!**** Alphebetical listing of routines follows.
!-----------------------------------------------------------------------
!1)   subroutine ludt(luf,ipvt,info,n,det)
!**** ludt computes the determinant of a matrix using the LU factors,
!**** the info flag, and the pivot vector returned by lufa.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!2)    subroutine lufa (a,info,ipvt,n,luf)
!****  lufa computes the singularity index, the piviot vector and the
!****  LU Factorzation of a square matrix. This is the most fundamental
!****  routine of the package.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!3)   subroutine luiv(luf,info,ipvt,n,b)
!**** luiv computes the inverse of a matrix using the LU factors
!**** ,the pivot vector, and the info flag returned by lufa.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!4)   subroutine lup1(a,n,anorm)
!**** lup1 calculates the P-1 norm of a matrix.
!**** It works directly from the definition;
!****           anorm=max(sum(abs(A(k,:)))) , k=1,2,...N
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!5)   subroutine lusl (luf,r,info,ipvt,n,x)
!**** lusl solves the real system A*x=r  using the LU factorization
!**** of the matrix A and the associated info flag and pivot vector
!**** returned by lufa.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!**** There is a test program called lutest.f in this same directory.
!**** This includes segments that illustrate the typical use of all
!**** five subroutines discribed above.
!-----------------------------------------------------------------------
!
      module lupack!begin the lupack module
!**** declare global variables.
!**** The next line can be edited if a different precision is desired.
      integer,parameter,private::rk=8! # of bytes in real numbers.
      real(kind=rk),private,parameter::ten=10,one=1,zero=0
      real(kind=rk),private::small=tiny(one)/epsilon(one)!small number
      contains 
!*******************BEGIN LUPACK ROUTINES*****************************
!
      subroutine ludt(luf,info,ipvt,n,det)
!-----------------------------------------------------------------------
!     ludt computes the determinant of a matrix using the LU factors,
!     info flag, and pivot vector returned by lufa.
!     This version dated 11/22/04 .
!     Author= Micah Dembo of Boston University
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! declare real input
      real(kind=rk),intent(in),dimension(n,n)::luf!upper factor of A
!-----------------------------------------------------------------------
! declare integer input
      integer,intent(in),dimension(n)::ipvt!pivot vector
      integer,intent(in):: info!singularity indicator flag
      integer,intent(in):: n!array dimensions
!-----------------------------------------------------------------------
! declare real output
      real(kind=rk),intent(out),dimension(2)::det(2)!determinant
!-----------------------------------------------------------------------
! declare integer output
!      none
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! summary of call sequence arguments
!
!
!     on entry
!
!
!        luf       real(kind=rk),dimension(n, n)
!                  Packed LU factors of the permuted  matrix P*A
!                  (output of lufa);
!                  Packing is in the standard form;i.e. LUF = LF+UF-I,
!
!        info    integer
!                the singularity index flag returned by lufai
!                (0 for nonsingular, >0 for singularity).
!
!        ipvt    integer(n)
!                pivot vector returned by  lufa.
!
!        n       integer
!                the order of the matrix  A .
!
!     on return
!
!        det     real(kind=rk),dimension(2)
!                The determinant of the matrix A. det(1) containes the 
!                mantessa of the determinent and det(2) the exponent.
!                Thus: |A| = det(1) * 10.0**det(2)
!                with  1.0 .le. abs(det(1)) .lt. 10.0. Note that
!                det(1) will be set .eq. 0.0  for singular matrix.
!-----------------------------------------------------------------------
!     subroutines and functions used:  none,
!     module global variables used, one,ten,zero
!-----------------------------------------------------------------------
!     Declare internal variables
      integer i
!-----------------------------------------------------------------------
! First executable statement:
        if(info.ne.0)then;det=zero;return;endif!singular matrix
!-----------------------------------------------------------------------
!The variable small is defined as tiny(one)/epsilon(one)
        det(1) = one;det(2) =zero;!initialize det(*)
!-----------------------------------------------------------------------
        do    i = 1, n!loop over diagonal elements of matrix a
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = luf(i,i)*det(1)
            do; if( abs(det(1)) .ge. one)exit
            det(1) = ten*det(1);det(2) = det(2) - one
            enddo
            do; if( abs(det(1)) .lt. ten)exit
            det(1) = det(1)/ten;det(2) = det(2) + one
            enddo
        enddo
!-----------------------------------------------------------------------
      return
      end subroutine ludt


      subroutine lufa (a,n,luf,info,ipvt)
!-----------------------------------------------------------------------
!**** A general matrix A can be written as the produce of two factors;
!**** A = Lf*Uf,
!**** where Lf is  unit lower triangular and Uf is upper triangular.
!****
!**** lufa computes Lf and Uf and for permutated matrix A'= P*A,
!**** where P is a row permutation matrix chosen so that pivot
!**** elements of A' are as large as possible in absolute value.
!****
!**** Results are returned in the standard "packed" form;
!**** LUF=[Lf-I]+Uf (where I is the init matrix). 
!**** lufa also returns a flag, info= an index>0 if A is  singular and,
!**** ipvt= a pivot vector equivalent to the aformentioned permutation
!**** matrix P.
!     This version dated 11/22/04 .
!     Author= Micah Dembo of Boston University
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! declare real input
      real(kind=rk),intent(in),dimension(n,n)::a!an array
!-----------------------------------------------------------------------
! declare integer input
      integer,intent(in):: n!rank of matrix a
!-----------------------------------------------------------------------
! declare real output
      real(kind=rk),intent(out),dimension(n,n)::luf!lower factor
!-----------------------------------------------------------------------
! declare integer output
      integer,intent(out):: info!singularity indicator flag
      integer,intent(out),dimension(n)::ipvt!pivot vector
!-----------------------------------------------------------------------
! summary of call sequence arguments
!
!
!     On Entry
!
!        a       real(kind=rk),dimension(n, n)
!                the matrix, A, that is to be factored.
!
!        n       integer
!                the rank of the matrix (number of rows) in  A .
!
!     On Return
!
!
!        luf       real(kind=rk),dimension(n, n)
!                  Packed LU factors of the permuted  matrix P*A ;
!                  packing is in the standard form;i.e. LUF = LF+UF-I,
!                  where LF= lower triangular factor of P*A,
!                  where UF= upper triangular factor of P*A,
!                  and where I= unit matrix.
!
!        info    integer
!                a singularity indicator flag.
!                = 0  normal value if matrix non-singular.
!                = K  if  U(K,K) very small. This is not an error
!                     condition for this subroutine, but it does
!                     indicate that lusl will divide by zero.
!
!        ipvt    integer, dimension(n)
!                an integer vector of pivot indices that specify P.
!-----------------------------------------------------------------------
!     subroutines and functions used:  none,
!     module global variables used, one,small
!-----------------------------------------------------------------------
!     Declare internal variables
      real(kind=rk),dimension(n)::temp!automatic workspace
      real(kind=rk):: dmax,xmax!internal parameters and variables
      real(kind=rk):: t!an internal variable
      integer i,j,k,p,m!internal variables
!-----------------------------------------------------------------------
      info = 0!initialize singularty indicator.
!-----------------------------------------------------------------------
      luf=a!initialize uf=a
!-----------------------------------------------------------------------
      do k = 1, n-1!loop over columns 
         if (mod(k,100).eq.0) print *,'lufa: k=',k
!-----------------------------------------------------------------------
!     find p = pivot index and dmax=pivot score of column k
      temp=abs(luf(:,k))
      p=(k-1)+sum(maxloc(temp(k:)));dmax=temp(p)
      ipvt(k) = p!store pivot index of row k in the pivot vector
!-----------------------------------------------------------------------
      if(p.ne.k)then ! interchange rows p and k of LUF.
       temp(:)=luf(p,:);luf(p,:)=luf(k,:);luf(k,:)=temp(:)
      endif
!-----------------------------------------------------------------------
      if (dmax.lt. small)then;info=k;cycle;endif!zero pivot check
!-----------------------------------------------------------------------
!     apply multiplier for terms below diagonal of col-k
      luf(k+1:,k)=luf(k+1:,k)/luf(k,k)
!-----------------------------------------------------------------------
      do j = k+1,n!loop over remaining submatrix
       luf(k+1:,j)=luf(k+1:,j)-luf(k+1:,k)*luf(k,j)
      enddo!over  j 
!-----------------------------------------------------------------------
      enddo!over index k
!-----------------------------------------------------------------------
! test the lower right corner
      ipvt(n) = n;if (abs(luf(n,n)) .lt. small) info = n
!-----------------------------------------------------------------------
      return
      end subroutine lufa

      
      subroutine luiv(luf,info,ipvt,n,b)
!-----------------------------------------------------------------------
!     luiv computes the inverse of a matrix using the LU factors
!     ,the pivot vector, and the info flag returned by lufa.
!     This version dated 11/22/04 .
!     Author= Micah Dembo of Boston University
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! declare real input
      real(kind=rk),intent(in),dimension(n,n)::luf!luf factor matrix
!-----------------------------------------------------------------------
! declare integer input
      integer,intent(in):: info!singularity flag from lufa
      integer,intent(in),dimension(n)::ipvt!pivot vector from lufa
      integer,intent(in):: n! number of erautions in linear system
!-----------------------------------------------------------------------
! declare real output
      real(kind=rk),intent(out),dimension(n,n)::b!=inverse(lf*uf)
!-----------------------------------------------------------------------
! declare integer output
!      none
!-----------------------------------------------------------------------
!
!     on entry
!
!
!        luf       real(kind=rk),dimension(n, n)
!                  Packed LU factors of matrix P*A (output of lufa)
!                  LUF= (LF-I)+UF .
!
!        info    integer
!                the singularity indicator flag (output from lufa).
!
!        ipvt    integer, dimension(n)
!                an integer vector of pivot indices which incodes the
!                permutation matrix in an efficent format.
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        b       real(kind=rk),dimension(N, N)
!                Inverse of the matrix lf*uf if info=0;
!                otherwise , lf*uf is singular) so b is set equal to 
!                the zero matrix as a default outcome.
!-----------------------------------------------------------------------
!     subroutines and functions used:  none,
!     module global variables used, small,one,ten,zero
!-----------------------------------------------------------------------
!declare internal variables
      real(kind=rk),dimension(n)::temp!automatic working storage
      real(kind=rk)::t!automatic working storage
      integer i,k,j,p!internal loop variables
!-----------------------------------------------------------------------
      b=zero
!-----------------------------------------------------------------------
      if(info.ne.0)return! since inverse(a) dosn't exist
!-----------------------------------------------------------------------
!     control reaches here only if inverse[A] exists.
!-----------------------------------------------------------------------
      do i=1,n;b(i,i)=one;enddo !set diagonal elements of B = one
      do  k = 1, n!loop over rows of matrix B
      p=ipvt(k)!get pivot index of row k
      if(p.ne.k)then !Interchange rows p and k of B .
      temp(:)=b(p,:);b(p,:)=b(k,:);b(k,:)=temp(:)
      endif
      enddo
!      the  right side matrix B is now the permutation matrix, P
!-----------------------------------------------------------------------
!     FIRST SOLVE  LF*Y = P using forward substitution for each column
      do  k = 1, n
      b(k,:)=b(k,:)!*(one/lf(k,k))
      do j=k+1,n; b(j,:)=b(j,:)-b(k,:)*luf(j,k) ;enddo
      enddo!over k
!-----------------------------------------------------------------------
!      NOW SOLVE  UF*X = Y using back substitution for each column.
       do  k = n,1,-1
       b(k,:) = b(k,:)*(one/luf(k,k))
       do j=1,k-1;b(j,:)=b(j,:)-b(k,:)*luf(j,k);enddo
       enddo!over k
!-----------------------------------------------------------------------
! B now contains the inverse of A
!-----------------------------------------------------------------------
      return
      end subroutine luiv


      subroutine lup1(a,n,anorm)
!-----------------------------------------------------------------------
!     lup1 calculates the P-1 norm of a matrix.
!     It works directly from the definition;
!                norm(A)==maxsum of abs value colums of A.
!     Here norm(*) signifies the 1-norm of the argument matrix. 
!     This version dated 11/22/04 .
!     Author= Micah Dembo of Boston University
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! declare real input
      real(kind=rk),intent(in),dimension(n,n)::a! matrix A
!-----------------------------------------------------------------------
! declare integer input
      integer,intent(in):: n!array dimensions
!-----------------------------------------------------------------------
! declare real output
      real(kind=rk),intent(out)::anorm!the p1 norm of A
!-----------------------------------------------------------------------
! declare integer output
!      none
!-----------------------------------------------------------------------
! summary of call sequence arguments
!
!
!     On Entry
!
!
!        a       real(kind=rk),DIMENSION(N, N)
!                Matrix for which the norm will
!                be  calculated.
!
!        n       integer
!                the order ofsquare matrix AB .
!
!
!     on return
!
!        anorm   real(kind=rk)
!                The P-1 norm (max absolute row sum) of  A .
!-----------------------------------------------------------------------
!     subroutines and functions used:  none,
!     module global variables used, one,zero
!-----------------------------------------------------------------------
!     Declare internal variables
      integer::k
      real(kind=rk),dimension(n):: acol!storage for column sums
!-----------------------------------------------------------------------
      do k=1,n;acol(k)=sum(abs(a(:,k)));enddo!load column sums
      anorm=maxval(acol)!take maximum column sum
!-----------------------------------------------------------------------
      return
      end subroutine lup1

      subroutine lusl (luf,r,info,ipvt,n,x)
!-----------------------------------------------------------------------
!     lusl solves the real system A*x=r  using the LU factorization
!     of the permuted matrix P*A and the associated info flag and
!     pivot vector returned by lufa.
!     This version dated 11/22/04 .
!     Author= Micah Dembo of Boston University
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! declare real input
      real(kind=rk),intent(in),dimension(n,n)::luf!packed lU factors
      real(kind=rk),intent(in),dimension(n)::r!right side vec
!-----------------------------------------------------------------------
! declare integer input
      integer,intent(in):: info!info flag from lufa
      integer,intent(in),dimension(n)::ipvt!pivot vector
      integer,intent(in):: n!array rank
!-----------------------------------------------------------------------
! declare real output
      real(kind=rk),intent(out),dimension(n):: x!solution vec(out)
!-----------------------------------------------------------------------
! declare integer output
!      none
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! summary of call sequence arguments
!
!     On Entry
!
!
!        luf       real(kind=rk),dimension(n, n)
!                  Packed LU factors of matrix P*A (output of lufa)
!                  LUF= (LF-I)+UF .
!
!        r       real(kind=rk) DIMENSION(N)
!                the right hand side vector of the system Ax=r.
!
!        info    integer
!                the singularity index flag (output from lufa).
!
!        ipvt    integer, dimension(n)
!                an integer vector of pivot indices which encode the
!                permutation matrix P in an efficent storage format.
!
!        n       INTEGER
!                the order of the matrix  a .
!
!     On Return
!
!        x       real(kind=rk) DIMENSION(N)
!                the solution vector x if info=0, else x =zero vector.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     subroutines and functions used:  none,
!     module global variables used: zero
!-----------------------------------------------------------------------
      real(kind=rk):: temp!internal variable
      integer k,p,j!internal variables
!-----------------------------------------------------------------------
      x=zero!initialize x to default solution
!-----------------------------------------------------------------------
      if(info.ne.0)return!the matrix is singluar so x=default solution
!-----------------------------------------------------------------------
!     Control reaches here only if the matrix is non-singular
!-----------------------------------------------------------------------
      x=r!initialize x = right-side vector
      do  k = 1, n !permute x using the pivot vector
      p= ipvt(k) ; temp = x(p) ; x(p) = x(k);x(k)=temp
      enddo
!     At this stage, x = P*r
!-----------------------------------------------------------------------
!     FIRST SOLVE  LF*x = P*r by forward substitution
      do  k = 1, n
      x(k+1:)=x(k+1:) - luf(k+1:,k)*x(k) 
      enddo
!-----------------------------------------------------------------------
!
!        NOW SOLVE  UF*x = y by back substitution
!
       do  k = n,1,-1;x(k) = x(k)/luf(k,k)
       x(:k-1)=x(:k-1) - luf(:k-1,k)*x(k) 
       enddo
      return
      end subroutine lusl

      end module lupack
!*********************END LUPACK ROUTINES*******************************
