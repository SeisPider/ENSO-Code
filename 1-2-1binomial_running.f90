SUBROUTINE running2(aa, kkk, bb, key)
!!                                                                                                                                    
!! Porpose    : binomial running.                                              
!!                                                                                                    
!!                                                                                                                                    
!! Input      : aa, original, 1-dimension, aa(key)
!!              kkk, length of running window. kkk =3 即是1-2-1二项式平滑
!!              key, length of aa                                                                                                                   
!!                                                                                                                                    
!! Output     : bb, smoothed data, 1-dimension, bb(key)                                                                                                                 
!!                                                                                                                                    
!! Restriction: ----                                                                                                                  
!!                                                                                                                                    
!! SR Call    : ----                                                                                                               
!!                                                                                                                                    
!! Remarks    : ----                                                                                                                  
!!                                                                                                                                    
!! Author     :                                                                                                
!!                                                                                                                                    
!! Version    :                                                                                               
!!                                                                                                                                    
!! Created    : 2016/08/29                                                                                                           
!!                                                                                                                                    
!! Modified   : 2016/08/29  Zhiping CHEN                                                                                             
!!                Created. 
!!              2016/08/29 Zhiping CHEN [valid]                                                                                                  
!!                                                                                                                                    
!! Fund       :                                                                                 
!!                                                                                                                                    
!! Copyright  :                                                                                              
!!                                                                                                                 
!!                                                                                                                                    
!* 
IMPLICIT NONE
  INTEGER  :: key, kkk, k, kw, k2, j, k0
  REAL*4   :: aa(key), bb(key)
  REAL*4   :: a(0:kkk*2+1), b(0:kkk*2+1)
  REAL*4   :: xx, yy, yb
  kw = kkk*2+1     !! 
  !!write(*,*) kw  !!; pause 5
  DO k = 0, kw
    a(k) = 0
	b(k) = 0
  END DO
  a(1) = 1
  b(1) = 1
  DO k = 1, kw
    k2 = (k+1)/2
	yy = 0
	DO j = 1, k
	  b(j) = a(j) + a(j-1)
	
	END DO
	
	DO j = 1, kw
	  a(j) = b(j)
	  b(j) = b(j)/(2**(k-1))
	  yy   = yy + b(j)
	END DO
	IF (MOD(k, 2).NE.0) THEN
	  xx = 0
	  DO k0 = 1, k
	    xx = xx + aa(k0)*b(k0)/yy
	  END DO
	  bb(k2) = xx
	  !!write(*,*) 'k2', k2, bb(k2)
	ENDIF
  END DO
  yb = yy
  DO k = kkk + 2, key - (kkk+1)
    xx = 0
	j  = 1
	DO k0 = k - kkk, k + kkk
	  xx = xx + aa(k0)*b(j)/yb
	  j  = j  + 1
	END DO
	bb(k) = xx
	!!write(*,*) k, bb(k)
   END DO
   DO k = 0, kw
     a(k) = 0
     b(k) = 0
   END DO
   a(1) = 1
   b(1) = 1
   DO k = 1, kw
     k2 = key - k/2
	 yy = 0
	 DO j = 1,k
	   b(j) = a(j) + a(j - 1)
	 END DO
	 DO j = 1, kw
	   a(j) = b(j)
	   b(j) = b(j)/(2**(k-1))
	   yy = yy + b(j)
	 END DO

    !! pause 15
	 IF (MOD(k,2) .NE. 0) THEN   !!若K能被2整除，则
	   xx = 0
	   DO k0 = 1, k
	     xx = xx + aa(key+1-k0)*b(k0)/yy
	   END DO
	   bb (k2) = xx
	  !! print*, 'k2', k2, bb(k2)
	 END IF
   END DO
   RETURN
 END 