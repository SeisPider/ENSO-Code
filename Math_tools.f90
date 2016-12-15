module math_tools
    
    contains
       
	Subroutine run_mean(datain,dataout,lenth,smonum)
    Implicit none
    Real,intent(in)::datain(:)
    Real,intent(out)::dataout(:)
    Integer,intent(in)::lenth,smonum
    
    Real::sum
    Integer::index,smo_index_min,smo_index_max,smo_index,smonum_half
    Do index=1,lenth
            smonum_half=floor(real(smonum)/2)
            smo_index_min=index-smonum_half
            smo_index_max=index-smonum_half+smonum-1
        if(smo_index_min.ge.1.and.smo_index_max.le.lenth)then
           sum=0
            Do smo_index=smo_index_min,smo_index_max
                sum=sum+datain(smo_index)
            Enddo
            dataout(index)=sum/smonum
        else
            dataout(index)=datain(index)
        Endif
    EndDo           
    End Subroutine
	
	! 求逆矩阵
subroutine inverse(A,IA)
  implicit none
  real    :: A(:,:), IA(:,:)
  real, allocatable :: B(:,:)
  integer :: i,j,N
  N = size(A,1)  
  allocate(B(N,N))
  ! 先把IA设定成单位矩阵
  forall(i=1:N,j=1:N,i==j) IA(i,j)=1.0
  forall(i=1:N,j=1:N,i/=j) IA(i,j)=0.0
  ! 保存原先的矩阵A, 使用B来计算
  B=A 
  ! 把B化成对角线矩阵(除了对角线外,都为0)
  call Upper(B,IA,N) ! 先把B化成上三角矩阵
  call Lower(B,IA,N) ! 再把B化成下三角矩阵
  ! 求解
  forall(i=1:N) IA(i,:)=IA(i,:)/B(i,i) 
  return
end subroutine
! 输出矩阵的子程序
subroutine output(matrix)
  implicit none
  real    :: matrix(:,:)
  integer :: m,n,i
  character(len=20) :: for='(??(1x,f6.3))'
  m = size(matrix,1)
  n = size(matrix,2)
  ! 用字符串来设定输出格式
  write( FOR(2:3), '(I2)' ) N
  do i=1,N
	write( *, FMT=FOR ) matrix(i,:)
  end do
  return
end subroutine output
! 求上三角矩阵的子程序
subroutine Upper(M,S,N)
  implicit none
  integer :: N
  real    :: M(N,N)
  real    :: S(N,N)
  integer :: I,J
  real :: E
  do I=1,N-1
    do J=I+1,N              
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine Upper
! 求下三角矩阵的子程序
subroutine Lower(M,S,N)
  implicit none
  integer :: N
  real    :: M(N,N)
  real    :: S(N,N)
  integer :: I,J
  real :: E
  do I=N,2,-1
    do J=I-1,1,-1           
      E=M(J,I)/M(I,I)
      M(J,1:N)=M(J,1:N)-M(I,1:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine Lower

Subroutine  vector_standard(vectorin,vectorout,lenth)
!USE Data_module
Implicit none
Integer,parameter  ::rk=4,ik=4
Real(rk),intent(in)::vectorin(:)
Real(rk),intent(out)::vectorout(:)
Integer(ik)::lenth
Real(rk)::sum,num,mean_spec,Std_spec
Integer(ik)::month_min=1,month_max,month_count
       month_max=lenth
	   sum=0;num=0;
	   Do month_count=month_min,month_max
	   sum=sum+vectorin(month_count);
	   num=num+1;
	   Enddo
	   mean_spec=sum/num;
	   
	   sum=0;num=0;
	   Do month_count=month_min,month_max
	   sum=sum+(vectorin(month_count)-mean_spec)**2
	   num=num+1;
	   Enddo
	   Std_spec=sqrt(sum/(num-1));
	   
	   Do month_count=month_min,month_max
	   vectorout(month_count)=(vectorin(month_count)-mean_spec)/Std_spec;
       Enddo
	   
End Subroutine vector_standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!找时间将皮尔逊相关的子程序
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine myPearson(vectorA,vectorB,Co_r,lenth)

Implicit none
Real,intent(in)::vectorA(:),vectorB(:)
Real,intent(out)::Co_r
Integer,intent(in)::lenth
Integer::count
Real::sumA2,sumB2,sumA,sumB,sumAB
sumA=0;sumB=0;
Do count=1,lenth
    sumA=sumA+vectorA(count)
    sumB=sumB+vectorB(count)
Enddo
sumA2=0;sumB2=0;sumAB=0;
Do count=1,lenth
    sumA2=sumA2+vectorA(count)**2
    sumB2=sumB2+vectorB(count)**2
    sumAB=sumAB+vectorA(count)*vectorB(count)
Enddo
Co_r=(lenth*sumAB-sumA*sumB)/sqrt(lenth*sumA2-sumA**2)/sqrt(lenth*sumB2-sumB**2)
End Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!本部分程序用以求取标准差
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Standard_deviation(vectorin,STD,length)
!USE Data_module

Implicit none
Integer,parameter  ::rk=4,ik=4
Integer(ik),intent(in)::length
Real(rk),intent(in)::vectorin(length)
Real(rk)::sum,mean
Real(rk),intent(out)::STD
Integer(ik)::length_count,num
  sum=0;num=0;mean=0;
  Do length_count=1,length
     sum=sum+vectorin(length_count)
     num=num+1
  Enddo
     mean=sum/num
     sum=0;num=0;
  Do length_count=1,length
     sum=sum+(vectorin(length_count)-mean)**2
     num=num+1
  Enddo
  STD=sqrt(sum/(num-1))
End Subroutine

    End module math_tools
	