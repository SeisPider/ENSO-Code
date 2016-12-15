Subroutine Signal_Correlation(x,n,Direction,Co_t,Del_t)
    USE math_tools
    Implicit none
    Character*200          ::Direction
    Integer*4,intent(in)   ::n
    Real*4                 ::vectorin(n),vectorout(n),Co_vector(n)
    Real*4                 ::x(n)
    Integer                ::mode_count,lenth,max_mode,position_count,latitude_count,longitude_count,month_count,t
    Real                   ::One_month_data(n),One_month_smoothed(n)
    Character              ::Time_C*14,Physical*1
    Real*4                 ::ONI_index(121),ONI(n)
    Real                   ::Co_t(25),Del_t(25),vector,Co_max,Del_max,Co_r
    Integer                ::month_max
    
    !********************************************************************************************************
    !本部分代码旨在载入ONI时间序列
    !********************************************************************************************************
    month_max=n
    Open(12,file=trim(Direction)//"\ONI_index.txt")
    Do month_count=1,121
        Read(12,*)ONI_index(month_count)
    Enddo
    close(12)
    !********************************************************************************************************
    !本部分代码旨在对每一点上的时间序列做1-2-1二项式滑动平均
    !********************************************************************************************************
    
    !  One_month_data(:)=x(:)
    !   call running2(One_month_data,3,One_month_smoothed,month_max)
    !   x(:)=One_month_smoothed(:)
    
    !******************************************************************************************************
    !下面旨在对上述滑动平均的结果与ONI做相关
    !******************************************************************************************************
        Co_vector(:)=x(:)    
        Do t=1,25
            ONI(:)=ONI_index(t:(t+month_max-1))
            call myPearson(Co_vector,ONI,Co_r,month_max)
            Co_t(t)=Co_r
            Del_t(t)=-13+t
        Enddo    
  
    End Subroutine