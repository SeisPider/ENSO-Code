Subroutine Signal_correlation_lag(ecof,Pressure,n,mnl,Corr_max_abs,f_up,Physical,max_mode_abs,Del_max_abs,Direction,Corr_max,max_mode,Del_max,Corr_min,min_mode,Del_min)
    USE math_tools
    Implicit none
    Integer                ::month_max
    Character*200          ::Direction
    Integer*4              ::n,mnl
    Real*4                 ::ecof(mnl,n),f_up
    Character              ::Pressure*4,Physical*1
    Real*4                 ::vectorin(n),vectorout(n),Co_vector(n)
    Integer                ::mode_count,lenth,month_count,mnl_count,t
    Real                   ::One_month_data(n),One_month_smoothed(n),ONI_index(121),ONI(n),Co_r
    Character              ::Time_C*14

    Real::Co_t(25),Del_t(25),Co_max,Co_mode(mnl),Del_mode(mnl)
    
    Integer::max_mode_abs,min_mode,max_mode
    Real   ::Corr_max_abs,Del_max_abs,Del_max,Corr_max,Corr_min,Del_min
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
    !本部分代码旨在对时间特征值向量做1-2-1二项式滑动平均
    !********************************************************************************************************
    
    Do mnl_count=1,mnl
       One_month_data(:)=ecof(mnl_count,:)
       call running2(One_month_data,3,One_month_smoothed,month_max)
       ecof(mnl_count,:)=One_month_smoothed(:)
    Enddo
    
    !******************************************************************************************************
    !下面旨在对上述滑动平均的结果与ONI做相关
    !******************************************************************************************************
    Do mode_count=1,mnl
        Co_vector(:)=ecof(mode_count,:)
            
         Do t=1,25
             ONI(:)=ONI_index(t:(t+month_max-1))
             call myPearson(Co_vector,ONI,Co_r,month_max)
             Co_t(t)=Co_r
             Del_t(t)=-13+t
         Enddo
        
        Co_max=0;Del_max=0;
        Do t=1,25
            if(abs(Co_t(t)).ge.abs(Co_max))then
                Co_max=Co_t(t);
                Del_max=Del_t(t);
            Endif
        Enddo
        Co_mode(mode_count)=Co_max;
        Del_mode(mode_count)=Del_max; 
    Enddo 
   
    Corr_max_abs=0;Del_max_abs=0;max_mode_abs=0;
	Corr_max=0;Del_max=0;max_mode=0;
	Corr_min=0;Del_min=0;min_mode=0;
    Do mode_count=1,mnl
        if((abs(Co_mode(mode_count)).ge.abs(Corr_max_abs)).and.(abs(Co_mode(mode_count)).le.1))then
		max_mode_abs=mode_count;
		Corr_max_abs=Co_mode(mode_count);
		Del_max_abs=Del_mode(mode_count);
		Endif
		if((Co_mode(mode_count).ge.Corr_max).and.(abs(Co_mode(mode_count)).le.1))then
		max_mode=mode_count;
		Corr_max=Co_mode(mode_count);
		Del_max=Del_mode(mode_count);
		Endif
		if((Co_mode(mode_count).le.Corr_min).and.(abs(Co_mode(mode_count)).le.1))then
		min_mode=mode_count;
		Corr_min=Co_mode(mode_count);
		Del_min=Del_mode(mode_count);
		Endif
    Enddo
	
    End Subroutine