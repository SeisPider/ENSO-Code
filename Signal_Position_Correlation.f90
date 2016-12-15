Subroutine Signal_Corr(Direction_in,Direction_out)
    Use math_tools
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 本部分代码旨在提取各个经纬度点上时变序列与提取出的ENSO信号间的相关场                    !
    ! 输入文件: 各异常场的月均值文件                                                          !
	!                路径在 E:\ONI-Pressure_Filter\Monthly_Anomaly\                             !                                  !
    ! 输出文件：各高度层上ENSO信号与原始异常场的相关场及其延迟信息                            !
	!                路径在 C:\Users\Seispider\Documents\Vertical_structure\Longitude-Height  !                                                        !
    ! 相关方法： Pearson相关                                                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                    变量定义区                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer,parameter   ::rk=4,ik=4,Pressure_num=40,Month_length=121,Month_Signal=97,Position_num=2701 &
	                      ,longitude_num=73,latitude_num=37
	Character,intent(in):: Direction_in*200,Direction_out*200
    Real*4              :: S_Co_position(position_num),S_Del_position(position_num) &
	                       ,T_Co_position(position_num),T_Del_position(position_num)
	INTEGER             :: length_vector_A,length_vector_B,longitude_count,latitude_count
	REAL                :: Coved_ecof(position_num,Month_Signal)
	CHARACTER           :: Pressure*4, Physical*1
	Character*4         :: Pressure_matrix(Pressure_num)
	Real                :: S_Signal(Month_Signal),T_Signal(Month_Signal) &
	                       ,S_Anoma(Position_num,Month_Signal),T_Anoma(Position_num,Month_Signal)
	Character           :: Raw_Month_day*14
    Character*400       :: S_File_Direction,T_File_Direction   !设置输入,输出文件路径
	Character           :: Head*111
	Real                :: S_f_up,T_f_up
	Integer             :: Pressure_count,Month_count,Position_count
	Logical             :: alive
	Data Pressure_matrix /"30","50","70","100","125",&
                       	"150","175","200","225","250",&
						"275","300","325","350","375",&
						"400","425","450","475","500",&
						"525","550","575","600","625",&
						"650","675","700","725","750",&
						"775","800","825","850","875",&
						"900","925","950","975","1000"/
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                         载入ENSO信号,ViP_information,异常场月均值文件                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do Pressure_count=1,Pressure_num
	S_File_Direction=trim(Direction_in)//"Glo_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\SMode_variation.txt"
	T_File_Direction=trim(Direction_in)//"Glo_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\TMode_variation.txt"
	Inquire(file=trim(S_File_Direction),exist=alive)
	if(.not.alive)then
	    cycle
	else
	Open(11,file=trim(S_File_Direction))
    Open(12,file=trim(T_File_Direction))
	Do Month_count=1,Month_Signal
	   Read(11,"(A14,f18.8)")Raw_Month_day,S_Signal(Month_count)
	   Read(12,"(A14,f18.8)")Raw_Month_day,T_Signal(Month_count)
	Enddo
	close(11);close(12);
	Open(13,file=trim(Direction_in)//"Glo_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\SVIP_information.txt")
	Open(14,file=trim(Direction_in)//"Glo_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\TVIP_information.txt")
	Read(13,"(A111,f16.6)")Head,S_f_up
	Read(14,"(A111,f16.6)")Head,T_f_up
	close(13);close(14);
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                     载入异常场月均值文件                                !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(15,file=trim(Direction_in)//"Monthly_Anomaly\SMonth_"//trim(Pressure_matrix(Pressure_count))//".txt")
    Open(16,file=trim(Direction_in)//"Monthly_Anomaly\TMonth_"//trim(Pressure_matrix(Pressure_count))//".txt")
	Do Position_count=1,Position_num
	   Read(15,*)(S_Anoma(Position_count,Month_count),Month_count=1,Month_Signal)
	   Read(16,*)(T_Anoma(Position_count,Month_count),Month_count=1,Month_Signal)
	Enddo
	close(15);close(16);
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!            对数据做滤波处理（频率为相关系数使ENSO信号与ONI达到最大的滤波频率）          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call bandpass_filter(Month_Signal,Month_Signal,Position_num,Coved_ecof,S_Anoma,S_f_up)
	S_Anoma=Coved_ecof
    call bandpass_filter(Month_Signal,Month_Signal,Position_num,Coved_ecof,T_Anoma,T_f_up)
	T_Anoma=Coved_ecof
	call Signal_Position_Correlation(S_Anoma,position_num,S_Signal,Month_Signal,S_Co_position,S_Del_position)
	call Signal_Position_Correlation(T_Anoma,position_num,T_Signal,Month_Signal,T_Co_position,T_Del_position)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                   输出文件至指定目录                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(17,file=trim(Direction_out)//"S_CP_"//trim(Pressure_matrix(Pressure_count))//"_Cor.txt")
	Open(18,file=trim(Direction_out)//"S_CP_"//trim(Pressure_matrix(Pressure_count))//"_Del.txt")
    Open(19,file=trim(Direction_out)//"T_CP_"//trim(Pressure_matrix(Pressure_count))//"_Cor.txt")
	Open(20,file=trim(Direction_out)//"T_CP_"//trim(Pressure_matrix(Pressure_count))//"_Del.txt")
	Do latitude_count=1,latitude_num
	   Do longitude_count=1,longitude_num
	   Position_count=(latitude_count-1)*longitude_num+longitude_count
	   write(17,"(2I8,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),S_Co_position(Position_count)
	   write(18,"(2I8,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),S_Del_position(Position_count)
	   write(19,"(2I8,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),T_Co_position(Position_count)
	   write(20,"(2I8,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),T_Del_position(Position_count)
	   Enddo
	Enddo
	close(17);close(18);close(19);close(20);
	Endif
	
	Enddo
End

Subroutine Signal_Position_Correlation(x,position,Signal,n,Co_position,Del_position)
    USE math_tools
    Implicit none
    Integer*4,intent(in)   ::n,position
    Real*4                 ::vectorin(n),vectorout(n),Co_vector(n)
    Real*4                 ::x(position,n),Signal(n)
    Integer                ::mode_count,lenth,max_mode,position_count,latitude_count,longitude_count,month_count,t
    Real                   ::One_month_data(n),One_month_smoothed(n)
    Character              ::Time_C*14,Physical*1
    Real*4                 ::Co_position(position),Del_position(position),ONI_index(121),ONI(n)
    Real                   ::Co_t(25),Del_t(25),vector,Co_max,Del_max,Co_r
    Integer                ::month_max
    
    !********************************************************************************************************
    !本部分代码旨在对每一点上的时间序列做1-2-1二项式滑动平均
    !********************************************************************************************************
    Month_max=n
    Do position_count=1,position
       One_month_data(:)=x(position_count,:)
       call running2(One_month_data,3,One_month_smoothed,month_max)
       x(position_count,:)=One_month_smoothed(:)
    Enddo
    
    !******************************************************************************************************
    !下面旨在对上述滑动平均的结果与ONI做相关
    !******************************************************************************************************
    Do position_count=1,position
        Co_vector(:)=x(position_count,:)
            
        Do t=1,25
    !******************************************************************************************************
	!下面部分给出t=1:13时的相关结果，表示Signal提前于时间序列
	!******************************************************************************************************
		   if(t<=13)then
		   call myPearson(Co_vector(1:(n-t+1)),Signal(t:n),Co_r,(n-t+1))
		   Co_t(t)=Co_r
           Del_t(t)=t-1 !说明Signal延迟于Time serries 的月份数
		   else
		   call myPearson(Co_vector((t-12):n),Signal(1:(n-(t-12)+1)),Co_r,(n-(t-12)+1))
		   Co_t(t)=Co_r
           Del_t(t)=-(t-13) !说明Signal延迟于Time serries 的月份数
		   Endif
		   
        Enddo
        
        Co_max=0;Del_max=0;
        Do t=1,25
            if(abs(Co_t(t)).ge.abs(Co_max))then
                Co_max=Co_t(t);
                Del_max=Del_t(t);
            Endif
        Enddo
        Co_position(position_count)=Co_max;
        Del_position(position_count)=Del_max;
       
    Enddo    
    End Subroutine
