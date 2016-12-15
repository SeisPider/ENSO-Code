Subroutine Latitude_mean_eddy(Direction_in,Direction_out,region)
    Use math_tools
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! �����ִ���ּ����ȡ����γ�ȵ㡢�߶ȵ��ϵ�γ���ֵ����������������EOF�ֽ⣬�����ָ���� !
	!           �ȷ�Χ�ϵ�Eddy��ֵ��Ϣ                                                        !                                                                                          !
    ! �����ļ�: ԭʼ�쳣�����ݵ��¾�ֵ�ļ�                                                    !
	!                ·���� E:\ONI-Pressure_????\Monthly_Anomaly\                             !
    ! ����ļ���ȫ��߶ȵ�γ���ֵ��߶ȵľ����ļ���EOF�ֽ��ļ�                               !
	!                ·���� C:\Users\Seispider\Documents\Vertical_structure\latitude_mean_EOF !                                                        !
    ! ��ط�����Pearson��ء�EOF�ֽ�                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                    ����������                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer,parameter   ::rk=4,ik=4,Pressure_num=40,Month_length=121,Month_Signal=97,Position_num=2701 &
	                      ,longitude_num=73,latitude_num=37
    Character,intent(in):: Direction_in*200,region*4
	
    Character*14        :: Time_series(Month_Signal)
	Integer             :: Pressure_count,Real_Pressure_num=0
	Integer             :: Position_count,Month_count,S_num,T_num
	Integer             :: latitude_min,latitude_max,longitude_min,longitude_max,longitude_count,&
	                       latitude_count,longitude_ALL,latitude_ALL,latitude_mean_position
	Real                :: Vector_standard_in(longitude_num),Vector_standard_out(longitude_num),Vector_height_in(Month_Signal),Vector_height_out(Month_Signal)					   
	Real                :: S_Month_Anoma(Pressure_num,Position_num,Month_Signal),T_Month_Anoma(Pressure_num,Position_num,Month_Signal)
	Real                :: S_latitude_Average(Pressure_num,latitude_num,Month_Signal),&
	                       T_latitude_Average(Pressure_num,latitude_num,Month_Signal),&
						   STD_in(Month_Signal),STD_out, &
						   S_longitude_STD(Pressure_num,longitude_num),T_longitude_STD(Pressure_num,longitude_num), &
						   S_Time_height_average(Pressure_num,Month_Signal),T_Time_height_average(Pressure_num,Month_Signal),&
						   S_latitude_mean((Position_num*latitude_num),Month_Signal),T_latitude_mean((Position_num*latitude_num),Month_Signal)
	Real                :: S_sum,T_sum,P,P0=1013.15!unit=hpa standard atmospheric pressure
	Logical             :: alive
	Character*200       :: Direction_out,S_File_Direction,T_File_Direction,ONI_direction   !��������,����ļ�·��
	
	Character*4         :: Pressure_matrix(Pressure_num)
	Data Pressure_matrix /"30","50","70","100","125",&
                       	"150","175","200","225","250",&
						"275","300","325","350","375",&
						"400","425","450","475","500",&
						"525","550","575","600","625",&
						"650","675","700","725","750",&
						"775","800","825","850","875",&
						"900","925","950","975","1000"/
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                         �����ļ�                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Real_Pressure_num=0;
    Do Pressure_count=1,Pressure_num
	   S_File_Direction=trim(Direction_in)//"SMonth_"//trim(Pressure_matrix(Pressure_count))//".txt"
	   T_File_Direction=trim(Direction_in)//"TMonth_"//trim(Pressure_matrix(Pressure_count))//".txt"
	   inquire(file=trim(S_File_Direction),exist=alive)
	   if(.not.alive)then
             cycle
       else
	   Real_Pressure_num=Real_Pressure_num+1
	   Open(11,file=trim(S_File_Direction))
	   Open(12,file=trim(T_File_Direction))
	   Do Position_count=1,Position_num
	      Read(11,*)(S_Month_Anoma(Pressure_count,Position_count,Month_count),Month_count=1,Month_Signal)
		  Read(12,*)(T_Month_Anoma(Pressure_count,Position_count,Month_count),Month_count=1,Month_Signal)
	   Enddo
       close(11);close(12);
	   Endif
	Enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                           ���׼��                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(trim(region).eq."Glo")then
	   latitude_min=1;latitude_max=37;longitude_min=1;longitude_max=73;
	elseif(trim(region).eq."Mid")then
	   latitude_min=7;latitude_max=31;longitude_min=1;longitude_max=73;
    elseif(trim(region).eq."Equ")then
	   latitude_min=13;latitude_max=25;longitude_min=1;longitude_max=73;
	elseif(trim(region).eq."Five")then
	   latitude_min=18;latitude_max=20;longitude_min=1;longitude_max=73;
	elseif(trim(region).eq."Nino")then
	   latitude_min=18;latitude_max=20;longitude_min=3;longitude_max=13;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !�˴�����������ѡ�����γɶ�ȫ������������variation��STD                                !
    !   region���û������������������ơ�                                                     !
    !   longitude_min:��ʾ���������ľ����������е��±���[-180 => 1;-175 => 2;-170 => 3;..] !
    !   longitude_max:��ʾ�������Ҳ�ľ����������е��±���[-180 => 1;-175 => 2;-170 => 3;..] !
    !   latitude_min :��ʾ�������²��γ���������е��±���[-90 => 1;-85 => 2;-80 => 3;..]    !
    !   latitude_max :��ʾ�������ϲ��γ���������е��±���[-90 => 1;-85 => 2;-80 => 3;..]    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  ��γ���ϵ�����������ά                                !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	latitude_mean_position=0;
	Do Pressure_count=1,Real_Pressure_num
	Do latitude_count=1,latitude_num
	   Do Month_count=1,Month_Signal
	   S_sum=0;T_sum=0;S_num=0;T_num=0;
	   Do longitude_count=1,longitude_num
	      Position_count=(latitude_count-1)*longitude_num+longitude_count
		  S_sum=S_sum+S_Month_Anoma(Pressure_count,Position_count,Month_count)
		  T_sum=T_sum+T_Month_Anoma(Pressure_count,Position_count,Month_count)
		  S_num=S_num+1;T_num=T_num+1;
	   Enddo
	   S_latitude_Average(Pressure_count,latitude_count,Month_count)=S_sum/S_num
	   T_latitude_Average(Pressure_count,latitude_count,Month_count)=T_sum/T_num
	   Enddo
	   latitude_mean_position=latitude_mean_position+1
	   Vector_height_in(:)=S_latitude_Average(Pressure_count,latitude_count,:)
	   call vector_standard(Vector_height_in,Vector_height_out,Month_Signal)
	   S_latitude_mean(latitude_mean_position,:)=Vector_height_out(:)	 

	   Vector_height_in(:)=T_latitude_Average(Pressure_count,latitude_count,:)
	   call vector_standard(Vector_height_in,Vector_height_out,Month_Signal)
	   T_latitude_mean(latitude_mean_position,:)=Vector_height_out(:)	   
	Enddo
	Enddo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                           ���γ���ֵ�������ļ�latitude_mean                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(13,file=trim(Direction_out)//"S_latitude_mean.txt")
	Open(14,file=trim(Direction_out)//"T_latitude_mean.txt")
	Do Position_count=1,latitude_mean_position
	   write(13,"<Month_Signal>f18.8")(S_latitude_mean(Position_count,Month_count),Month_count=1,Month_Signal)
	   write(14,"<Month_Signal>f18.8")(T_latitude_mean(Position_count,Month_count),Month_count=1,Month_Signal)
	Enddo
	close(13);close(14);
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!         ����ͨ������Filter_EOF_Correlationģ���γ���ֵ��������ͨ�˲���EOF�ֽ�        !                            
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    
    
    
    
    
    
    
    
   
    
End Subroutine

Subroutine  Filter_EOF_Correlation(Field,Position_num,Month_Signal,Direction_out)
    USE math_tools  
    Implicit none
    Integer,parameter     ::rk=4,ik=4
    Integer,parameter     ::latitude_max=37,longitude_max=73
    Integer(ik),parameter ::Fre_grid=8
    Real(rk),parameter    ::Fre_min=0.06,Fre_max=0.083
    Real(rk)              ::Coved_ecof(position_num,Month_Signal)
    Integer,parameter     ::m=Position_num,n=Month_Signal,mnl=Month_Signal,ks=1
    real(rk)              ::Field(Position_num,Month_Signal),f_up
    real                  ::x(m,n),egvt(m,mnl),ecof(mnl,n),er(mnl,4)
    Character             ::Pressure*4,Direction*200,Physical*1
    Integer               ::i,j,Fre_count,Fre_num
    Real                  ::Frequency_variation(Fre_grid,2)
    
    Integer               ::max_mode_abs,min_mode,max_mode
    Real                  ::Corr_max_abs,Del_max_abs,Del_max,Corr_max,Corr_min,Del_min
    
    Integer               ::month_count,month_max=97
    Character             ::Time_C*14
    Real*4                ::Co_position(position),Del_position(position)
    
    Integer               ::latitude_count,longitude_count,position_count
    
    Print*,"*******************************Notification****************************"
    Print*,"  Step A: Extraction,Measurement and Correlation of global field       "
    Print*,"***********************************************************************"
    
    !**************************************************************************************
    !�����ִ���ּ��ȷ��Ƶ���˲����������ֵ
    !**************************************************************************************
    Do Fre_count=1,Fre_grid
        Frequency_variation(Fre_count,1)=Real((Fre_count-1))*(Fre_max-Fre_min)/Real(Fre_grid)+Fre_min
		call bandpass_filter(Month_Signal,Month_Signal,position_num,Coved_ecof,Pressure,Filed,Frequency_variation(Fre_count,1))
		x=Coved_ecof
        call eof(m,n,mnl,x,ks,er,egvt,ecof)
		call Signal_correlation_lag(ecof,Pressure,n,mnl,Corr_max_abs,f_up,Physical,max_mode_abs,Del_max_abs,Direction,Corr_max,max_mode,Del_max,Corr_min,min_mode,Del_min)
		Frequency_variation(Fre_count,2)=Corr_max_abs
    Enddo
    Open(11,file=trim(Direction)//"\Glo_Corr_result\"//trim(Pressure)//"\"//trim(Physical)//"Frequency_Correlation.txt")
    Do Fre_count=1,Fre_grid
        write(11,"(2f18.8)")Frequency_variation(Fre_count,1),Frequency_variation(Fre_count,2)
    Enddo
    close(11)
	
    Corr_max_abs=0.0;
	Do Fre_count=1,Fre_grid
	   if(abs(Frequency_variation(Fre_count,2)).ge.abs(Corr_max_abs))then
	      Corr_max_abs=Frequency_variation(Fre_count,2);
		  f_up=Frequency_variation(Fre_count,1);
	    Endif
	Enddo 	
    !**************************************************************************************
    !�����ִ���ּ�ڶ�ԭʼ���쳣����������ͨ�˲�
    !���������
    !    length_vector_A: ʱ�䳤��
    !           position: �ռ��λ��Ŀ
    !           Pressure: ���������ѹ�߶Ȳ�
    !               f_up: Ƶ���˲������������Ƶ��
    !���������
    !         Coved_ecof: �˲�����쳣��
    !**************************************************************************************
    call bandpass_filter(length_vector_A,length_vector_B,position,Coved_ecof,Pressure,Filed,f_up,Physical)
    Open(200,file=trim(Direction)//"\Glo_Bandfiltered_data\"//trim(Physical)//trim(Pressure)//"_F.txt")
    Do Position_count=1,Position
        write(200,"(<length_vector_A>f18.8)")(Coved_ecof(Position_count,month_count),month_count=1,month_max)
    Enddo
    close(200)
    write(*,*)"Bandfilting of "//trim(Pressure)//" Fini."
    !**************************************************************************************
    !�����ִ���ּ�ڶԴ�ͨ�˲�����쳣��������EOF�ֽ�
    !���������
    !                  m: �ռ�����Ŀ
    !                  n: ʱ������Ŀ
    !                mnl: ���ֽ���������Ŀ
    !                  x: �˲�����쳣��
    !                 ks: EOF�ֽ�ģʽ
    !           Pressure: ���������ѹ�߶Ȳ�
    !���������
    !                 er: EOF�ֽ�����ɷ���ռ����
    !               egvt: ��ſռ�������
    !               ecof: ���ʱ������������
    !**************************************************************************************
    x=Coved_ecof
    call eof(m,n,mnl,x,ks,er,egvt,ecof)

        open(2,file=trim(Direction)//"\Glo_EOF_result\"//trim(Physical)//trim(Pressure)//"er.txt")
        do i=1,mnl
        write(2,"(4f30.8)") (er(i,j),j=1,4)
        enddo
        close(2)
        open(3,file=trim(Direction)//"\Glo_EOF_result\"//trim(Physical)//trim(Pressure)//"ecof.txt")
        do j=1,n
        write(3,"(<mnl>f20.6)") (ecof(i,j),i=1,mnl)
        enddo
        close(3)
        open(4,file=trim(Direction)//"\Glo_EOF_result\"//trim(Physical)//trim(Pressure)//"egvt.txt")
        do i=1,m
        write(4,"(<mnl>f20.6)") (egvt(i,j),j=1,mnl)
        enddo
        close(4)
    write(*,*)"EOF deformation of "//trim(Pressure)//" Fini."
    !**************************************************************************************
    !�����ִ���ּ�ڶ�EOF�ֽ���ʱ������ֵ������Է���
    !���������
    !           Pressure: ���������ѹ�߶Ȳ�
    !               ecof: ʱ������������
    !                  n: ʱ�䳤��
    !                mnl: ���ֽ���������Ŀ
    !           position: �ռ��λ��Ŀ
    !               f_up: Ƶ���˲������������Ƶ��
    !**************************************************************************************
    call Signal_correlation_lag(ecof,Pressure,n,mnl,Corr_max_abs,f_up,Physical,max_mode_abs,Del_max_abs,Direction,Corr_max,max_mode,Del_max,Corr_min,min_mode,Del_min)
    Open(13,file=trim(Direction)//"\Glo_Corr_result\"//trim(Pressure)//"\"//trim(Physical)//"VIP_information.txt")
    write(13,"(3(I5,2f16.6),f16.6)")max_mode_abs,Corr_max_abs,Del_max_abs,max_mode,Corr_max,Del_max,min_mode,Corr_min,Del_min,f_up
    close(13)
    
    Open(15,file=trim(Direction)//"\Glo_Corr_result\"//trim(Pressure)//"\"//trim(Physical)//"Mode_variation.txt")
    Open(14,file=trim(Direction)//"\Time_Month_day.txt")
    Do month_count=1,month_max
        Read(14,"(A14)")Time_C
        write(15,"(A14,f18.8)")Time_C,ecof(max_mode_abs,month_count)
    Enddo
    close(14)
    close(15)
    write(*,*)"Mode correlation with ONI of "//trim(Pressure)//" Fini."
    !**************************************************************************************
    !�����ִ���ּ�ڶ��˲�����쳣��������Է���
    !���������
    !           Pressure: ���������ѹ�߶Ȳ�
    !                  x: �˲�����쳣��
    !                  n: ʱ�䳤��
    !           position: �ռ��λ��Ŀ
    !    
    !**************************************************************************************
    call Position_Correlation(x,position,n,Pressure,Physical,Direction,Co_position,Del_position)
    
    Open(15,file=trim(Direction)//"\Glo_Corr_result\"//trim(Pressure)//"\"//trim(Physical)//"Del_position.txt")
    Open(14,file=trim(Direction)//"\Glo_Corr_result\"//trim(Pressure)//"\"//trim(Physical)//"Corr_position.txt")
    Do latitude_count=1,latitude_max
    Do longitude_count=1,longitude_max
        position_count=(latitude_count-1)*longitude_max+longitude_count
        write(14,"(2I6,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),Co_position(position_count)
        write(15,"(2I6,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),Del_position(position_count)
    Enddo
    Enddo
    close(14)
    close(15)
    write(*,*)"All measurement of "//trim(Pressure)//" Fini."
    
    Print*,"*******************************Notification****************************"
    Print*,"  Step A: Extraction,Measurement and Correlation of global field Fini.      "
    Print*,"***********************************************************************"
    End 



