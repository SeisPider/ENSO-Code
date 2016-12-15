Subroutine Longitude_Height(Direction_in,region)
    Use math_tools
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 本部分代码旨在提取各个经度点、高度点上序列随时间变化的STD，定位出垂直上物理场的EOF分解及!
	!           提取出相关的空间典型场和时间模态                                              !
    ! 输入文件: 各异常场的月均值文件                                                          !
	!                路径在 E:\ONI-Pressure_????\Monthly_Anomaly\                             !                                  !
    ! 输出文件：各区域的STD随经度-高度的分布，EOF分解的典型场                                 !
	!                路径在 C:\Users\Seispider\Documents\Vertical_structure\Longitude-Height  !                                                        !
    ! 相关方法： Pearson相关                                                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                    变量定义区                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer,parameter   ::rk=4,ik=4,Pressure_num=40,Month_length=121,Month_Signal=97,Position_num=2701 &
	                      ,longitude_num=73,latitude_num=37
	Character,intent(in):: Direction_in*200,region*4
	Character*4         :: Pressure_matrix(Pressure_num)
	Character*14        :: Time_series(Month_Signal)
	Integer             :: Pressure_count,Real_Pressure_num=0
	Integer             :: Position_count,Month_count,S_num,T_num
	Integer             :: latitude_min,latitude_max,longitude_min,longitude_max,longitude_count,&
	                       latitude_count,longitude_ALL,latitude_ALL
	Real                :: Vector_standard_in(longitude_num),Vector_standard_out(longitude_num),Vector_height_in(Month_Signal),Vector_height_out(Month_Signal)					   
	Real                :: S_Month_Anoma(Pressure_num,Position_num,Month_Signal),T_Month_Anoma(Pressure_num,Position_num,Month_Signal)
	Real                :: S_longitude_Average(Pressure_num,longitude_num,Month_Signal),&
	                       T_longitude_Average(Pressure_num,longitude_num,Month_Signal),&
						   STD_in(Month_Signal),STD_out, &
						   S_longitude_STD(Pressure_num,longitude_num),T_longitude_STD(Pressure_num,longitude_num), &
						   S_Time_height_average(Pressure_num,Month_Signal),T_Time_height_average(Pressure_num,Month_Signal)
	Real                :: S_sum,T_sum,P,P0=1013.15!unit=hpa standard atmospheric pressure
	Logical             :: alive
	Character*200       :: Direction_out,S_File_Direction,T_File_Direction,ONI_direction   !设置输入,输出文件路径
	Data Pressure_matrix /"30","50","70","100","125",&
                       	"150","175","200","225","250",&
						"275","300","325","350","375",&
						"400","425","450","475","500",&
						"525","550","575","600","625",&
						"650","675","700","725","750",&
						"775","800","825","850","875",&
						"900","925","950","975","1000"/
	Direction_out  ="C:\Users\Seispider\Documents\Vertical_structure\Longitude-Height\Longitude_Height_STD_Var\"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                         读入文件                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Real_Pressure_num=0
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
	!                                           求标准差                                     !
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
    !此处可添加任意多选项以形成对全球任意区域做variation和STD                                !
    !   region是用户自行命名的区域名称。                                                     !
    !   longitude_min:表示该区域左侧的经度在序列中的下标数[-180 => 1;-175 => 2;-170 => 3;..] !
    !   longitude_max:表示该区域右侧的经度在序列中的下标数[-180 => 1;-175 => 2;-170 => 3;..] !
    !   latitude_min :表示该区域下侧的纬度在序列中的下标数[-90 => 1;-85 => 2;-80 => 3;..]    !
    !   latitude_max :表示该区域上侧的纬度在序列中的下标数[-90 => 1;-85 => 2;-80 => 3;..]    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!          做经度上叠加消除纬度维,时间维上求STD,对同一高度层上水汽STD做标准化            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do Pressure_count=1,Real_Pressure_num
	Do longitude_count=longitude_min,longitude_max
	   Do Month_count=1,Month_Signal
	   S_sum=0;T_sum=0;S_num=0;T_num=0;
	   Do latitude_count=latitude_min,latitude_max
	      Position_count=(latitude_count-1)*longitude_num+longitude_count
		  S_sum=S_sum+S_Month_Anoma(Pressure_count,Position_count,Month_count)
		  T_sum=T_sum+T_Month_Anoma(Pressure_count,Position_count,Month_count)
		  S_num=S_num+1;T_num=T_num+1;
	   Enddo
	   S_longitude_Average(Pressure_count,longitude_count,Month_count)=S_sum/S_num
	   T_longitude_Average(Pressure_count,longitude_count,Month_count)=T_sum/T_num
	   Enddo
	   STD_in(:)=S_longitude_Average(Pressure_count,longitude_count,:)
	   call Standard_deviation(STD_in,STD_out,Month_Signal)
	   S_longitude_STD(Pressure_count,longitude_count)=STD_out
	   
	   STD_in(:)=T_longitude_Average(Pressure_count,longitude_count,:)
	   call Standard_deviation(STD_in,STD_out,Month_Signal)
	   T_longitude_STD(Pressure_count,longitude_count)=STD_out
	
	Enddo
	Vector_standard_in(:)=S_longitude_STD(Pressure_count,:)
	call vector_standard(Vector_standard_in,Vector_standard_out,longitude_num)
	S_longitude_STD(Pressure_count,:)=Vector_standard_out(:)
	Enddo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!          做经纬度上叠加消除经纬度，对各个高度层上同一时间值序列做标准化                !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do Pressure_count=1,Real_Pressure_num
	Do Month_count=1,Month_Signal
	   S_sum=0;T_sum=0;S_num=0;T_num=0;
	   Do latitude_count=latitude_min,latitude_max
	   Do longitude_count=longitude_min,longitude_max
	   Position_count=(latitude_count-1)*longitude_num+longitude_count
	   S_sum=S_sum+S_Month_Anoma(Pressure_count,Position_count,Month_count)
	   T_sum=T_sum+T_Month_Anoma(Pressure_count,Position_count,Month_count)
	   S_num=S_num+1;T_num=T_num+1;
	   Enddo
	   Enddo
	   S_Time_height_average(Pressure_count,Month_count)=S_sum/S_num;
	   T_Time_height_average(Pressure_count,Month_count)=T_sum/T_num;
	Enddo
	Vector_height_in(:)=S_Time_height_average(Pressure_count,:)
	call vector_standard(Vector_height_in,Vector_height_out,Month_Signal)
	S_Time_height_average(Pressure_count,:)=Vector_height_out(:)

	Enddo  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                          输出STD结果                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(13,file=trim(Direction_out)//trim(region)//"_S_ln_L2H_N_STD.txt")
	Open(14,file=trim(Direction_out)//trim(region)//"_S_ch_L2H_N_STD.txt")
	Open(15,file=trim(Direction_out)//trim(region)//"_T_ln_L2H_STD.txt")
	Open(16,file=trim(Direction_out)//trim(region)//"_T_ch_L2H_STD.txt")
	Do Pressure_count=1,Real_Pressure_num
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!将气压值转换为其对数值：H~-ln(p/p0)其中p0为标准大气压
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Read(Pressure_matrix(Pressure_count),*)P
	P=-log(P/P0)
   	Do longitude_count=longitude_min,longitude_max
	   longitude_ALL=(longitude_count-1)*5
	   if(longitude_ALL.lt.180)then
	      longitude_ALL=longitude_ALL+180
	   else
	      longitude_ALL=longitude_ALL-180
	   Endif
	   write(13,"(I10,f10.6,f18.8)")longitude_ALL,P,S_longitude_STD(Pressure_count,longitude_count)
	   write(14,"(I10,A7,f18.8)")longitude_ALL,Pressure_matrix(Pressure_count),S_longitude_STD(Pressure_count,longitude_count)
	   write(15,"(I10,f10.6,f18.8)")longitude_ALL,P,T_longitude_STD(Pressure_count,longitude_count)
	   write(16,"(I10,A7,f18.8)")longitude_ALL,Pressure_matrix(Pressure_count),T_longitude_STD(Pressure_count,longitude_count)
	Enddo
	Enddo
    close(13);close(14);close(15);close(16);
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                       输出时间变化结果                                 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(13,file=trim(Direction_out)//trim(region)//"_S_ln_L2H_N_Var.txt")
	Open(14,file=trim(Direction_out)//trim(region)//"_S_ch_L2H_N_Var.txt")
	Open(15,file=trim(Direction_out)//trim(region)//"_T_ln_L2H_Var.txt")
	Open(16,file=trim(Direction_out)//trim(region)//"_T_ch_L2H_Var.txt")
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                          载入时间序列                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(17,file=trim(Direction_out)//"Time_Month_day.txt")
	Do Month_count=1,Month_Signal
	   Read(17,"(A14)")Time_series(Month_count)
	Enddo
	close(17)
	Do Pressure_count=1,Real_Pressure_num
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!将气压值转换为其对数值：H~-ln(p/p0)其中p0为标准大气压
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Read(Pressure_matrix(Pressure_count),*)P
	P=-log(P/P0)
	Do Month_count=1,Month_Signal
	   write(13,"(A14,f10.6,f18.8)")Time_series(Month_count),P,S_Time_height_average(Pressure_count,Month_count)
	   write(14,"(A14,A7,f18.8)")Time_series(Month_count),Pressure_matrix(Pressure_count),S_Time_height_average(Pressure_count,Month_count)
	   write(15,"(A14,f10.6,f18.8)")Time_series(Month_count),P,T_Time_height_average(Pressure_count,Month_count)
	   write(16,"(A14,A7,f18.8)")Time_series(Month_count),Pressure_matrix(Pressure_count),T_Time_height_average(Pressure_count,Month_count)
	Enddo
	Enddo
	close(13);close(14);close(15);close(16);
	
	
	End Subroutine