Subroutine Normalize(Direction_in,region)
    Use math_tools
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           本部分代码旨在对数据做归一化及数据平移工作                                !
    ! 输入文件:  原始EOF分解得到的最大模态的文件。路径在 E:\ONI-Pressure?\*_Corr_Result\  !
    ! 输出文件： 时间模态平移前，平移后的结果。路径在C:\Users\Seispider\Documents\Vertical!
    !           _structure\Normalize                                                      !
    ! 归一方法： x'= (x-mu)/delta ! mu是样本均值，delta是样本一倍标准差                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  变量定义区                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer,parameter::rk=4,ik=4,Pressure_num=40,Month_length=121,Month_Signal=97
	Character,intent(in)::Direction_in*200,region*4
	Character*200 :: Direction_out   !设置输入,输出文件路径
	Character*4   :: Pressure_matrix(Pressure_num)
	Character*14  :: Time_Month_day_Ext(Month_length),Raw_Month_day(Month_Signal)
	Integer       :: Pressure_count,Month_count,S_mode,T_mode
	Real          :: S_corr,T_corr,S_del,T_del
	Real          :: S_Signal(Month_Signal),T_Signal(Month_Signal),ONI_Month_length(Month_length)
	Real          :: Vector_Signal_in(Month_Signal),Vector_Signal_out(Month_Signal),Vector_ONI_in(Month_length),Vector_ONI_out(Month_length)
    Integer       :: Time_min,Time_max
    Logical       :: alive
	Data Pressure_matrix /"30","50","70","100","125",&
                       	"150","175","200","225","250",&
						"275","300","325","350","375",&
						"400","425","450","475","500",&
						"525","550","575","600","625",&
						"650","675","700","725","750",&
						"775","800","825","850","875",&
						"900","925","950","975","1000"/
    Direction_out  ="C:\Users\Seispider\Documents\Vertical_structure\Normalize\filter\"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                           载入原始的时间坐标文件并输出ONI文件                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(10,file=trim(Direction_out)//"Time_Month_day_Ext.txt")
	Open(9,file=trim(Direction_out)//"ONI_index.txt")
	Open(8,file=trim(Direction_out)//"ONI_Month_length.txt")
	Do  Month_count=1,Month_length
	    Read(10,"(A14)")Time_Month_day_Ext(Month_count)
		Read(9,*)ONI_Month_length(Month_count)
	Enddo
	close(9)
	close(10)
	Vector_ONI_in=ONI_Month_length(:)
	call vector_standard(Vector_ONI_in,Vector_ONI_out,Month_length)
	ONI_Month_length(:)=Vector_ONI_out(:)
	Do Month_count=1,Month_length
	    write(8,"(A14,f18.8)")Time_Month_day_Ext(Month_count),ONI_Month_length(Month_count)
    Enddo
    close(8)	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  对未滤波的数据做处理                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Direction_in   ="E:\ONI-Pressure2\"
	Do Pressure_count=1,Pressure_num
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                            对各层数据做归一化和偏移处理                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  读入文件                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!region = "Glo"
    inquire(file=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\SMode_variation.txt",exist=alive)
    if(.not.alive)then
        cycle
    else
	Open(11,file=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\SMode_variation.txt")
	Open(12,file=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\SVIP_information.txt")
	Open(13,file=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\TMode_variation.txt")
	Open(14,file=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\TVIP_information.txt")
	Do Month_count=1,Month_Signal
	   Read(11,"(A14,f18.8)")Raw_Month_day(Month_count),S_Signal(Month_count)
	   Read(13,"(A14,f18.8)")Raw_Month_day(Month_count),T_Signal(Month_count)
	Enddo
	Read(12,"(I5,2f16.6)")S_mode,S_corr,S_del
	Read(14,"(I5,2f16.6)")T_mode,T_corr,T_del
	close(11)
	close(12)
	close(13)
	close(14)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  做数据归一化                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	Vector_Signal_in(:)=S_Signal(:)
	call vector_standard(Vector_Signal_in,Vector_Signal_out,Month_Signal)
	S_Signal(:) = Vector_Signal_out(:) 
	
	Vector_Signal_in(:)=T_Signal(:)
	call vector_standard(Vector_Signal_in,Vector_Signal_out,Month_Signal)
	T_Signal(:) = Vector_Signal_out(:)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                  做数据输出                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Open(15,file=trim(Direction_out)//trim(region)//"_F_"//trim(Pressure_matrix(Pressure_count))//"_S_UM.txt")
	Open(16,file=trim(Direction_out)//trim(region)//"_F_"//trim(Pressure_matrix(Pressure_count))//"_S_M.txt")
	Open(17,file=trim(Direction_out)//trim(region)//"_F_"//trim(Pressure_matrix(Pressure_count))//"_T_UM.txt")
	Open(18,file=trim(Direction_out)//trim(region)//"_F_"//trim(Pressure_matrix(Pressure_count))//"_T_M.txt")
    Do Month_count=1,Month_Signal
       write(15,"(A14,f18.8)")Raw_Month_day(Month_count),S_Signal(Month_count)
       write(16,"(A14,f18.8)")Time_Month_day_Ext(Month_count+12+int(S_del)),S_Signal(Month_count)
       write(17,"(A14,f18.8)")Raw_Month_day(Month_count),T_Signal(Month_count)
       write(18,"(A14,f18.8)")Time_Month_day_Ext(Month_count+12+int(T_del)),T_Signal(Month_count)
    Enddo	   
    Endif
    Enddo
    
   

End Subroutine
