Subroutine Corr_Field(Direction_in,region,Filt_Par)
    Use math_tools
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! �����ִ���ּ�ڶ���س����ӳ�ʱ�䳡���ӳ�-���ϵ����ͼ���˲�Ƶ��-���ϵ����ͼ��EOF�ռ䳡!
    ! �����ļ�:  ��س��ļ����ӳ�ʱ�䳡�ļ����˲�Ƶ��-���ϵ����ͼ���ź�����                 !
	!                ·���� E:\ONI-Pressure?\*_Corr_Result\                                  !
	!            EOF�ռ䳡                                                                   ��
	!                ·���� E:\ONI-Pressure?\*_EOF_Result\                                   !
    ! ����ļ��� ��س����ӳ�ʱ�䳡���ӳ�-���ϵ����ͼ���˲�Ƶ��-���ϵ����ͼ��EOF�ռ䳡     !
	!                ·���� C:\Users\Seispider\Documents\Vertical_structure\Cor_Field        !                                                        !
    ! ��ط����� Pearson���                                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                    ����������                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer,parameter::rk=4,ik=4,Pressure_num=40,Month_length=121,Month_Signal=97,t_num=25
	Character,intent(in)::Direction_in*200,region*4,Filt_Par*10
	Character*200 :: Direction_out,File_Direction,ONI_direction   !��������,����ļ�·��
	Character*4   :: Pressure_matrix(Pressure_num)
	Character*400 :: Copy_command
	Character     :: Raw_Time_Month*14
	Real          :: S_egvt(2701,Month_Signal),T_egvt(2701,Month_Signal),S_Signal(Month_Signal),T_Signal(Month_Signal)
	Real          :: Co_t(25),Del_t(25)
	Integer       :: S_mode,T_mode,Position_count,Position_num,mode_count,Pressure_count,t_count,Month_count,n
	Integer       :: latitude_min,latitude_max,longitude_min,longitude_max,longitude_count,latitude_count,longitude_ALL,latitude_ALL
	Real          :: Del_Revert
	Logical       :: alive
	Data Pressure_matrix /"30","50","70","100","125",&
                       	"150","175","200","225","250",&
						"275","300","325","350","375",&
						"400","425","450","475","500",&
						"525","550","575","600","625",&
						"650","675","700","725","750",&
						"775","800","825","850","875",&
						"900","925","950","975","1000"/
    Direction_out  ="C:\Users\Seispider\Documents\Vertical_structure\Cor_Field\"//trim(Filt_Par)//"\"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                    �����ļ���                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do Pressure_count=1,Pressure_num
        
	File_Direction=trim(Direction_in)//trim(region)//"_Corr_Result\"//trim(Pressure_matrix(Pressure_count))//"\"
	inquire(file=trim(File_Direction)//"SCorr_position.txt",exist=alive)
	if(.not.alive)then
        cycle
    else
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                    �����ļ���Ŀ���ļ��У���س����ӳٳ���Ƶ��-���ϵ����ͼ             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    Copy_command = "copy  "//trim(File_Direction)//"SCorr_position.txt  "//trim(Direction_out)//&
		                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Corr.txt"
						call system(Copy_command)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!    Copy_command = "copy  "//trim(File_Direction)//"SDel_position.txt  "//trim(Direction_out)//&                   !
	!	                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Del.txt" !
	!					call system(Copy_command)                                                                       !
    !2016-12-10��ԭʼ�ĸ�ֵ��ʾ�ӳ٣���ֵ��ʾ��ǰ                                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    Copy_command = "copy  "//trim(File_Direction)//"SFrequency_Correlation.txt  "//trim(Direction_out)//&
		                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Fre.txt"	
                        call system(Copy_command)						
	    Copy_command = "copy  "//trim(File_Direction)//"TCorr_position.txt  "//trim(Direction_out)//&
		                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Corr.txt"
						call system(Copy_command)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !    Copy_command = "copy  "//trim(File_Direction)//"TDel_position.txt  "//trim(Direction_out)//&                   !
	!	                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Del.txt" !
	!					call system(Copy_command)                                                                       !
    !2016-12-10��ԭʼ�ĸ�ֵ��ʾ�ӳ٣���ֵ��ʾ��ǰ                                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    Copy_command = "copy  "//trim(File_Direction)//"TFrequency_Correlation.txt  "//trim(Direction_out)//&
		                trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Fre.txt"
                        call system(Copy_command)						

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                   ��ȡ����Ӧ�Ŀռ䳡                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
        Open(11,file=trim(File_Direction)//"SVIP_information.txt")	
		Open(12,file=trim(File_Direction)//"TVIP_information.txt")
        Read(11,"(I5)")S_mode
	    Read(12,"(I5)")T_mode
        Open(13,file=trim(Direction_in)//trim(region)//"_EOF_Result\"//"S"//trim(Pressure_matrix(Pressure_count))//"egvt.txt")
        Open(14,file=trim(Direction_in)//trim(region)//"_EOF_Result\"//"T"//trim(Pressure_matrix(Pressure_count))//"egvt.txt")
        if(trim(region).eq."Glo")then
		     Position_num=2701
	    elseif(trim(region).eq."Mid")then
		     Position_num=1825
		elseif(trim(region).eq."Equ")then
		     Position_num=949
		elseif(trim(region).eq."Five")then
		     Position_num=219
	    elseif(trim(region).eq."Nino")then
		     Position_num=33
	    Endif	 
		Do Position_count=1,Position_num		
		   Read(13,"(<S_mode>f20.6)")(S_egvt(Position_count,mode_count),mode_count=1,S_mode)
		   Read(14,"(<T_mode>f20.6)")(T_egvt(Position_count,mode_count),mode_count=1,T_mode)
		Enddo
		close(11);close(12);close(13);close(14);
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                               �����ȡ����Ӧ�Ŀռ䳡                                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Open(15,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Spat.txt")
		Open(16,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Spat.txt")
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
	    Endif
        Position_count=0
        Do latitude_count=latitude_min,latitude_max
            Do longitude_count=longitude_min,longitude_max
            Position_count=Position_count+1
			write(15,"(2I6,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),S_egvt(Position_count,S_mode)
			write(16,"(2I6,f18.8)")((longitude_count-1)*5-180),((latitude_count-1)*5-90),T_egvt(Position_count,S_mode)
		    Enddo
		Enddo
		close(15);close(16);
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                    ����ENSO�źŲ���ONI����صõ��ӳ�-���ϵ����ͼ                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Open(17,file=trim(File_Direction)//"SMode_variation.txt")
		Open(18,file=trim(File_Direction)//"TMode_variation.txt")
        Do Month_count=1,Month_Signal
            Read(17,"(A14,f18.8)")Raw_Time_Month,S_Signal(Month_count)
			Read(18,"(A14,f18.8)")Raw_Time_Month,T_Signal(Month_count)
		Enddo
		close(17)
		close(18)
		n = Month_Signal
        ONI_direction="C:\Users\Seispider\Documents\Vertical_structure\Cor_Field"
		call Signal_Correlation(S_Signal,n,ONI_direction,Co_t,Del_t)
		Open(19,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Del_Corr.txt")
        Do 	t_count=1,t_num
            write(19,"(2f12.6)")(-Del_t(t_count)),Co_t(t_count)
        Enddo
        close(19)
		call Signal_Correlation(T_Signal,n,ONI_direction,Co_t,Del_t)
		Open(20,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Del_Corr.txt")
        Do 	t_count=1,t_num
            write(20,"(2f12.6)")(-Del_t(t_count)),Co_t(t_count)
        Enddo
        close(20)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                       ���ӳٳ���ʱ�䷴��                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Open(21,file=trim(File_Direction)//"SDel_position.txt")
	Open(22,file=trim(File_Direction)//"TDel_position.txt")
	Open(23,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"S_Del.txt")
	Open(24,file=trim(Direction_out)//trim(region)//"_"//trim(Pressure_matrix(Pressure_count))//"_"//trim(Filt_Par)//"_"//"T_Del.txt")
	Do Position_count=1,Position_num
	    Read(21,"(2I6,f18.8)")longitude_ALL,latitude_ALL,Del_Revert
		Del_Revert=(-Del_Revert)
		write(23,"(2I6,f18.8)")longitude_ALL,latitude_ALL,Del_Revert
	    Read(22,"(2I6,f18.8)")longitude_ALL,latitude_ALL,Del_Revert
		Del_Revert=(-Del_Revert)
		write(24,"(2I6,f18.8)")longitude_ALL,latitude_ALL,Del_Revert
	Enddo
    close(21);close(22);close(23);close(24);	
	Endif
	Enddo
	
End Subroutine