Program Vertical_structure
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           �����ִ���ּ��������ģ��ĵ��ù��������������ݴ�����              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Character::Direction_in*200,region*4                                             !
    !Direction_in="E:\ONI-Pressure1\"                                                 !
    !region="Glo"                                                                     !
	!call Normalize(Direction_in,region)                                              !
    !region="Mid"                                                                     !
	!call Normalize(Direction_in,region)                                              !
    !region="Equ"                                                                     ! 
	!call Normalize(Direction_in,region)                                              !
    !region="Five"                                                                    !
	!call Normalize(Direction_in,region)                                              !
    !region="Nino"                                                                    !
	!call Normalize(Direction_in,region)                                              !
    !Direction_in="E:\ONI-Pressure2\"                                                 !
    !region="Glo"                                                                     !
	!call Normalize(Direction_in,region)                                              !
    !region="Mid"                                                                     !
	!call Normalize(Direction_in,region)                                              !
    !region="Equ"                                                                     !
	!call Normalize(Direction_in,region)                                              !
    !region="Five"                                                                    !
	!call Normalize(Direction_in,region)                                              !
    !region="Nino"                                                                    ! 
	!call Normalize(Direction_in,region)                                              !
	!����ʹ��ʱ��:2016-12-09                                                          !
	!����Ԥ�ڹ���:����˲��Ͳ��˲��źŵĹ�һ����ʱ��ƫ��                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!             �����ģ��ּ�ڵ���ģ�鸴����س����ӳٳ���Ƶ��-���ϵ����ͼ         !
	!                           �����ӳ����ϵ����ͼ����ȡEOF���ͳ�                   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character::Direction_in*200,region*4,Filt_Par*10,Direction_out*200
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !region="Glo";Filt_Par="Filter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\";  !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Mid";Filt_Par="Filter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\";  !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
    !region="Equ";Filt_Par="Filter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\";  !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Five";Filt_Par="Filter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\"; !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Nino";Filt_Par="Filter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\"; !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Glo";Filt_Par="UFilter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\"; !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Mid";Filt_Par="UFilter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\"; !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
    !region="Equ";Filt_Par="UFilter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\"; !
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Five";Filt_Par="UFilter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\";!
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
	!region="Nino";Filt_Par="UFilter";Direction_in="E:\ONI-Pressure_"//trim(Filt_Par)//"\";!
	!call Corr_Field(Direction_in,region,Filt_Par)                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Direction_in="E:\ONI-Pressure_Filter\Monthly_Anomaly\";region="Glo";                  !
    !call Longitude_Height(Direction_in,region)                                            !
    !Direction_in="E:\ONI-Pressure_Filter\Monthly_Anomaly\";region="Mid";                  !
    !call Longitude_Height(Direction_in,region)                                            !
    !Direction_in="E:\ONI-Pressure_Filter\Monthly_Anomaly\";region="Equ";                  !
    !call Longitude_Height(Direction_in,region)                                            !
    !Direction_in="E:\ONI-Pressure_Filter\Monthly_Anomaly\";region="Five";                 !
    !call Longitude_Height(Direction_in,region)                                            !
    !Direction_in="E:\ONI-Pressure_Filter\Monthly_Anomaly\";region="Nino";                 !
    !call Longitude_Height(Direction_in,region)                                            !
    !�����ִ���ּ�ڵ����ļ���ȡ��ͬ����ľ�ֵ��ʱ��仯����Ϣ�������ȡ��߶ȵ���������ʱ���!
    !          ����STD��2016-12-10                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Direction_out="C:\Users\Seispider\Documents\Vertical_structure\Signal_correlation\"   !
    !Direction_in="E:\ONI-Pressure_Filter\"                                                !
    !call Signal_Corr(Direction_in,Direction_out)                                          !
	!�����ִ���ּ�����ø���ȡ�������ź����˲���ĵ��ͳ���Ƥ��ѷ��أ��õ���س��ֲ����ӳٳ�!
	!          ���ֲ����˲�����ѡ������˲�Ƶ�ʡ�: 2016-12-13
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    End Program
