Program Vertical_structure
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           本部分代码旨在做各分模块的调用工作并完成相关数据处理工作              !
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
	!代码使用时间:2016-12-09                                                          !
	!代码预期功能:完成滤波和不滤波信号的归一化和时间偏移                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!             下面的模块旨在调用模块复制相关场，延迟场，频率-相关系数线图         !
	!                           计算延迟相关系数线图并提取EOF典型场                   !
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
    !本部分代码旨在调用文件提取不同区域的均值随时间变化的信息及各经度、高度点上数据随时间变!
    !          话的STD：2016-12-10                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Direction_out="C:\Users\Seispider\Documents\Vertical_structure\Signal_correlation\"   !
    !Direction_in="E:\ONI-Pressure_Filter\"                                                !
    !call Signal_Corr(Direction_in,Direction_out)                                          !
	!本部分代码旨在利用各提取出来的信号与滤波后的典型场做皮尔逊相关，得到相关场分布，延迟场!
	!          场分布。滤波参数选择最佳滤波频率。: 2016-12-13
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    End Program
