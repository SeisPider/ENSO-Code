SUBROUTINE bandpass_filter(length_vector_A,length_vector_B,position_num,Coved_ecof,Filed,f_up)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: rk=4, ik=4
    REAL,PARAMETER          :: Pi=4*atan(1.0)
    REAL                    :: f_up
    INTEGER, INTENT(IN)     :: length_vector_A,length_vector_B, position_num
    REAL                    :: vector_A(length_vector_A,2), vector_B(length_vector_B,2)
    REAL, INTENT(OUT)       :: Coved_ecof(position_num,length_vector_A)
    REAL, INTENT(IN)        :: Filed(position_num,length_vector_A)
    INTEGER                 :: time_count, i, position_count
    REAL                    :: Convo_ed_vector((length_vector_A+length_vector_B-1),2)

 
    DO time_count=1,length_vector_B
        vector_B(time_count,1) = time_count
        vector_B(time_count,2) = sin(2*Pi*f_up*time_count)/(Pi*time_count)
    ENDDO
    

    DO position_count=1,position_num
        DO time_count=1,length_vector_A
            vector_A(time_count,1) = time_count
            vector_A(time_count,2) = Filed(position_count,time_count)
        ENDDO
        CALL Convolution(vector_A,length_vector_A,vector_B,length_vector_B,Convo_ed_vector)
        Coved_ecof(position_count,1:length_vector_A) = Convo_ed_vector(1:length_vector_A,2)
    ENDDO  
END SUBROUTINE
    
SUBROUTINE Convolution(vector_A,length_vector_A,vector_B,length_vector_B,Convo_ed_vector)
    IMPLICIT NONE
    !Integer,parameter     ::rk=selected_real_kind(15,307),ik=selected_int_kind(9)
    INTEGER, PARAMETER      :: rk = 4, ik = 4
    INTEGER(ik), INTENT(IN) :: length_vector_A, length_vector_B
    REAL(rk), INTENT(IN)    :: vector_A(length_vector_A,2), vector_B(length_vector_B,2)
    REAL(rk), INTENT(OUT)   :: Convo_ed_vector((length_vector_A+length_vector_B-1),2)
    REAL(rk)                :: sum, Element_A, Element_B
	INTEGER(ik)             :: Total_count
    INTEGER(ik)             :: Time_count, vector_Ak, vector_BNK, Argument_max_A, Argument_min_A, Argument_max_B, Argument_min_B
    INTEGER(ik)             :: Argument_max, Argument_min, Argument
    Argument_max_A=vector_A(1,1); Argument_min_A=vector_A(1,1);
    DO Time_count = 1, length_vector_A
	    IF(vector_A(Time_count,1).ge.Argument_max_A) THEN
	        Argument_max_A = vector_A(Time_count,1)
    	ENDIF
	    IF(vector_A(Time_count,1).le.Argument_min_A) THEN
	        Argument_min_A = vector_A(Time_count,1)
	    ENDIF
	ENDDO
	
    Argument_max_B=vector_B(1,1);Argument_min_B=vector_B(1,1);
    DO Time_count=1,length_vector_B
	    IF(vector_B(Time_count,1).ge.Argument_max_B) THEN
	        Argument_max_B = vector_B(Time_count,1)
	    ENDIF
	    IF(vector_B(Time_count,1).le.Argument_min_B) THEN
	        Argument_min_B = vector_B(Time_count,1)
	    ENDIF
	ENDDO
    Argument_max = Argument_max_A+Argument_max_B
	Argument_min = Argument_min_A+Argument_min_B
	
	Time_count   = 1
	DO Argument = Argument_min,Argument_max
	    Convo_ed_vector(Time_count,1) = Argument
        Time_count                    = Time_count + 1
	ENDDO
	
	Total_count=1;
	DO Argument=Argument_min,Argument_max
	   sum=0;
	   DO vector_Ak = Argument_min_A,Argument_max_A
	      
	      vector_BNK = Argument-vector_Ak 
		  
		  DO Time_count=1,length_vector_A
		      IF(vector_Ak.eq.vector_A(Time_count,1)) THEN
		          Element_A = vector_A(Time_count,2)
                  CYCLE
              ENDIF
		  ENDDO
		  
		  DO Time_count=1,length_vector_B
		      IF(vector_BNk.eq.vector_B(Time_count,1)) THEN
		          Element_B = vector_B(Time_count,2)
                  CYCLE
              ENDIF
		  ENDDO
		  sum = sum + Element_A*Element_B
	      Element_A=0; Element_B=0; 
       ENDDO
       Convo_ed_vector(Total_count,2) = sum
       Total_count                    = Total_count + 1
    ENDDO
END SUBROUTINE Convolution