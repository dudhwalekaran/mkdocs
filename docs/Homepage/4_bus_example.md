``` data title="busno_fdlf.dat"
1        %Slack bus number
0.001    %Loadflow convergence tolerance
4        %Number of buses in the system
4        %Number of lines
0        %Number of transformers
1        %Number of PV buses
1        %For accounting Q-limits set this bit 1, otherwise 0  
0        %For accounting Load modelling set this bit 1, otherwise 0  
0        %For accounting Frequency Dependency Load modelling set this bit 1, otherwise 0  
4        %Number of load buses
0        %Number of shunts
1.0      %Slack bus voltage magnitude   
0        %NO of HVDC links
```


``` data title="busno_nr.dat"
1        %Slack bus number
0.001    %Loadflow convergence tolerance
4        %Number of buses in the system
4        %Number of lines
0        %Number of transformers
1        %Number of PV buses
1        %For accounting Q-limits set this bit 1, otherwise 0  
4        %Number of load buses
0        %Number of shunts
1.0      %Slack bus voltage magnitude   
```


``` data title="lfl_NR_sc0.dat"
  1 	     1.000000 	    0.000000	    1.884566	    2.667856	 0.500000	 0.309900 
  2 	     0.941383 	   -0.533208	    0.000000	    0.000000	 1.700000	 1.053500 
  3 	     0.942312 	   -1.631389	    0.000000	    0.000000	 2.000000	 1.239400 
  4 	     0.951521 	    2.636509	    3.180000	    0.400001	 0.800000	 0.495800 
```

``` data title="nt.dat"
1 2 0.01008 0.05040 0.1025
1 3 0.00744 0.03720 0.0775
2 4 0.00744 0.03720 0.0775
3 4 0.01272 0.06360 0.1275
```


``` data title="pvpq.dat"
      4   1.02     3.18
      1   0.5      0.3099
      4   0.8      0.4958
      2   1.7      1.0535
      3   2        1.2394  
```

``` data title="Qlim.dat"
4 0 0.4
```


``` data title="report_NR_sc0.dat"
                   Detailed Report of Load flow (NR method)
                   ----------------------------
Iter. No = 1 Max. Real power mismatch at bus = 4, Max. mismatch = 0.22860224
Iter. No = 2 Max. Real power mismatch at bus = 4, Max. mismatch = 0.00430169
__________________________________________________________________________
Converged loadflow iteration No = 3, Max. mismatch = 0.00000147
----------------------------------
 Q-limits at PV-buses accounted
 ------------------------------------------------------------
 Loadflow results:
---------------------------------------------------------------------------------------

Bus No     VbO          thetaO           PGO             QGO         PLO         QLO 
---------------------------------------------------------------------------------------
  1 	     1.000000 	    0.000000	    1.884566	    2.667856	 0.500000	 0.309900 
  2 	     0.941383 	   -0.533208	    0.000000	    0.000000	 1.700000	 1.053500 
  3 	     0.942312 	   -1.631389	    0.000000	    0.000000	 2.000000	 1.239400 
  4 	     0.951521 	    2.636509	    3.180000	    0.400001	 0.800000	 0.495800 
---------------------------------------------------------------------------------------
Line flows:
----------------------------------------------------------------------------------
                       Line flows                                Line flows 
                   _____________________                   _______________________
 From   To         P-flow         Q-flow    From   To      P-flow           Q-flow
----------------------------------------------------------------------------------
  1 	   2	       0.3910	       1.0344	   2	    1	      -0.3775	      -1.0640
  1 	   3	       0.9936	       1.3236	   3	    1	      -0.9725	      -1.2910
  2 	   4	      -1.3225	       0.0105	   4	    2	       1.3372	      -0.0064
  3 	   4	      -1.0275	       0.0516	   4	    3	       1.0428	      -0.0894
----------------------------------------------------------------------------------
Total real power losses in the system =  0.064567
Total reactive power losses in the system = -0.030743
```



``` data title="trans.tbl"
F Qlim_data.dat                                                                                                                                                                                                   	Qlim_data.dat
F busno_fdlf.dat                                                                                                                                                                                                  	busno_fdlf.dat
F busno_nr.dat                                                                                                                                                                                                    	busno_nr.dat
F lfl_NR_sc0.dat                                                                                                                                                                                                  	lfl_NR_sc0.dat
F nt.dat                                                                                                                                                                                                          	nt.dat
F pvpq.dat                                                                                                                                                                                                        	pvpq.dat
F report_NR_sc0.dat                                                                                                                                                                                               	report_NR_sc0.dat
```