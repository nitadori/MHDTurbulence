	  //direction dependent w,u,v
	Const[mudn] = Plefte[nden]; // rho
	Const[muvw] = Plefte[nve1]*Plefte[nden]; // rho v_x
	Const[muvu] = Plefte[nve2]*Plefte[nden]; // rho v_y
	Const[muvv] = Plefte[nve3]*Plefte[nden]; // rho v_z
        Const[muet] = Plefte[nene]*Plefte[nden]  // e_i
	              +0.5e0*Plefte[nden]*(                  
                          +Plefte[nve1]*Plefte[nve1]                 
                          +Plefte[nve2]*Plefte[nve2]                 
                          +Plefte[nve3]*Plefte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Plefte[nbm1]*Plefte[nbm1]                 
                          +Plefte[nbm2]*Plefte[nbm2]                 
                          +Plefte[nbm3]*Plefte[nbm3]); // + B^2/2

	Const[mubw] = Plefte[nbm1]; // b_x
	Const[mubu] = Plefte[nbm2]; // b_y
	Const[mubv] = Plefte[nbm3]; // b_z
	Const[mubp] = Plefte[nbps]; // psi
	for(int n=0; n<ncomp; n++){
	  Const[must+n] = Plefte[nden]*Plefte[nst+n]; // composition
	}
	
        double  ptl = Plefte[npre] + ( Plefte[nbm1]*Plefte[nbm1]
                                      +Plefte[nbm2]*Plefte[nbm2]
				      +Plefte[nbm3]*Plefte[nbm3])/2.0e0;
	//direction dependent, nve2 or nbm2
	Const[mfdn] = Plefte[nden]*Plefte[nve2];
	Const[mfvw] = Plefte[nden]*Plefte[nve1]*Plefte[nve2] 
	                           -Plefte[nbm1]*Plefte[nbm2];
	Const[mfvu] = Plefte[nden]*Plefte[nve2]*Plefte[nve2]
      	                     + ptl -Plefte[nbm2]*Plefte[nbm2];// p diagnonal
        Const[mfvv] = Plefte[nden]*Plefte[nve3]*Plefte[nve2]
	                          - Plefte[nbm3]*Plefte[nbm2];
        Const[mfet] = (Const[muet]+ptl)*Plefte[nve2]
                           -( Plefte[nbm1]*Plefte[nve1]
                             +Plefte[nbm2]*Plefte[nve2]
			     +Plefte[nbm3]*Plefte[nve3])*Plefte[nbm2];

	// direction dependent 3, 2, 1, 2
	Const[mfbu] =  0.0;
	Const[mfbv] =  Plefte[nbm3]*Plefte[nve2]
	              - Plefte[nve3]*Plefte[nbm2];
	Const[mfbw] =  Plefte[nbm1]*Plefte[nve2]
	              - Plefte[nve1]*Plefte[nbm2];
	Const[mfbp] = 0.0e0;  // psi
	
	for(int n=0; n<ncomp; n++){
	  Const[mfst+n] = Plefte[nden]*Plefte[nst+n]*Plefte[nve2]; // composition
	}
	double css = Plefte[ncsp]*Plefte[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Plefte[nbm1]*Plefte[nbm1]  
                +Plefte[nbm2]*Plefte[nbm2]  
		+Plefte[nbm3]*Plefte[nbm3] )/Plefte[nden];

         Const[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Plefte[nbm2]*Plefte[nbm2]  //direction dependent
                                            /Plefte[nden])   
			      )/2.0e0);
         Const[mvel] = Plefte[nve2];//direction dependent
         Const[mpre] = ptl;
	 
