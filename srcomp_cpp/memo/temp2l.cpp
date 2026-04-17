	  //direction dependent w,u,v
	Cons[mudn] = Plefte[nden]; // rho
	Cons[muvw] = Plefte[nve1]*Plefte[nden]; // rho v_x
	Cons[muvu] = Plefte[nve2]*Plefte[nden]; // rho v_y
	Cons[muvv] = Plefte[nve3]*Plefte[nden]; // rho v_z
        Cons[muet] = Plefte[nene]*Plefte[nden]  // e_i
	              +0.5e0*Plefte[nden]*(                  
                          +Plefte[nve1]*Plefte[nve1]                 
                          +Plefte[nve2]*Plefte[nve2]                 
                          +Plefte[nve3]*Plefte[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Plefte[nbm1]*Plefte[nbm1]                 
                          +Plefte[nbm2]*Plefte[nbm2]                 
                          +Plefte[nbm3]*Plefte[nbm3]); // + B^2/2

	Cons[mubw] = Plefte[nbm1]; // b_x
	Cons[mubu] = Plefte[nbm2]; // b_y
	Cons[mubv] = Plefte[nbm3]; // b_z
	Cons[mubp] = Plefte[nbps]; // psi
	for(int n=0; n<ncomp; n++){
	  Cons[must+n] = Plefte[nden]*Plefte[nst+n]; // composition
	}
	
        double  ptl = Plefte[npre] + ( Plefte[nbm1]*Plefte[nbm1]
                                      +Plefte[nbm2]*Plefte[nbm2]
				      +Plefte[nbm3]*Plefte[nbm3])/2.0e0;
	//direction dependent, nve2 or nbm2
	Cons[mfdn] = Plefte[nden]*Plefte[nve2];
	Cons[mfvw] = Plefte[nden]*Plefte[nve1]*Plefte[nve2] 
	                           -Plefte[nbm1]*Plefte[nbm2];
	Cons[mfvu] = Plefte[nden]*Plefte[nve2]*Plefte[nve2]
      	                     + ptl -Plefte[nbm2]*Plefte[nbm2];// p diagnonal
        Cons[mfvv] = Plefte[nden]*Plefte[nve3]*Plefte[nve2]
	                          - Plefte[nbm3]*Plefte[nbm2];
        Cons[mfet] = (Cons[muet]+ptl)*Plefte[nve2]
                           -( Plefte[nbm1]*Plefte[nve1]
                             +Plefte[nbm2]*Plefte[nve2]
			     +Plefte[nbm3]*Plefte[nve3])*Plefte[nbm2];

	// direction dependent 3, 2, 1, 2
	Cons[mfbu] =  0.0;
	Cons[mfbv] =  Plefte[nbm3]*Plefte[nve2]
	              - Plefte[nve3]*Plefte[nbm2];
	Cons[mfbw] =  Plefte[nbm1]*Plefte[nve2]
	              - Plefte[nve1]*Plefte[nbm2];
	Cons[mfbp] = 0.0e0;  // psi
	
	for(int n=0; n<ncomp; n++){
	  Cons[mfst+n] = Plefte[nden]*Plefte[nst+n]*Plefte[nve2]; // composition
	}
	double css = Plefte[ncsp]*Plefte[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Plefte[nbm1]*Plefte[nbm1]  
                +Plefte[nbm2]*Plefte[nbm2]  
		+Plefte[nbm3]*Plefte[nbm3] )/Plefte[nden];

         Cons[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Plefte[nbm2]*Plefte[nbm2]  //direction dependent
                                            /Plefte[nden])   
			      )/2.0e0);
         Cons[mvel] = Plefte[nve2];//direction dependent
         Cons[mpre] = ptl;
	 
