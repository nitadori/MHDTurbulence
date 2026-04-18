	  //direction dependent w,u,v
	Cons[mudn] = Prim[nden]; // rho
	Cons[muvw] = Prim[nve1]*Prim[nden]; // rho v_x
	Cons[muvu] = Prim[nve2]*Prim[nden]; // rho v_y
	Cons[muvv] = Prim[nve3]*Prim[nden]; // rho v_z
        Cons[muet] = Prim[nene]*Prim[nden]  // e_i
	              +0.5e0*Prim[nden]*(                  
                          +Prim[nve1]*Prim[nve1]                 
                          +Prim[nve2]*Prim[nve2]                 
                          +Prim[nve3]*Prim[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Prim[nbm1]*Prim[nbm1]                 
                          +Prim[nbm2]*Prim[nbm2]                 
                          +Prim[nbm3]*Prim[nbm3]); // + B^2/2

	Cons[mubw] = Prim[nbm1]; // b_x
	Cons[mubu] = Prim[nbm2]; // b_y
	Cons[mubv] = Prim[nbm3]; // b_z
	Cons[mubp] = Prim[nbps]; // psi
	for(int n=0; n<ncomp; n++){
	  Cons[must+n] = Prim[nden]*Prim[nst+n]; // composition
	}
	
        double  ptl = Prim[npre] + ( Prim[nbm1]*Prim[nbm1]
                                      +Prim[nbm2]*Prim[nbm2]
				      +Prim[nbm3]*Prim[nbm3])/2.0e0;
	//direction dependent, nve2 or nbm2
	Cons[mfdn] = Prim[nden]*Prim[nve2];
	Cons[mfvw] = Prim[nden]*Prim[nve1]*Prim[nve2] 
	                           -Prim[nbm1]*Prim[nbm2];
	Cons[mfvu] = Prim[nden]*Prim[nve2]*Prim[nve2]
      	                     + ptl -Prim[nbm2]*Prim[nbm2];// p diagnonal
        Cons[mfvv] = Prim[nden]*Prim[nve3]*Prim[nve2]
	                          - Prim[nbm3]*Prim[nbm2];
        Cons[mfet] = (Cons[muet]+ptl)*Prim[nve2]
                           -( Prim[nbm1]*Prim[nve1]
                             +Prim[nbm2]*Prim[nve2]
			     +Prim[nbm3]*Prim[nve3])*Prim[nbm2];

	// direction dependent 3, 2, 1, 2
	Cons[mfbu] =  0.0;
	Cons[mfbv] =  Prim[nbm3]*Prim[nve2]
	              - Prim[nve3]*Prim[nbm2];
	Cons[mfbw] =  Prim[nbm1]*Prim[nve2]
	              - Prim[nve1]*Prim[nbm2];
	Cons[mfbp] = 0.0e0;  // psi
	
	for(int n=0; n<ncomp; n++){
	  Cons[mfst+n] = Prim[nden]*Prim[nst+n]*Prim[nve2]; // composition
	}
	double css = Prim[ncsp]*Prim[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Prim[nbm1]*Prim[nbm1]  
                +Prim[nbm2]*Prim[nbm2]  
		+Prim[nbm3]*Prim[nbm3] )/Prim[nden];

         Cons[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Prim[nbm2]*Prim[nbm2]  //direction dependent
                                            /Prim[nden])   
			      )/2.0e0);
         Cons[mvel] = Prim[nve2];//direction dependent
         Cons[mpre] = ptl;
	 
