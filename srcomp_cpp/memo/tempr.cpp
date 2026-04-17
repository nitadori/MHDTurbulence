	Cons[mudn] = Prim[nden]; // rho
	Cons[muvu] = Prim[nve1]*Prim[nden]; // rho v_x
	Cons[muvv] = Prim[nve2]*Prim[nden]; // rho v_y
	Cons[muvw] = Prim[nve3]*Prim[nden]; // rho v_z
        Cons[muet] = Prim[nene]*Prim[nden]  // e_i
	              +0.5e0*Prim[nden]*(                  
                          +Prim[nve1]*Prim[nve1]                 
                          +Prim[nve2]*Prim[nve2]                 
                          +Prim[nve3]*Prim[nve3])  // + rho v^2/2
                      +0.5e0*             (                 
                          +Prim[nbm1]*Prim[nbm1]                 
                          +Prim[nbm2]*Prim[nbm2]                 
                          +Prim[nbm3]*Prim[nbm3]); // + B^2/2

	Cons[mubu] = Prim[nbm1]; // b_x
	Cons[mubv] = Prim[nbm2]; // b_y
	Cons[mubw] = Prim[nbm3]; // b_z
	Cons[mubp] = Prim[nbps]; // psi
	for(int n=0; n<ncomp;n++){
	  Cons[must+n] = Prim[nden]*Prim[nst+n]; // composition
	}
	// total pressure
        double  ptl = Prim[npre] + ( Prim[nbm1]*Prim[nbm1]
			+Prim[nbm2]*Prim[nbm2]
			+Prim[nbm3]*Prim[nbm3])/2.0e0;

	Cons[mfdn] = Prim[nden]*Prim[nve1];
	Cons[mfvu] = Prim[nden]*Prim[nve1]*Prim[nve1] 
	                      + ptl-Prim[nbm1]*Prim[nbm1];
	Cons[mfvv] = Prim[nden]*Prim[nve2]*Prim[nve1]
	                           -Prim[nbm2]*Prim[nbm1];
        Cons[mfvw] = Prim[nden]*Prim[nve3]*Prim[nve1]
	                           -Prim[nbm3]*Prim[nbm1];
        Cons[mfet] = (Cons[muet]+ptl)*Prim[nve1]
                           -( Prim[nbm1]*Prim[nve1]
                             +Prim[nbm2]*Prim[nve2]
			     +Prim[nbm3]*Prim[nve3])*Prim[nbm1];

	Cons[mfbu] =  0.0e0;
	Cons[mfbv] =  Prim[nbm2]*Prim[nve1]
	               -Prim[nve2]*Prim[nbm1];
	Cons[mfbw] =  Prim[nbm3]*Prim[nve1]
	               -Prim[nve3]*Prim[nbm1];
	Cons[mfbp] = 0.0e0;  // psi
     
	for(int n=0; n<ncomp;n++){
	  Cons[mfst+n] = Prim[nden]*Prim[nst+n]*Prim[nve1]; // composition
	}
	double css = Prim[ncsp]*Prim[ncsp];
        double cts =  css // c_s^2*c_a^2;
	     + ( Prim[nbm1]*Prim[nbm1]  
                +Prim[nbm2]*Prim[nbm2]  
		+Prim[nbm3]*Prim[nbm3] )/Prim[nden];

         Cons[mcsp] = sqrt((cts +sqrt(cts*cts
                                  -4.0e0*css*Prim[nbm1]*Prim[nbm1] 
                                            /Prim[nden])   
			      )/2.0e0);
         Cons[mvel] = Prim[nve1]; //direction dependent
         Cons[mpre] = ptl;
