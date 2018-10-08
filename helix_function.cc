int helix(vector <double>& outvect, double px, double py, double pz, double x0, double y0, double z0, double B, int charge, double target_z, double s)
{
  if (target_z==-9999 && s==-9999) {cout<<"Error!!!   Please clarify what to calculate"<<endl; return -1; }
  
  double p_perp = sqrt(px*px + pz*pz);  // The magnet field point to positive y axis
  double sin_lambda = py/sqrt(p_perp*p_perp + py*py); 
  double cos_lambda = p_perp/sqrt(p_perp*p_perp + py*py); 
  double R = 1000*p_perp/(0.3*B);  // R(mm), p(GeV), B(T)
  double vRx = -charge*pz/p_perp;  // vR is a norm vecter pointing from intial position to the center of the circle
  double vRz = charge*px/p_perp;
  double Cx = x0 + R*vRx;  // (Cx, Cz) are the coordinates of center of the circle  
  double Cz = z0 + R*vRz;
  double phi=0;  // -vR determines phi
  if (-vRx > 0) phi = std::asin(-vRz);
    else if (-vRz >0) phi = Pi - std::asin(-vRz);
    else phi =  -Pi - std::asin(-vRz);

  if (s != -9999) {
    outvect[0] = Cx + R*cos(phi+charge*s*cos_lambda/R);  //x 
    outvect[1] = y0 + s*sin_lambda;  //y
    outvect[2] = Cz + R*sin(phi+charge*s*cos_lambda/R);  //z
    outvect[3] = -charge*p_perp*sin(phi+charge*s*cos_lambda/R); //px
    outvect[4] = py;  //py
    outvect[5] = charge*p_perp*cos(phi+charge*s*cos_lambda/R); //pz
  }
	
  if (target_z != -9999) {
    if (target_z-Cz > R) 
      {
	 cout<<"Warning!!!  target_z is outside the circle. target_z="<<target_z<<",  Cz="<<Cz<<",  R="<<R<<endl; 
	 return -1;
      }
    double target_s=0;
    if (charge==1)  target_s = (std::asin((target_z - Cz)/R) - phi)*R/(charge*cos_lambda);
    else if (-vRz >0) target_s = ((Pi-std::asin((target_z - Cz)/R)) - phi)*R/(charge*cos_lambda);
    else target_s = ((-Pi-std::asin((target_z - Cz)/R)) - phi)*R/(charge*cos_lambda);
    if (target_s<0) {cout<<"Error!!!  Negtive target_s !!!"<<endl; return -1; }
    outvect[0] = Cx + R*cos(phi+charge*target_s*cos_lambda/R); 
    outvect[1] = y0 + target_s*sin_lambda;
    outvect[2] = target_z;
    outvect[3] = -charge*p_perp*sin(phi+charge*target_s*cos_lambda/R); 
    outvect[4] = py;
    outvect[5] = charge*p_perp*cos(phi+charge*target_s*cos_lambda/R); 
  }

  return 0; 
}
