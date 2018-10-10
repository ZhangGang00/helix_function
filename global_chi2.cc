#define Pi 3.1415926

double target_z[8]={0,50,910,960,1010,1060,1942,1992}; // z coordinates of each layer
double points[8][2]={0}; // space points (x,z)
int charge=1;

void helix_fit(int& npar, double* deriv, double& f, double par[], int flag)
{
    double R = par[0];
    double x0 = par[1];
    double y0 = par[2];
    double z0 = 0;
    double phi = par[3];
    if (Pi < phi) phi = phi - 2*Pi;  // keep the same definition of angle in the following code
    double sin_lambda = std::sin(par[4]);
    double cos_lambda = std::cos(par[4]);
    double Cx = x0 - R*std::cos(phi);  // (Cx, Cz) are the coordinates of center of the circle 
    double Cz = z0 - R*std::sin(phi);

    double target_s=0,target_x=0,target_y=0;
    double chi2=0, error=0.08/sqrt(12);
    for (int i=0; i<8; i++) 
      {
        if (charge==1)  target_s = (std::asin((target_z[i] - Cz)/R) - phi)*R/(charge*cos_lambda);
        else if (std::sin(phi) >0) target_s = ((Pi-std::asin((target_z[i] - Cz)/R)) - phi)*R/(charge*cos_lambda);
        else target_s = ((-Pi-std::asin((target_z[i] - Cz)/R)) - phi)*R/(charge*cos_lambda);
        if (target_s<0) 
        { 
           cout<<"Error!!!  Negtive target_s !!!"<<endl; 
           chi2=1e20; 
           if (charge==-1) cout<<setprecision(9)<<R<<" "<<x0<<" "<<y0<<" "<<phi<<" "<<cos_lambda<<"  Cx="<<Cx<<"  Cz="<<Cz<<"  target_z="<<target_z[i]<<"  (Pi-std::asin((target_z[i] - Cz)/R))="<<(Pi-std::asin((target_z[i] - Cz)/R))<<endl; 
             else cout<<setprecision(9)<<R<<" "<<x0<<" "<<y0<<" "<<phi<<" "<<cos_lambda<<"  Cx="<<Cx<<"  Cz="<<Cz<<"  target_z="<<target_z[i]<<"  std::asin((target_z[i] - Cz)/R)="<<std::asin((target_z[i] - Cz)/R)<<endl;
        }
       target_x = Cx + R*cos(phi+charge*target_s*cos_lambda/R); 
       target_y = y0 + target_s*sin_lambda;

       chi2 += ((points[i][0] - target_x)*(points[i][0] - target_x) + (points[i][1] - target_y)*(points[i][1] - target_y)) / (error*error);	    
     }

   f = chi2;
}


void global_chi2(TString input) 
{
    
        const int npar = 5;              // the number of parameters
        TMinuit minuit(npar);
        minuit.SetFCN(helix_fit);
	// minuit.DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i])
	minuit.DefineParameter(0,"R",1e6,1e4,1e4,1e7);
	minuit.DefineParameter(1,"x0",0,1,-100,100);
	minuit.DefineParameter(2,"y0",0,1,-100,100);
	if (charge==1) minuit.DefineParameter(3,"phi0",0,0.01,-Pi/2,Pi/2);
          else minuit.DefineParameter(3,"phi0",3*Pi/4,0.01,Pi/2,3*Pi/2);
	minuit.DefineParameter(4,"lambda",0,0.01,-Pi,Pi);
	minuit.SetErrorDef(1);             
        minuit.SetMaxIterations(1000000);
        minuit.SetPrintLevel(1); //-1 no output, 1 standard output
	minuit.Migrad();
	double outpar[npar], err[npar];
        for (int i=0; i<npar; i++){
           minuit.GetParameter(i,outpar[i],err[i]);
           cout<<"par["<<i<<"] "<<outpar[i]<<"+/-"<<err[i]<<endl;
        }
        // Print results
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
        cout<<"amin="<<amin<<"  edm="<<edm<<"  errdef="<<errdef<<"  nvpar="<<nvpar<<"  nparx="<<nparx<<"  icstat="<<icstat<<endl;
              
	// plot the post-fit helix, using helix function
	vector <double> outvect(6,0);  
	int index=0;
	double si=0,xi=0,yi=0,zi=0;
	double B=0.5;
	double p_perp=0.3*B*outpar[0]/1000.0;
	double px=-charge*p_perp*sin(outpar[3]),py=p_perp*tan(outpar[4]),pz=charge*p_perp*cos(outpar[3]),x00=outpar[1],y00=outpar[2],z00=0;
        while (zi>=z00 && zi<zMax) {
	  helix(outvect,px,py,pz,x00,y00,z00,B,charge,-9999,si);
	  xi = outvect[0];
	  yi = outvect[1];
	  zi = outvect[2];
	  cout<<"pl3D2: si="<<si<<"  xi="<<xi<<"  yi="<<yi<<"  zi="<<zi<<"  px="<<outvect[3]<<"  py="<<outvect[4]<<"  pz="<<outvect[5]<<endl;
	  pl3D2->SetPoint(index,zi,xi,yi);
	  pl_xz->SetPoint(index,zi,xi);
	  pl_yz->SetPoint(index,zi,yi);
	  si += 1;
	  index++;
	}

}

