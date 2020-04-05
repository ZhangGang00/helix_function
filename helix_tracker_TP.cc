#include "FaserTracker/TrackCandidate.hh"
#include "FaserTracker/TrackFinder.hh"
#include "FaserTracker/TrackFitter.hh"
#include "FaserTracker/common_includes.hh"

#include "FaserTrackerGeometry.hh"
#include "FaserTrackerEvent.hh"
#include "FaserTrackerSpacePoint.hh"
#include "FaserTrackerTruthParticle.hh"

#include <cstdlib>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TGraph.h"

static const int N_plane=9; 
//static double target_z[8]={0, 50, 909.04, 959.04, 1009.04, 1059.04, 1942.02, 1992.02};
static double target_z[9]={8, 58, 108, 934, 984, 1034, 1892, 1942, 1992};
static std::vector<TVector3> selected_track_point;
static int charge=-1;
static const double B=0.6;
static const int target_pdgID=13;

//------------------------------------------------------------------------------

void usage(string executable) {
    cout << "Usage: " << executable << " [options]" << "\n"
         << "Options:\n"
         << "  -h/--help                      print help\n"
         << "  -i/--input FILE_NAME           path to input ROOT file\n"
         << "  -n/--nEntries N_ENTRIES        number of entries from input to process\n"
         << "\n";
}

//------------------------------------------------------------------------------

void helix_fit(int& npar, double* deriv, double& f, double par[], int flag)
{
    double R = par[0];
    double x0 = par[1];
    double y0 = par[2];
    double z0 = target_z[0];
    double phi = par[3];
    if (TMath::Pi() < phi) phi = phi - 2*TMath::Pi();  // keep the same definition of angle in the following code
    double sin_lambda = std::sin(par[4]);
    double cos_lambda = std::cos(par[4]);
    double Cx = x0 - R*std::cos(phi);  // (Cx, Cz) are the coordinates of center of the circle 
    double Cz = z0 - R*std::sin(phi);

    double target_s=0,target_x=0,target_y=0;
    double chi2=0, error_x=0.08/(sqrt(24)*std::cos(0.02)), error_y=0.08/(sqrt(24)*std::sin(0.02));
    int ipoint=0;
    for (int i=0; i<selected_track_point.size(); i++) 
      {
    if (charge==1)  target_s = (std::asin((selected_track_point[i].Z() - Cz)/R) - phi)*R/(charge*cos_lambda);
    else if (std::sin(phi) >0) target_s = ((TMath::Pi()-std::asin((selected_track_point[i].Z() - Cz)/R)) - phi)*R/(charge*cos_lambda);
    else target_s = ((-TMath::Pi()-std::asin((selected_track_point[i].Z() - Cz)/R)) - phi)*R/(charge*cos_lambda);

    if (target_s<-0.001) 
    { 
      cout<<"Error!!!  Negtive target_s !!!"<<std::endl; 
      chi2=1e10; 
      if (charge==1) cout<<std::setprecision(9)<<"R="<<R<<" x0="<<x0<<" y0="<<y0<<" phi0="<<phi<<" cos(lambda)="<<cos_lambda<<"  Cx="<<Cx<<"  Cz="<<Cz<<"  target_z="<<selected_track_point[i].Z()<<"  std::asin((target_z - Cz)/R)="<<std::asin((selected_track_point[i].Z() - Cz)/R)<<std::endl;
      else if (std::sin(phi) >0) cout<<std::setprecision(9)<<"R="<<R<<" x0="<<x0<<" y0="<<y0<<" phi0="<<phi<<" cos(lambda)="<<cos_lambda<<"  Cx="<<Cx<<"  Cz="<<Cz<<"  target_z="<<selected_track_point[i].Z()<<"  (TMath::Pi()-std::asin((target_z - Cz)/R))="<<(TMath::Pi()-std::asin((selected_track_point[i].Z() - Cz)/R))<<std::endl; 
      else cout<<std::setprecision(9)<<"R="<<R<<" x0="<<x0<<" y0="<<y0<<" phi0="<<phi<<" cos(lambda)="<<cos_lambda<<"  Cx="<<Cx<<"  Cz="<<Cz<<"  target_z="<<selected_track_point[i].Z()<<"  (-TMath::Pi()-std::asin((target_z - Cz)/R))="<<(-TMath::Pi()-std::asin((selected_track_point[i].Z() - Cz)/R))<<std::endl;
    }

    target_x = Cx + R*cos(phi+charge*target_s*cos_lambda/R); 
    target_y = y0 + target_s*sin_lambda;

         chi2 += (selected_track_point[i].X() - target_x)*(selected_track_point[i].X() - target_x)/(error_x*error_x) + (selected_track_point[i].Y() - target_y)*(selected_track_point[i].Y() - target_y)/(error_y*error_y);	   
      }

    f = chi2;
}

int helix(vector<double>& outvect, double px, double py, double pz, double x0, double y0, double z0, double B, int charge, double target_z, double s)
{
  if (target_z==-9999 && s==-9999) {cout<<"Error!!!   Please clarify what to calculate"<<std::endl; return -1; }
  
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
    else if (-vRz >0) phi = TMath::Pi() - std::asin(-vRz);
    else phi =  -TMath::Pi() - std::asin(-vRz);

  if (s != -9999) {
    outvect[0] = Cx + R*cos(phi+charge*s*cos_lambda/R);  //x 
    outvect[1] = y0 + s*sin_lambda;  //y
    outvect[2] = Cz + R*sin(phi+charge*s*cos_lambda/R);  //z
    //cout<<"SSS  CxCz"<<setprecision(9)<<Cx<<"  "<<Cz<<"  R="<<R<<"  phi="<<phi<<"  "<<charge*s*cos_lambda/R<<"  "<<sin(phi+charge*s*cos_lambda/R)<<"  "<<Cz + R*sin(phi+charge*s*cos_lambda/R)<<"  "<<outvect[2]<<std::endl;
    outvect[3] = -charge*p_perp*sin(phi+charge*s*cos_lambda/R); //px
    outvect[4] = py;  //py
    outvect[5] = charge*p_perp*cos(phi+charge*s*cos_lambda/R); //pz
  }
	
  if (target_z != -9999) {
    if (target_z-Cz > R) 
      {
	 cout<<"Warning!!!  target_z is outside the circle. target_z="<<target_z<<",  Cz="<<Cz<<",  R="<<R<<std::endl; 
	 return -1;
      }
    double target_s=0;
    if (charge==1)  target_s = (std::asin((target_z - Cz)/R) - phi)*R/(charge*cos_lambda);
    else if (-vRz >0) target_s = ((TMath::Pi()-std::asin((target_z - Cz)/R)) - phi)*R/(charge*cos_lambda);
    else target_s = ((-TMath::Pi()-std::asin((target_z - Cz)/R)) - phi)*R/(charge*cos_lambda);
    if (target_s<-0.001) {cout<<"Error!!!  Negtive target_s !!!"<<std::endl; return -1; }
    outvect[0] = Cx + R*cos(phi+charge*target_s*cos_lambda/R); 
    outvect[1] = y0 + target_s*sin_lambda;
    outvect[2] = target_z;
    outvect[3] = -charge*p_perp*sin(phi+charge*target_s*cos_lambda/R); 
    outvect[4] = py;
    outvect[5] = charge*p_perp*cos(phi+charge*target_s*cos_lambda/R); 
  }

  return 1; 
}

int main(int argc, char ** argv) {
    string input = "";
    long   nEntries = -1;

    int iOpt = 1;
    while (iOpt < argc) {
        string sw = argv[iOpt];
        if      (sw=="-h" || sw=="--help"      ) { usage(argv[0]); return 0; }
        else if (sw=="-i" || sw=="--input"     ) { input = argv[++iOpt]; }
        else if (sw=="-n" || sw=="--nEntries"  ) { nEntries = stoi(argv[++iOpt]); }
        else {
            cout << "Unknown switch: " << sw << "\n";
            usage(argv[0]);
            return 1;
        }
        ++iOpt;
    }
    cout << "Flags:\n"
         << "------\n"
         << "  nEntries   = " << nEntries   << "\n"
         << "  input      = " << input      << "\n"
         << "\n";


    TFile inputFile {input.c_str()};

    TTree * eventTree = (TTree*) inputFile.Get("events");
    if (eventTree == nullptr) {
        cout << "ERROR  Unable to load event tree: `events`\n"
             << "       Exiting.\n"
             << "\n";
        return 1;
    }
    cout << "INFO  Loaded event tree `events` from "
         << input << " with " << eventTree->GetEntries() << " entries.\n";
#if 0
    TTree * geoTree = (TTree*) inputFile.Get("geo");
    if (geoTree == nullptr) {
        cout << "ERROR  Unable to load event tree: `geo`\n"
             << "       Exiting.\n"
             << "\n";
        return 1;
    }
    long nGeoEntries = geoTree->GetEntries();
    cout << "INFO  Loaded input geometry tree `geo` from "
         << input << " with " << nGeoEntries << " entries.\n";
    if (nGeoEntries != 1) {
        cout << "ERROR  The geometry tree `geo` has " << nGeoEntries << " entries.\n"
             << "       It should have exactly one entry. Exiting.\n"
             << "\n";
        return 1;
    }

    FaserTrackerGeometry * geo = nullptr;
    geoTree->SetBranchAddress("geo", &geo);
    geoTree->GetEntry(0);
#endif

    TFile *outf = new TFile("/afs/cern.ch/user/g/gang/FASER/faser_tracker/helix_tracker_TP_output.root","recreate");
    TTree *mytree = new TTree("tracks","");
    double truth_px,truth_py,truth_pz,truth_E,truth_R,truth_phi0,truth_lambda;
    std::vector<double> truth_x;
    std::vector<double> truth_y;
    std::vector<double> truth_z;
    double reco_px,reco_py,reco_pz,reco_E,reco_R,reco_phi0,reco_lambda;
    std::vector<double> reco_x;
    std::vector<double> reco_y;
    std::vector<double> reco_z;
    double chi_square;
    int N_track_points,pdgID;
    std::vector<double> space_points_x;
    std::vector<double> space_points_y;
    std::vector<double> space_points_z;
    std::vector<double> track_points_x;
    std::vector<double> track_points_y;
    std::vector<double> track_points_z;
    std::vector<int> track_points_rowID;
    mytree->Branch("truth_px", &truth_px, "truth_px/D");
    mytree->Branch("truth_py", &truth_py, "truth_py/D");
    mytree->Branch("truth_pz", &truth_pz, "truth_pz/D");
    mytree->Branch("truth_E", &truth_E, "truth_E/D");
    mytree->Branch("truth_R", &truth_R, "truth_R/D");
    mytree->Branch("truth_phi0", &truth_phi0, "truth_phi0/D");
    mytree->Branch("truth_lambda", &truth_lambda, "truth_lambda/D");
    mytree->Branch("truth_x", &truth_x);
    mytree->Branch("truth_y", &truth_y);
    mytree->Branch("truth_z", &truth_z);
    mytree->Branch("reco_px", &reco_px, "reco_px/D");
    mytree->Branch("reco_py", &reco_py, "reco_py/D");
    mytree->Branch("reco_pz", &reco_pz, "reco_pz/D");
    mytree->Branch("reco_E", &reco_E, "reco_E/D");
    mytree->Branch("reco_R", &reco_R, "reco_R/D");
    mytree->Branch("reco_phi0", &reco_phi0, "reco_phi0/D");
    mytree->Branch("reco_lambda", &reco_lambda, "reco_lambda/D");
    mytree->Branch("reco_x", &reco_x);
    mytree->Branch("reco_y", &reco_y);
    mytree->Branch("reco_z", &reco_z);
    mytree->Branch("space_points_x", &space_points_x);
    mytree->Branch("space_points_y", &space_points_y);
    mytree->Branch("space_points_z", &space_points_z);
    mytree->Branch("track_points_x", &track_points_x);
    mytree->Branch("track_points_y", &track_points_y);
    mytree->Branch("track_points_z", &track_points_z);
    mytree->Branch("chi_square", &chi_square, "chi_square/D");
    mytree->Branch("N_track_points",&N_track_points,"N_track_points/I");
    mytree->Branch("pdgID",&pdgID,"pdgID/I");
    mytree->Branch("track_points_rowID",&track_points_rowID);

    FaserTrackerEvent * event = nullptr;
    eventTree->SetBranchAddress("events", &event);
    double hit[N_plane][1000]={0}; //hit[i][0] is the hits in corresponding layer, hit[i][j>0] are coordinates
    double hit_n1[N_plane][1000]={0};
    double hit_n2[N_plane][1000]={0};
    double hit_frac1[N_plane][1000]={0};
    double hit_frac2[N_plane][1000]={0};
    double hit_pos1[N_plane][1000]={0};
    double hit_pos2[N_plane][1000]={0};
    int hit_rowID[N_plane][5]={0};
    int N_hit_plane=0;
    TVector3 shift_pos={0,0,0};  // For the area without B filed
    double subtract_z=0;
    std::vector <double> outvect(6,0);

    int trackID_minus=1;
    for (int iEntry=0; iEntry<eventTree->GetEntries(); ++iEntry) 
    //for (int iEntry=0; iEntry<50; ++iEntry) 
      {
        eventTree->GetEntry(iEntry);
        cout << "INFO  Processing event " << event->eventNumber << "\n";
#if 1	
	bool cal_stat=1, find_truth=0;

	truth_x.clear(); truth_y.clear(); truth_z.clear();
        reco_x.clear(); reco_y.clear(); reco_z.clear(); 	
        space_points_x.clear(); space_points_y.clear(); space_points_z.clear();
        track_points_x.clear(); track_points_y.clear(); track_points_z.clear();
	track_points_rowID.clear();
        selected_track_point.clear();
        
	N_hit_plane=0;
	for (int i=0; i<N_plane; i++) 
	  { 
	    for (int j=0; j<1000; j++) 
	      {
		 hit[i][j]=0; 
		 hit_pos1[i][j]=0; 
		 hit_pos2[i][j]=0; 
                 hit_n1[i][j]=0;
                 hit_n2[i][j]=0;
                 hit_frac1[i][j]=0;
                 hit_frac2[i][j]=0;
	      }
	  }

        for (const FaserTrackerSpacePoint * sp : event->spacePoints) 
	{
            const TVector3 & pos = sp->globalPos;

            cout << "      Read space point with global position (" << pos.X() << "," << pos.Y() << "," << pos.Z() << ")\n";

	    for (int i=0; i<N_plane; i++) { 
	       if (std::fabs(pos.Z()-target_z[i]) < 20) 
	        {
		   hit[i][3*(int)hit[i][0]+1]=pos.X();
		   hit[i][3*(int)hit[i][0]+2]=pos.Y();
		   hit[i][3*(int)hit[i][0]+3]=pos.Z();
		   if (hit[i][0]==0) N_hit_plane++;
		   int N_th=0; 
		   double total_E1=0, total_E2=0;
                   for (uint icl : sp->analogClusterIndices) {
		      FaserTrackerCluster * cl = event->analogClusters[icl];
                      for (uint idigit : cl->digitIndices) {
			 FaserTrackerDigit * digit = event->digits[idigit];
                         for (uint iHit : digit->truthHitIndices) 
                         //for (uint iHit : sp->truthHitIndices) 
			   {
                             FaserTrackerTruthHit * th = event->truthHits[iHit];
                             const TVector3 & thpos = th->globalPos;
                             cout<<"truth hits:   trackID = "<<th->trackID<<"  coordinate = "<<thpos.X()<<"  "<<thpos.Y()<<"  "<<thpos.Z()<<"\n";

			     if (th->trackID==trackID_minus)
			      {
				hit_n1[i][(int)hit[i][0]]++;
				hit_pos1[i][3*(int)hit[i][0]+1] += th->energy * thpos.X(); 
				hit_pos1[i][3*(int)hit[i][0]+2] += th->energy * thpos.Y(); 
				hit_pos1[i][3*(int)hit[i][0]+3] += th->energy * thpos.Z(); 
				N_th++;
				total_E1 += th->energy; 
			      }
			     
	                   }
	              }
		   }
		 
		   if (N_th > 0) { 
                   hit_frac1[i][(int)hit[i][0]] = hit_n1[i][(int)hit[i][0]]/(double)N_th;
                   hit_frac2[i][(int)hit[i][0]] = hit_n2[i][(int)hit[i][0]]/(double)N_th;
		   } else {
                   hit_frac1[i][(int)hit[i][0]] = -9999;
                   hit_frac2[i][(int)hit[i][0]] = -9999;
		   }
		   if (total_E1 > 0) {
		   hit_pos1[i][3*(int)hit[i][0]+1] /= total_E1; 
		   hit_pos1[i][3*(int)hit[i][0]+2] /= total_E1; 
		   hit_pos1[i][3*(int)hit[i][0]+3] /= total_E1; 
		   } else {
		   hit_pos1[i][3*(int)hit[i][0]+1] = -9999; 
		   hit_pos1[i][3*(int)hit[i][0]+2] = -9999; 
		   hit_pos1[i][3*(int)hit[i][0]+3] = -9999; 
		   }
		   
		   hit[i][0]++;
		} 
	    }
        }
    
        if (N_hit_plane<4) continue;  // discard the tracks with less than 4 space points

        for (int i=0; i<N_plane; i++)  
	  {
	     if (hit[i][0] > 0) {
		int temp1=0, temp2=0;
		for (int j=1; j<hit[i][0]; j++) {
                   if ((hit_n1[i][j] > hit_n1[i][temp1]) || (hit_n1[i][j]==hit_n1[i][temp1] && hit_frac1[i][j] > hit_frac1[i][temp1])) 
		    {
		     temp1=j;
		    }
                   
		} 
#if 1
		for (int j=0; j<hit[i][0]; j++) { cout<<"plane = "<<i<<"   (x,y,z) = ("<<hit[i][3*j+1]<<", "<<hit[i][3*j+2]<<", "<<hit[i][3*j+3]<<")"<<"  hit_n1 = "<<hit_n1[i][j]<<"  "<<"hit_n2 = "<<hit_n2[i][j]<<"     hit_frac1 = "<<hit_frac1[i][j]<<"  "<<"hit_frac2 = "<<hit_frac2[i][j]<<"\n"; }
		cout<<"temp1 = "<<temp1<<"   temp2 = "<<temp2<<"\n";
		cout<<"hit_n1 = "<<hit_n1[i][temp1]<<"     "<<"hit_n2 = "<<hit_n2[i][temp2]<<"\n";
		cout<<"hit_frac1 = "<<hit_frac1[i][temp1]<<"      "<<"hit_frac2 = "<<hit_frac2[i][temp2]<<"\n";
		cout<<"hit_pos1 = ("<<hit_pos1[i][3*temp1+1]<<", "<<hit_pos1[i][3*temp1+2]<<", "<<hit_pos1[i][3*temp1+3]<<")     "<<"hit_pos2 = ("<<hit_pos2[i][3*temp2+1]<<", "<<hit_pos2[i][3*temp2+2]<<", "<<hit_pos2[i][3*temp2+3]<<")\n";
		cout<<"distance = "<<sqrt(pow(hit_pos1[i][3*temp1+1]-hit_pos2[i][3*temp2+1],2)+pow(hit_pos1[i][3*temp1+2]-hit_pos2[i][3*temp2+2],2))<<"\n\n";
#endif
	       if (hit_n1[i][temp1]>0) {
		target_z[i]=hit[i][3*temp1+3];
	        space_points_x.push_back(hit[i][3*temp1+1]);
	        space_points_y.push_back(hit[i][3*temp1+2]);
	        space_points_z.push_back(hit[i][3*temp1+3]);
		TVector3 tempvec(hit[i][3*temp1+1], hit[i][3*temp1+2], hit[i][3*temp1+3]);
		selected_track_point.push_back(tempvec); 
		track_points_x.push_back(hit[i][3*temp1+1]);
		track_points_y.push_back(hit[i][3*temp1+2]);
		track_points_z.push_back(hit[i][3*temp1+3]);
	       }
	     }
	  }

	N_track_points = selected_track_point.size();
	if (N_track_points < 4) continue;

	for (const FaserTrackerTruthParticle * tp : event->truthParticles) 
	  {
	    if (tp->pdgCode==target_pdgID && tp->parentID==0)
	     {
	       pdgID = tp->pdgCode;
	       find_truth=1;
               const TVector3 & ver = tp->vertex;
	       const TLorentzVector & fmom = tp->fourMomentum;
	       for (int iz=0; iz<N_plane; iz++) 
	         {
	          if (helix(outvect,fmom.Px()/1000.0,fmom.Py()/1000.0,fmom.Pz()/1000.0,ver.X(),ver.Y(),ver.Z(),B,charge,target_z[iz],-9999) == -1) cal_stat=0;
	          truth_x.push_back(outvect[0]);  truth_y.push_back(outvect[1]);  truth_z.push_back(outvect[2]);
                  if (iz==0) {truth_px=outvect[3]; truth_py=outvect[4]; truth_pz=outvect[5]; truth_E=fmom.E()/1000.0;}
	         }
	       int truth_charge=-1;
	       if (pdgID==target_pdgID) truth_charge=1;
	       truth_charge=charge;
	       double truth_p_perp = sqrt(truth_px*truth_px + truth_pz*truth_pz);  // The magnet field point to positive y axis
	       double truth_sin_lambda = truth_py/sqrt(truth_p_perp*truth_p_perp + truth_py*truth_py);
	       truth_lambda = std::asin(truth_sin_lambda); 
	       truth_R = 1000*truth_p_perp/(0.3*B);  // R(mm), p(GeV), B(T)
	       double truth_vRx = -truth_charge*truth_pz/truth_p_perp;  // vR is a norm vecter pointing from intial position to the center of circle
	       double truth_vRz = truth_charge*truth_px/truth_p_perp;
	       if (-truth_vRx > 0) truth_phi0 = std::asin(-truth_vRz);
                 else if (-truth_vRz >0) truth_phi0 = TMath::Pi() - std::asin(-truth_vRz);
                        else truth_phi0 =  -TMath::Pi() - std::asin(-truth_vRz);
	     }
	  }
	if (!find_truth) continue;
	if (!cal_stat) continue;


	//if (0.5*(selected_track_point[0].X()+selected_track_point[N_hit_plane-1].X()) < selected_track_point[N_hit_plane/2].X()) charge=1;
	//  else charge=-1;
            
        const int npar = 5;              // the number of parameters
        TMinuit minuit(npar);
        minuit.SetFCN(helix_fit);
	// minuit.DefineParameter(i, parName[i].c_str(), par_initial[i], stepSize[i], minVal[i], maxVal[i])
#if 0
	minuit.DefineParameter(0,"R",1e6,1e4,1e4,1e8);
	minuit.DefineParameter(1,"x0",0,0.01,-100,100);
	minuit.DefineParameter(2,"y0",0,0.01,-100,100);
	if (charge==1) minuit.DefineParameter(3,"phi0",0,0.0001,-TMath::Pi()/2,TMath::Pi()/2);
          else minuit.DefineParameter(3,"phi0",3*TMath::Pi()/4,0.0001,TMath::Pi()/2,3*TMath::Pi()/2);
	minuit.DefineParameter(4,"lambda",0,0.0001,-TMath::Pi()/2,TMath::Pi()/2);
#endif
// use truth info as initial parameters
#if 1
	minuit.DefineParameter(0,"R",truth_R,1e4,1e4,1e8);
	minuit.DefineParameter(1,"x0",truth_x[0],0.01,-100,100);
	minuit.DefineParameter(2,"y0",truth_y[0],0.01,-100,100);
	if (charge==1) minuit.DefineParameter(3,"phi0",truth_phi0,0.0001,-TMath::Pi()/2,TMath::Pi()/2);
          else minuit.DefineParameter(3,"phi0",truth_phi0,0.0001,TMath::Pi()/2,3*TMath::Pi()/2);
	minuit.DefineParameter(4,"lambda",truth_lambda,0.0001,-TMath::Pi()/2,TMath::Pi()/2);
#endif
	minuit.SetErrorDef(1);             
        minuit.SetMaxIterations(1000000);
        minuit.SetPrintLevel(1); //-1 no output, 1 standard output
	int fit_status=minuit.Migrad();
	if (fit_status == 4) cout<<"\n"<<"  WARNING!!!  Migrad not converged !!!"<<"\n"<<"\n";
        // Print results
	double outpar[npar], err[npar];
        for (int i=0; i<npar; i++){
           minuit.GetParameter(i,outpar[i],err[i]);
           cout<<"par["<<i<<"] "<<outpar[i]<<"+/-"<<err[i]<<std::endl;
        }
        Double_t fmin,edm,errdef;
        Int_t nvpar,nparx,istat;
        minuit.mnstat(fmin,edm,errdef,nvpar,nparx,istat); 
	// istat: 0= not calculated at all, 1= approximation only, not accurate, 2= full matrix, but forced positive-definite, 3= full accurate covariance matrix
        cout<<"fmin="<<fmin<<"  edm="<<edm<<"  errdef="<<errdef<<"  nvpar="<<nvpar<<"  nparx="<<nparx<<"  istat="<<istat<<std::endl;
	chi_square = fmin;

	double mass=0.105;
	double p_perp=0.3*B*outpar[0]/1000.0;
        reco_px=-charge*p_perp*sin(outpar[3]);
	reco_py=p_perp*std::tan(outpar[4]);
	reco_pz=charge*p_perp*cos(outpar[3]);
	reco_E=std::sqrt(reco_px*reco_px+reco_py*reco_py+reco_pz*reco_pz+mass*mass);
	reco_R=outpar[0];
	reco_phi0=outpar[3];
	reco_lambda=outpar[4];
	reco_x.push_back(outpar[1]);  reco_y.push_back(outpar[2]);  reco_z.push_back(target_z[0]);
	for (int iz=1; iz<N_plane; iz++) {
	   if (helix(outvect,reco_px,reco_py,reco_pz,reco_x[0],reco_y[0],reco_z[0],B,charge,target_z[iz],-9999) == -1) cal_stat=0;
	   reco_x.push_back(outvect[0]);  reco_y.push_back(outvect[1]);  reco_z.push_back(outvect[2]);
	}
	if (!cal_stat) continue;

	mytree->Fill();
#endif
    }

    outf->cd();
    mytree->Write();
    outf->Close();
    //delete mytree;
    delete outf;
}

