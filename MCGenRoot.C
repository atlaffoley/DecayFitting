//To be run within ROOT, creating histograms of each component of your decay
//To save the histograms to a root file after running, uncomment the last line

#include<stdio.h>
#include<ctype.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<sys/ioctl.h>
#include<sys/mtio.h>
#include<math.h>
#include<time.h>
#include<signal.h>
#include<stdlib.h>
#include<errno.h>
#include<string.h>
#include<unistd.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sys/time.h>
#include <TH1D.h>
using namespace std;
    
//Timing variable to change for your experiment
double bgtime = 2.0;
double beamon = 1.5;  //3 for T1/2, 20 for BR
double decaytime = 1.5; //15 for T1/2, 3 for BR

//Other timing variables, don't change unless absolutely necessary
double timestep = 1e-5; //1 us time during implant
double cycletime = bgtime + beamon + decaytime;
double bin = 0.01;

TH1D *Total = new TH1D ("Total", "Sum of all components", cycletime/bin, 0, cycletime);
TH1D *BG = new TH1D ("BG", "Background", cycletime/bin, 0, cycletime);
TH1D *A = new TH1D ("A", "Parent Component", cycletime/bin, 0, cycletime);
TH1D *B = new TH1D ("B", "Daughter Component (beta)", cycletime/bin, 0, cycletime);
TH1D *C = new TH1D ("C", "Daughter Component (beta-n)", cycletime/bin, 0, cycletime);
TH1D *D = new TH1D ("D", "Daughter Component (beta-2n)", cycletime/bin, 0, cycletime);
TH1D *E = new TH1D ("E", "GrandDaughter Component (beta-beta)", cycletime/bin, 0, cycletime);
TH1D *F = new TH1D ("F", "GrandDaughter Component (beta-beta-n)", cycletime/bin, 0, cycletime);
TH1D *G = new TH1D ("G", "GrandDaughter Component (beta-beta-2n)", cycletime/bin, 0, cycletime);
TH1D *H = new TH1D ("H", "GrandDaughter Component (beta-beta-3n)", cycletime/bin, 0, cycletime);
TH1D *I = new TH1D ("I", "GreatGrandDaughter Component (beta-beta-beta)", cycletime/bin, 0, cycletime);
TH1D *J = new TH1D ("J", "GreatGrandDaughter Component (beta-beta-beta-n)", cycletime/bin, 0, cycletime);
TH1D *K = new TH1D ("K", "GreatGrandDaughter Component (beta-beta-beta-2n)", cycletime/bin, 0, cycletime);
TH1D *L = new TH1D ("L", "GreatGrandDaughter Component (beta-beta-beta-3n)", cycletime/bin, 0, cycletime);
TH1D *M = new TH1D ("M", "GreatGreatGrandDaughter Component (beta-beta-beta-beta)", cycletime/bin, 0, cycletime);


void WriteToFile() {
    TFile out_file("MCoutput.root", "RECREATE");

    Total->Write();
    BG->Write();
    A->Write();
    B->Write();
    C->Write();
    D->Write();
    E->Write();
    F->Write();
    G->Write();
    H->Write();
    I->Write();
    J->Write();
    K->Write();
    L->Write();
    M->Write();
}

void MCGenRoot() {
	FILE *file, *myFile;
	char line[256], filename[128];

	vector<double>::iterator it, it2, it3;

	double timer;
	static unsigned long int outbufferlong[1024];

	unsigned int seed = 124819;//Seeds 0-9 used for the 10 "subruns"
	srand48(seed);

    //Intensity variables to edit
    double AbeamI = 3e6; //measured ISAC yield, delivered isotope A
    double BbeamI = 0; //measured ISAC yield, delivered isotope B
    double CbeamI = 0; //measured ISAC yield, delivered isotope C
    double DbeamI = 0; //measured ISAC yield, delivered isotope D
    double EbeamI = 0; //measured ISAC yield, delivered isotope E
    double FbeamI = 0; //measured ISAC yield, delivered isotope F
    double GbeamI = 0; //measured ISAC yield, delivered isotope G
    double HbeamI = 0; //measured ISAC yield, delivered isotope H
    double IbeamI = 0; //measured ISAC yield, delivered isotope I
    double JbeamI = 0; //measured ISAC yield, delivered isotope J
    double KbeamI = 0; //measured ISAC yield, delivered isotope K
    double LbeamI = 0; //measured ISAC yield, delivered isotope L
    double MbeamI = 0; //measured ISAC yield, delivered isotope M
    
    double bgrate_b = 150; // c/s beta background per cycle 
    int totcycles = 1; // increase as needed to get enough stats
    
    //Designed for the decay of 34Mg, which has lots of beta-n and beta-2n branches going down to great-great-granddaughter before stability
    //In practice, we can set lots of our beta-n and beta-2n probabilities to 0 and "kill" after a certain number of generations
    //If the half-life is set to <=0, we will consider it to be stable
    double par_half = 0.020; // parent half-life (A)
    double daub_half = 0.06; // beta daughter half-life (B = A->beta)
    double daubn_half = 0.15; // beta-n daughter half-life (C = A->beta+n)
    double daub2n_half = 0.3; // beta-2n daughter half-life (D = A->beta+2n)
    double gdb_half  = 2.; // 2xbeta granddaughter half-life (E = B->beta)
    double gdbn_half = 1.11; // 2xbeta-1xn granddaughter half-life (F = B->beta+n = C->beta)
    double gdb2n_half = 15; // 2xbeta-2xn granddaughter half-life (G = C->beta-n = D->beta)
    double gdb3n_half = 1.0; // 2xbeta-3xn granddaughter half-life (H = D->beta-n)
    double ggdb_half  = 1.0; // 3xbeta greatgranddaughter half-life (I = E->beta)
    double ggdbn_half = 1.0; //3xbeta-1xn greatgranddaughter half-life (J = F->beta)
    double ggdb2n_half = 1.0; //3xbeta-2xn greatgranddaughter half-life (K = G->beta)
    double ggdb3n_half = 1.0; //3xbeta-3xn greatgranddaughter half-life (L = H->beta)
    double gggdb_half = 1.0; //4xbeta greatgreatgranddaughter half-life (M = I->beta

    double pn_par = 0.25; 
    double p2n_par = 0.25;
    double pn_daub = 0.25; 
    double p2n_daub = 0.25;
    double pn_daubn = 0.25;
    double p2n_daubn = 0.00;

	static unsigned long int comp_hits[1024], sum_hits[1024];
    int holder;
    double time_par, time_daub, time_gdb;
    double time_g1, time_g2, time_g3;
    double BR_g1, BR_g2;
    
    double value;
   
    //Fill in background counts
    int NBG_b = bgrate_b * totcycles * (bgtime + beamon + decaytime); // number of bg counts
    for(int i = 0; i < NBG_b; i++) {
        value = drand48()*cycletime;
        BG->Fill(value);
        Total->Fill(value);
    }

    for(double timer=0; timer<=beamon; timer+=timestep) {
        printf("Working ... beam on time %lf of %lf (~%d%%)\r", timer,beamon,(int)(timer/beamon*100));

        //Loop over initial beam (parent, A)
        for(int m = AbeamI*totcycles*timestep; m>0; m--) {
            //Initialize times to -1 so they go in underflow bin if we don't modify them
            time_g1 = -1;
            time_g2 = -1;
            time_g3 = -1;
            BR_g1 = 1;
            BR_g2 = 1;

            time_g1 = -log(drand48())*par_half/log(2.0); //randomize how long until generation 1 (A) decays
            time_g1 += bgtime + timer; //offset by where we actually are in time
            A->Fill(time_g1);
            Total->Fill(time_g1);
            
            BR_g1 = drand48(); //determine where we go to
            BR_g2 = drand48(); //determine where we go to
           
            //1st generation of decay (daughters)
            if(BR_g1 < (1 - pn_par - p2n_par)) {//just beta decay
                time_g2 = -log(drand48())*daub_half/log(2.0);
                time_g2 += time_g1;
                B->Fill(time_g2);
                Total->Fill(time_g2);
                //2nd generation of decay (granddaughters)
                if(BR_g2 < (1 - pn_daub - p2n_daub)) {//just beta decay
                    time_g3 = -log(drand48())*gdb_half/log(2.0);
                    time_g3 += time_g2;
                    E->Fill(time_g3);
                    Total->Fill(time_g3);
                }
                else if(BR_g2 < (1 - p2n_daub) && BR_g2 > (1 - pn_daub - p2n_daub)) { //beta-n decay
                    time_g3 = -log(drand48())*gdbn_half/log(2.0);
                    time_g3 += time_g2;
                    F->Fill(time_g3);
                    Total->Fill(time_g3);
                }
                else {
                    time_g3 = -log(drand48())*gdb2n_half/log(2.0);
                    time_g3 += time_g2;
                    G->Fill(time_g3);
                    Total->Fill(time_g3);
                }

            }
            else if(BR_g1 < (1 - p2n_par) && BR_g1 > (1 - pn_par - p2n_par)) { //beta-n decay
                time_g2 = -log(drand48())*daubn_half/log(2.0);
                time_g2 += time_g1;
                C->Fill(time_g2);
                Total->Fill(time_g2);
                //2nd generation of decay (granddaughters)
                if(BR_g2 < (1 - pn_daubn - p2n_daubn)) {//just beta decay
                    time_g3 = -log(drand48())*gdbn_half/log(2.0);
                    time_g3 += time_g2;
                    F->Fill(time_g3);
                    Total->Fill(time_g3);
                }
                else if(BR_g2 < (1 - p2n_daubn) && BR_g2 > (1 - pn_daubn - p2n_daubn)) { //beta-n decay
                    time_g3 = -log(drand48())*gdb2n_half/log(2.0);
                    time_g3 += time_g2;
                    G->Fill(time_g3);
                    Total->Fill(time_g3);
                }
                //we are not following beta-2n decay of the beta-n daugher (3 fewer nucleons than parent)
            }
            else { //beta-2n decay
                time_g2 = -log(drand48())*daub2n_half/log(2.0);
                time_g2 += time_g1;
                D->Fill(time_g2);
                Total->Fill(time_g2);
                //2nd generation of decay (granddaughters)
                //We are not considering any beta-n or beta-2n from the beta-2n decay of the parent
                time_g3 = -log(drand48())*gdb2n_half/log(2.0);
                time_g3 += time_g2;
                G->Fill(time_g3);
                Total->Fill(time_g3);
            }
        }
    }
    cout << endl;

    WriteToFile();
}//end main
