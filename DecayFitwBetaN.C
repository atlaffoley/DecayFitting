//Double_t corr = 0.0001/8560.;
Double_t corr = 0.000;
Double_t pbn  = 0.58;     // parent (34Mg) beta-n branching
Double_t pb2n  = 0.05;    // parent (34Mg) beta-2n branching
Double_t dbn  = 0.26;     // daughter (34Al) beta-n branching
Double_t db2n  = 0.01;    // daughter (34Al) beta-2n branching
Double_t dnbn = 0.085;     // beta-n branching of the first beta-n daughter(33Al)
Double_t dnb2n = 0.00;    // beta-2n branching of the first beta-n daughter(33Al)

////////// Parent Functions
Double_t Implantation(const double & intens, const double & lambda, const double & time){
   return intens*(1. - TMath::Exp(-lambda*time));   // what is intens, beam intensity?? do we see grow in parent (beam) in this function return?
}

////////// Daughter grow-in during beam on
Double_t GrowInImplantationDaughter(const double & parent_intens, const double & parent_lambda, const double & daughter_lambda, const double & time){
   double term1 = parent_intens*daughter_lambda/(daughter_lambda-parent_lambda)*(TMath::Exp(-daughter_lambda*time)-TMath::Exp(-parent_lambda*time));  
   double term2 = parent_intens*(1.-TMath::Exp(-daughter_lambda*time));  

   return term1 + term2;   
}

///////// Daughter grow-in during beam off
Double_t GrowInDecayDaughter(const double & parent_initial_amount, const double & parent_lambda, const double & daughter_lambda, const double & time){
   return parent_initial_amount*daughter_lambda/(daughter_lambda-parent_lambda)*(TMath::Exp(-parent_lambda*time)-TMath::Exp(-daughter_lambda*time));
}  
   
//////////// Granddaughter Functions
Double_t GrowInImplantationGrandD(const double & parent_intens, const double & parent_lambda, const double & daughter_lambda, const double & grandD_lambda, const double & time){
   double term1 = parent_intens*parent_lambda*grandD_lambda*TMath::Exp(-daughter_lambda*time)/(grandD_lambda-daughter_lambda)/(daughter_lambda-parent_lambda);
   double term2 = parent_intens*daughter_lambda*grandD_lambda/(daughter_lambda-parent_lambda)/(parent_lambda-grandD_lambda)*TMath::Exp(-parent_lambda*time);
   double term3 = parent_intens;
   double term4 = parent_intens*parent_lambda*daughter_lambda/(daughter_lambda-grandD_lambda)/(grandD_lambda-parent_lambda)*TMath::Exp(-grandD_lambda*time);

   return term1 + term2 + term3 + term4;   
}

Double_t GrowInDecayGrandD(const double & parent_intens, const double & daughter_start_amount, const double & parent_lambda, const double & daughter_lambda, const double & grandD_lambda, const double & time){
      double term1 = daughter_start_amount*grandD_lambda/(grandD_lambda-daughter_lambda)*TMath::Exp(-daughter_lambda*time);
      double term2 = parent_intens*grandD_lambda*daughter_lambda/(grandD_lambda-parent_lambda)/(daughter_lambda-parent_lambda)*TMath::Exp(-parent_lambda*time);
      double term3 = parent_intens*grandD_lambda*daughter_lambda/(daughter_lambda-grandD_lambda)/(daughter_lambda-parent_lambda)*TMath::Exp(-daughter_lambda*time);
      double term4 = (parent_intens*daughter_lambda*grandD_lambda/(grandD_lambda-parent_lambda)/(grandD_lambda-daughter_lambda)+daughter_start_amount*grandD_lambda/(daughter_lambda-grandD_lambda))*TMath::Exp(-grandD_lambda*time);
         
      return term1 + term2 + term3 + term4;
}

/////////// Background function
Double_t Background(Double_t *x, Double_t *par){
   return par[2];
}

///////////// parent (34Mg) components
Double_t Parent_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda    = TMath::Log(2.0)/parent_hl;
   
   if(t>=beam_on_time && t<beam_off_time){
      return Implantation(parent_intens, parent_lambda, t-beam_on_time);
   }
   if(t>=beam_off_time){
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);  //parent_start_amount is like the A_0
      return parent_start_amount*TMath::Exp(-parent_lambda*(t-beam_off_time));  // simple exponential decay of the parent after the beam is off
   }
   return 0.0;
}

//////////// Daughter (34Al) Components
Double_t Daughter_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter_hl   = par[5];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter_lambda = TMath::Log(2.0)/daughter_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of daughter from decaying of parent during implantation
      double daughter_component = 0.0;
      daughter_component   += GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, t-beam_on_time);
      daughter_component   *= (1-pbn-pb2n);
      return daughter_component;
   }

   if(t>=beam_off_time){
      double daughter_component = 0.0;
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);  // Does the numerical term come from pn probablity/
      double daughter_start_amount = GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, beam_off_time-beam_on_time);
      
      daughter_component += daughter_start_amount*TMath::Exp(-daughter_lambda*(t-beam_off_time)); // Decay of daughter present after beam on
      daughter_component += GrowInDecayDaughter(parent_start_amount, parent_lambda, daughter_lambda, t-beam_off_time); //Grow-in and decay of daughter from parent amount after beam off
      daughter_component *= (1-pbn-pb2n);
      return daughter_component;
   }
   return 0.0;
}

////////////////  Granddaughter (34Si) Components
Double_t GrandD_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter_hl   = par[5];
   double grandD_hl   	= par[6];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter_lambda = TMath::Log(2.0)/daughter_hl;
   double grandD_lambda = TMath::Log(2.0)/grandD_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of daughter from decaying of parent during implantation
      double grandD_component = 0.0;
      grandD_component += GrowInImplantationGrandD(parent_intens, parent_lambda, daughter_lambda, grandD_lambda, t-beam_on_time);
      grandD_component *= (1-pbn-pb2n)*(1-dbn-db2n);
      return grandD_component;
   }

   if(t>=beam_off_time){
      double grandD_component = 0.0;
      double daughter_start_amount = (1-pbn-pb2n)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, beam_off_time-beam_on_time);
      double grandD_start_amount = (1-pbn-pb2n)*(1-dbn-db2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      grandD_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      grandD_component += (1-dbn-db2n)*GrowInDecayGrandD(parent_intens, daughter_start_amount, parent_lambda, daughter_lambda, grandD_lambda, t-beam_off_time);
      return grandD_component;
   }
   return 0.0;
}

//////////// Beta-n Daughter (33Al) Components
//This is the same as Daughter_Component with a different value for daughter_hl and scaling by pbn rather than 1-pbn-pb2n
Double_t BetaN_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter_hl   = par[7];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter_lambda = TMath::Log(2.0)/daughter_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of daughter from decaying of parent during implantation
      double daughter_component = 0.0;
      daughter_component   += GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, t-beam_on_time);
      daughter_component   *= pbn;
      return daughter_component;
   }

   if(t>=beam_off_time){
      double daughter_component = 0.0;
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);  // Does the numerical term come from pn probablity/
      double daughter_start_amount = GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, beam_off_time-beam_on_time);
      
      daughter_component += daughter_start_amount*TMath::Exp(-daughter_lambda*(t-beam_off_time)); // Decay of daughter present after beam on
      daughter_component += GrowInDecayDaughter(parent_start_amount, parent_lambda, daughter_lambda, t-beam_off_time); //Grow-in and decay of daughter from parent amount after beam off
      daughter_component *= pbn;
      return daughter_component;
   }
   return 0.0;
}

//////////// Beta-2n Daughter (32Al) Components
//This is the same as Daughter_Component with a different value for daughter_hl and scaling by pb2n rather than 1-pbn-pb2n
Double_t Beta2N_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter_hl   = par[9];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter_lambda = TMath::Log(2.0)/daughter_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of daughter from decaying of parent during implantation
      double daughter_component = 0.0;
      daughter_component   += GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, t-beam_on_time);
      daughter_component   *= pb2n;
      return daughter_component;
   }

   if(t>=beam_off_time){
      double daughter_component = 0.0;
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);  // Does the numerical term come from pn probablity/
      double daughter_start_amount = GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, beam_off_time-beam_on_time);
      
      daughter_component += daughter_start_amount*TMath::Exp(-daughter_lambda*(t-beam_off_time)); // Decay of daughter present after beam on
      daughter_component += GrowInDecayDaughter(parent_start_amount, parent_lambda, daughter_lambda, t-beam_off_time); //Grow-in and decay of daughter from parent amount after beam off
      daughter_component *= pb2n;
      return daughter_component;
   }
   return 0.0;
}

////////////////  two beta N (33Si)  Components.
//This is the same as the grandD_component but with some additions (multiple paths can bring you to the same place)
Double_t TwoBetaN_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter1_hl   = par[5];
   double grandD_hl   	= par[8];
   double daughter2_hl   = par[7];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter1_lambda = TMath::Log(2.0)/daughter1_hl;
   double daughter2_lambda = TMath::Log(2.0)/daughter2_hl;
   double grandD_lambda = TMath::Log(2.0)/grandD_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of two beta-n daughter from decaying of parent during implantation
      double twobetaN_component = 0.0;
      twobetaN_component += (1-pbn-pb2n)*dbn*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter1_lambda, grandD_lambda, t-beam_on_time);
      twobetaN_component += (pbn)*(1-dnbn-dnb2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter2_lambda, grandD_lambda, t-beam_on_time);
      return twobetaN_component;
   }

   if(t>=beam_off_time){
      double twobetaN_component = 0.0;
      double daughter1_start_amount = (1-pbn-pb2n)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter1_lambda, beam_off_time-beam_on_time);
      double grandD_start_amount = (1-pbn-pb2n)*(dbn)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter1_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      twobetaN_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      twobetaN_component += dbn*GrowInDecayGrandD(parent_intens, daughter1_start_amount, parent_lambda, daughter1_lambda, grandD_lambda, t-beam_off_time);

      double daughter2_start_amount = (pbn)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter2_lambda, beam_off_time-beam_on_time);
      grandD_start_amount = (pbn)*(1-dnbn-dnb2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter2_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      twobetaN_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      twobetaN_component += (1-dnbn-dnb2n)*GrowInDecayGrandD(parent_intens, daughter2_start_amount, parent_lambda, daughter2_lambda, grandD_lambda, t-beam_off_time);

      return twobetaN_component;
   }
   return 0.0;
}

////////////////  two beta N (33Si)  Components.
//This is the same as the grandD_component but with some additions (multiple paths can bring you to the same place)
Double_t TwoBeta2N_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter1_hl   = par[5];
   double grandD_hl   	= par[10];
   double daughter2_hl   = par[7];
   double daughter3_hl   = par[9];

   double t = x[0];
   //t = x[0]/(1-x[0]*corr);

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter1_lambda = TMath::Log(2.0)/daughter1_hl;
   double daughter2_lambda = TMath::Log(2.0)/daughter2_hl;
   double daughter3_lambda = TMath::Log(2.0)/daughter3_hl;
   double grandD_lambda = TMath::Log(2.0)/grandD_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins of two beta-n daughter from decaying of parent during implantation
      double twobetaN_component = 0.0;
      twobetaN_component += (1-pbn-pb2n)*db2n*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter1_lambda, grandD_lambda, t-beam_on_time);
      twobetaN_component += (pbn)*(dnbn)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter2_lambda, grandD_lambda, t-beam_on_time);
      twobetaN_component += (pb2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter3_lambda, grandD_lambda, t-beam_on_time);
      return twobetaN_component;
   }

   if(t>=beam_off_time){
      double twobetaN_component = 0.0;
      double daughter1_start_amount = (1-pbn-pb2n)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter1_lambda, beam_off_time-beam_on_time);
      double grandD_start_amount = (1-pbn-pb2n)*(db2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter1_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      twobetaN_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      twobetaN_component += (db2n)*GrowInDecayGrandD(parent_intens, daughter1_start_amount, parent_lambda, daughter1_lambda, grandD_lambda, t-beam_off_time);

      double daughter2_start_amount = (pbn)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter2_lambda, beam_off_time-beam_on_time);
      grandD_start_amount = (pbn)*(dnbn)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter2_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      twobetaN_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      twobetaN_component += (dnbn)*GrowInDecayGrandD(parent_intens, daughter2_start_amount, parent_lambda, daughter2_lambda, grandD_lambda, t-beam_off_time);

      double daughter3_start_amount = (pb2n)*GrowInImplantationDaughter(parent_intens, parent_lambda, daughter3_lambda, beam_off_time-beam_on_time);
      grandD_start_amount = (pb2n)*GrowInImplantationGrandD(parent_intens, parent_lambda, daughter3_lambda, grandD_lambda, beam_off_time-beam_on_time);
      
      twobetaN_component += grandD_start_amount*TMath::Exp(-grandD_lambda*(t-beam_off_time));
      twobetaN_component += GrowInDecayGrandD(parent_intens, daughter3_start_amount, parent_lambda, daughter3_lambda, grandD_lambda, t-beam_off_time);

      return twobetaN_component;
   }
   return 0.0;
} 

Double_t FitFunction(Double_t *x, Double_t *par){
   return Background(x,par)
          + Parent_Component(x,par)
          + Daughter_Component(x,par)
	      + GrandD_Component(x,par)
	      + BetaN_Component(x,par)
	      + Beta2N_Component(x,par)
	      + TwoBetaN_Component(x,par)
	      + TwoBeta2N_Component(x,par);
}

void MakeResidualPlot(TH1* hist_to_fit, TF1* fit_fn){
   gStyle->SetOptStat(0);
   auto resid_canvas = new TCanvas("resid_canvas", "Residuals of Fit");
	int size = hist_to_fit->GetNbinsX();
	double x[size], y[size];

	for(int i=0; i<size; ++i) {
		if(i<151) {
			x[i] = 0;
			y[i] = 0;
		}
		else {
			x[i] = hist_to_fit->GetXaxis()->GetBinCenter(i);
			y[i] = (hist_to_fit->GetBinContent(i) - fit_fn->Eval(x[i]))/TMath::Sqrt(hist_to_fit->GetBinContent(i));
		}
	}

	TGraph *tg = new TGraph(size,x,y);

	tg->Draw("A*");
	//TRatioPlot *rp1 = new TRatioPlot(hist_to_fit,"",fit_result.Get());
   /*rp1->Draw();
   rp1->GetLowerRefYaxis()->SetTitle("ratio");
   rp1->GetUpperRefYaxis()->SetTitle("entries");*/
	//cout <<"Or here?" << endl;
   //resid_canvas->Update();
	//cout<<"and finally?"<<endl;
}

TF1* DrawComponent(TF1* fit, int color, Double_t(*fit_func)(Double_t *, Double_t*),const char *name){
   double range_low, range_high;
   static int line_style = 2;
   fit->GetRange(range_low,range_high);
   TF1 *draw_func = new TF1(name,fit_func,range_low,range_high,11);
   for(int i = 0; i < fit->GetNpar();++i){
      draw_func->SetParameter(i,fit->GetParameter(i));
   }
   draw_func->SetLineStyle(line_style++);
   draw_func->SetNpx();
   draw_func->SetLineColor(color);
   draw_func->Print();
   draw_func->Draw("same");

   return draw_func;
}

Double_t IntegrateComponent(TF1* fit, Double_t(*fit_func)(Double_t *, Double_t *)){
   double range_low, range_high;
   fit->GetRange(range_low,range_high);
range_low = 2.0;	// what is this range?
range_high = 5.0;	// what is this range?
   TF1 *integral_func = new TF1("temp",fit_func,range_low,range_high,11);	
   for(int i = 0; i < fit->GetNpar(); ++i){
      integral_func->SetParameter(i,fit->GetParameter(i));
   }
   
   return integral_func->Integral(range_low,range_high);

}

Double_t IntegrateComponentErr(TF1* fit, TFitResultPtr fit_result, const Int_t &first_par,const Int_t &last_par){
   double range_low, range_high;
   fit->GetRange(range_low,range_high);
range_low = 2.0;
range_high = 5.0;
   TF1 *integral_func = static_cast<TF1*>(fit->Clone());
   integral_func->Print();
   TMatrixDSym copy_result = fit_result->GetCovarianceMatrix();
   for(int i = 0; i < fit->GetNpar(); ++i){
      integral_func->SetParameter(i,fit->GetParameter(i));
      //zero diagonals of covariance matrix for non par terms
      //if(((i < first_par) || (i > last_par)) && ((i == 2) || (i==4) || (i==6))) {
      if((i < first_par) || (i > last_par)) {
         copy_result(i,i) = 0.0;
         integral_func->SetParameter(i,0.0);
         //std::cout << integral_func->GetParameter(i) <<std::endl;
      }
      else{
         //std::cout << integral_func->GetParameter(i) <<std::endl;
      }
   }
  // copy_result.Print();
   return fit->IntegralError(range_low,range_high,integral_func->GetParameters(),copy_result.GetMatrixArray());

}

double RateFromIntegral(const double &integral, const double &half_life, const double &beam_on_t, const double &beam_off_t, const double &stop_t){

   double decay_rate = TMath::Log(2.)/half_life;
   double from_beam_on = (beam_off_t - beam_on_t) + (TMath::Exp(-decay_rate*beam_off_t) - TMath::Exp(-decay_rate*beam_on_t))/decay_rate;
   double from_beam_off = (1.-TMath::Exp(-decay_rate*(beam_off_t - beam_on_t)))*(TMath::Exp(-decay_rate*stop_t)-TMath::Exp(-decay_rate*beam_off_t))/decay_rate;

   return integral/(from_beam_on - from_beam_off);

}

double DecayFitwBetaN(TH1* raw_data){

	//Create a new histogram and apply dead time corrections
	TH1D* hist_to_fit = new TH1D("hist_to_fit", "Dead time corrected data", 500,0,5.);
	for(int i=1; i<=500; i++) {
		hist_to_fit->SetBinContent(i, raw_data->GetBinContent(i)/(1.-raw_data->GetBinContent(i)*corr));
	}

   double bin_width = hist_to_fit->GetBinWidth(1);

   int nBins = hist_to_fit->GetNbinsX();
   double t_bg;
   double t_on;
   double t_off;
   double t_stop;

   //Cycle used from run 15515 and onwards
   t_bg = 1.5;
   t_on = 2.0;
   t_off = 3.5;
   //t_stop = 3.5;
   t_stop = 5.0;

   double parent_hl   = 0.02;    // 34Mg, in sec
   double daughter_hl = 0.0563;    // 34Al
   double grandD_hl   = 2.77;      // 34Si
   double betaN_hl    = 0.0417;    // 33Al
   double twobetaN_hl = 6.11;      // 33Si
   double beta2N_hl    = 0.0319;    // 32Al
   double twobeta2N_hl = 153*365; // 32Si

   double bg_bins = TMath::Floor(t_on/bin_width) - TMath::Floor(t_bg/bin_width) - 1;
   double bg_level = hist_to_fit->Integral(TMath::Floor(t_bg/bin_width)+1,TMath::Floor(t_on/bin_width)-1);
	//bg_level/=(1.-corr); //Apply DT correction

   //spectrum is binned in s
   TF1 * cycle_fit_function = new TF1("Total Fit",FitFunction,t_bg,t_stop,11); 
   cycle_fit_function->SetParName(0,"t_on");
   cycle_fit_function->SetParName(1,"t_off");
   cycle_fit_function->SetParName(2,"bg_intensity");
   cycle_fit_function->SetParName(3,"parent_hl");
   cycle_fit_function->SetParName(4,"parent_inten");
   cycle_fit_function->SetParName(5,"daughter_hl");
   cycle_fit_function->SetParName(6,"grandD_hl");
   cycle_fit_function->SetParName(7,"betaN_hl");
   cycle_fit_function->SetParName(8,"twobetaN_hl");
   cycle_fit_function->SetParName(9,"beta2N_hl");
   cycle_fit_function->SetParName(10,"twobeta2N_hl");

   //Parameter guesses

   cycle_fit_function->FixParameter(0,t_on);
   cycle_fit_function->FixParameter(1,t_off);
   cycle_fit_function->FixParameter(2,bg_level/bg_bins); //simulated bg rate*cycles*bin time
   cycle_fit_function->FixParameter(3,parent_hl);
   //cycle_fit_function->FixParameter(4,3e4);
   cycle_fit_function->FixParameter(5,daughter_hl);
   cycle_fit_function->FixParameter(6,grandD_hl);
   cycle_fit_function->FixParameter(7,betaN_hl);
   cycle_fit_function->FixParameter(8,twobetaN_hl);
   cycle_fit_function->FixParameter(9,beta2N_hl);
   cycle_fit_function->FixParameter(10,twobeta2N_hl);


   TFitResultPtr fit_result = hist_to_fit->Fit(cycle_fit_function,"N0RS");//Q is quiet mode, no output of variables (used when fitting each cycle!)
   hist_to_fit->Draw();
   std::cout << "Chi2/NDF = " << fit_result->Chi2()/fit_result->Ndf() << std::endl;

   cycle_fit_function->SetLineColor(kRed);
   cycle_fit_function->SetNpx(5000);
   cycle_fit_function->Draw("same");
   TF1* parent_func = DrawComponent(cycle_fit_function, kBlue, Parent_Component, "^{34}Mg");
   TF1* daughter_func = DrawComponent(cycle_fit_function, kBlack, Daughter_Component, "^{34}Al Daughter");
   TF1* daughtern_func = DrawComponent(cycle_fit_function, kRed, BetaN_Component, "^{33}Al Daughter");
   TF1* daughter2n_func = DrawComponent(cycle_fit_function, kGreen, Beta2N_Component, "^{32}Al Daughter");
   TF1* grandDn_func = DrawComponent(cycle_fit_function, kGreen, TwoBetaN_Component, "^{33}Si Granddaughter");
   TF1* grandD2n_func = DrawComponent(cycle_fit_function, kRed, TwoBeta2N_Component, "^{32}Si Granddaughter");
   TF1* grandD_func = DrawComponent(cycle_fit_function, kBlack, GrandD_Component, "^{34}Si Granddaughter");

   TF1* bg_func = DrawComponent(cycle_fit_function, kGreen, Background, "Background");
   TLegend * leg = gPad->BuildLegend();
   leg->SetBorderSize(1);

   //Now we want to integrate to get the total number of betas.
   Double_t total_integral_fit  = IntegrateComponent(cycle_fit_function,FitFunction)/bin_width;
   Double_t background_integral_fit  = IntegrateComponent(cycle_fit_function,Background)/bin_width;
   Double_t parent_integral_fit      = IntegrateComponent(cycle_fit_function,Parent_Component)/bin_width;
   Double_t daughter_integral_fit    = IntegrateComponent(cycle_fit_function,Daughter_Component)/bin_width;
   Double_t grandD_integral_fit      = IntegrateComponent(cycle_fit_function,GrandD_Component)/bin_width;
   Double_t betaN_integral_fit       = IntegrateComponent(cycle_fit_function,BetaN_Component)/bin_width;
   Double_t twobetaN_integral_fit    = IntegrateComponent(cycle_fit_function,TwoBetaN_Component)/bin_width;

   Double_t total_integral_fit_err       = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t background_integral_fit_err  = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t parent_integral_fit_err      = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t daughter_integral_fit_err    = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t grandD_integral_fit_err      = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t betaN_integral_fit_err       = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;
   Double_t twobetaN_integral_fit_err    = IntegrateComponentErr(cycle_fit_function,fit_result,0,7)/bin_width;

   std::cout.precision(10);
   std::cout << "Integrals" << std::endl;
   std::cout << "*************************************************" << std::endl;
   std::cout << "total:    " << total_integral_fit << " +/- " << total_integral_fit_err << std::endl;
   std::cout << "background: " << background_integral_fit << " +/- " << background_integral_fit_err	 << std::endl;
   std::cout << "parent: "     << parent_integral_fit     << " +/- " << parent_integral_fit_err     	 << std::endl;
   std::cout << "daughter: "   << daughter_integral_fit   << " +/- " << daughter_integral_fit_err   	 << std::endl;
   std::cout << "granddaughter: "   << grandD_integral_fit   << " +/- " << grandD_integral_fit_err   	 << std::endl;
   std::cout << "beta-n: "   << betaN_integral_fit   << " +/- " << betaN_integral_fit_err   		 << std::endl;
   std::cout << "twobeta-n: "   << twobetaN_integral_fit   << " +/- " << twobetaN_integral_fit_err   		 << std::endl;
   std::cout << "*************************************************" << std::endl;


   //std::cout << "background: "  << background_integral_fit  << " +/- " << background_integral_fit_err << std::endl;
   //MakeResidualPlot(hist_to_fit,cycle_fit_function);

  // new TCanvas;
  // cycle_fit_function->Draw();
  
  return fit_result->Chi2()/fit_result->Ndf();
}
