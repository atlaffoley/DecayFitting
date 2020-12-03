
Double_t Implantation(const double & intens, const double & lambda, const double & time){
   return intens*(1. - TMath::Exp(-lambda*time));
}

Double_t GrowInImplantationDaughter(const double & parent_intens, const double & parent_lambda, const double & daughter_lambda, const double & time){
   double term1 = parent_intens*daughter_lambda/(daughter_lambda-parent_lambda)*(TMath::Exp(-daughter_lambda*time)-TMath::Exp(-parent_lambda*time));
   double term2 = parent_intens*(1.-TMath::Exp(-daughter_lambda*time));

   return term1 + term2;
}

Double_t GrowInDecayDaughter(const double & parent_initial_amount, const double & parent_lambda, const double & daughter_lambda, const double & time){
         return parent_initial_amount*daughter_lambda/(daughter_lambda-parent_lambda)*(TMath::Exp(-parent_lambda*time)-TMath::Exp(-daughter_lambda*time));
}

Double_t Background(Double_t *x, Double_t *par){
   return par[2];
}

Double_t Parent_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];

   double t = x[0];

   double parent_lambda    = TMath::Log(2.0)/parent_hl;
   
   if(t>=beam_on_time && t<beam_off_time){
      return Implantation(parent_intens, parent_lambda, t-beam_on_time);
   }
   if(t>=beam_off_time){
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);
      return parent_start_amount*TMath::Exp(-parent_lambda*(t-beam_off_time));
   }
   return 0.0;
}

Double_t Daughter_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double parent_hl     = par[3];
   double parent_intens = par[4];
   double daughter_hl   = par[5];

   double t = x[0];

   double parent_lambda   = TMath::Log(2.0)/parent_hl;
   double daughter_lambda = TMath::Log(2.0)/daughter_hl;

   if(t>=beam_on_time && t<beam_off_time){
      //Now grow-ins from decaying during implantation
      double daughter_component = 0.0;
      daughter_component   += GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, t-beam_on_time);
      return daughter_component;
   }

   if(t>=beam_off_time){
      double daughter_component = 0.0;
      double parent_start_amount = Implantation(parent_intens,parent_lambda,beam_off_time - beam_on_time);
      double daughter_start_amount = GrowInImplantationDaughter(parent_intens, parent_lambda, daughter_lambda, beam_off_time-beam_on_time);
      
      daughter_component += daughter_start_amount*TMath::Exp(-daughter_lambda*(t-beam_off_time));
      daughter_component += GrowInDecayDaughter(parent_start_amount, parent_lambda, daughter_lambda, t-beam_off_time);
      //daughter_component += parent_start_amount*daughter_lambda/(daughter_lambda-parent_lambda)*(TMath::Exp(-parent_lambda*(t-beam_off_time)) - TMath::Exp(-daughter_lambda*(t-beam_off_time)));
      return daughter_component;
   }

   return 0.0;

}

Double_t Contam1_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double contam_hl     = par[5];
   double contam_intens = par[6];
    
   double t = x[0];

   double contam_lambda   = TMath::Log(2.0)/contam_hl;
   
   if(t>=beam_on_time && t<beam_off_time){
      return Implantation(contam_intens, contam_lambda, t-beam_on_time);
   }

   if(t>=beam_off_time){
      double contam_start_amount = Implantation(contam_intens,contam_lambda,beam_off_time - beam_on_time);
      return contam_start_amount*TMath::Exp(-contam_lambda*(t-beam_off_time));
   }

   return 0.0;
}

Double_t Contam2_Component(Double_t *x, Double_t *par){
   double beam_on_time  = par[0];
   double beam_off_time = par[1];
   double contam_hl     = par[7];
   double contam_intens = par[8];
    
   double t = x[0];

   double contam_lambda   = TMath::Log(2.0)/contam_hl;
   
   if(t>=beam_on_time && t<beam_off_time){
      return Implantation(contam_intens, contam_lambda, t-beam_on_time);
   }

   if(t>=beam_off_time){
      double contam_start_amount = Implantation(contam_intens,contam_lambda,beam_off_time - beam_on_time);
      return contam_start_amount*TMath::Exp(-contam_lambda*(t-beam_off_time));
   }

   return 0.0;
}

Double_t FitFunction(Double_t *x, Double_t *par){

   return Background(x,par)
          + Parent_Component(x,par)
          + Daughter_Component(x,par)
          + Contam1_Component(x,par)
          + Contam2_Component(x,par);
}

void MakeResidualPlot(TH1* hist_to_fit, TFitResultPtr fit_result){
   gStyle->SetOptStat(0);
   auto resid_canvas = new TCanvas("resid_canvas", "Residuals of Fit");
   TRatioPlot *rp1 = new TRatioPlot(hist_to_fit,"",fit_result.Get());
   /*rp1->Draw();
   rp1->GetLowerRefYaxis()->SetTitle("ratio");
   rp1->GetUpperRefYaxis()->SetTitle("entries");*/
   resid_canvas->Update();
}

TF1* DrawComponent(TF1* fit, int color, Double_t(*fit_func)(Double_t *, Double_t*),const char *name){
   double range_low, range_high;
   static int line_style = 2;
   fit->GetRange(range_low,range_high);
   TF1 *draw_func = new TF1(name,fit_func,range_low,range_high,9);
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
   TF1 *integral_func = new TF1("temp",fit_func,range_low,range_high,9);
   for(int i = 0; i < fit->GetNpar(); ++i){
      integral_func->SetParameter(i,fit->GetParameter(i));
   }
   
   return integral_func->Integral(range_low,range_high);

}

Double_t IntegrateComponentErr(TF1* fit, TFitResultPtr fit_result, const Int_t &first_par,const Int_t &last_par){
   double range_low, range_high;
   fit->GetRange(range_low,range_high);
   TF1 *integral_func = static_cast<TF1*>(fit->Clone());
   integral_func->Print();
   TMatrixDSym copy_result = fit_result->GetCovarianceMatrix();
   for(int i = 0; i < fit->GetNpar(); ++i){
      integral_func->SetParameter(i,fit->GetParameter(i));
      //zero diagonals of covariance matrix for non par terms
      if(((i < first_par) || (i > last_par)) && ((i == 2) || (i==4) || (i==6))) {
         copy_result(i,i) = 0.0;
         integral_func->SetParameter(i,0.0);
      //   std::cout << integral_func->GetParameter(i) <<std::endl;
      }
      else{
      //   std::cout << integral_func->GetParameter(i) <<std::endl;
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

void DecayFit(TH1* hist_to_fit){

//The first thing we need to do is create and equation.
   //Ar -> 1 rate + 1 half-life
   //Cl (daughter) -> 1 half-life
   //Cl (contam) -> 1 rate
   //Contam1 (34mCl) -> 1 rate + 1 half-life
   //Background -> 1 rate
   //beam on time -> 1 par
   //beam off time -> 1 par
   //9 total parameters.
   
   double bin_width = hist_to_fit->GetBinWidth(10);

   double t_bg = 0.0;
   double t_on = 1.0;
   double t_off = 21.0;
   double t_stop = 24.0;

   double parent_hl = 0.8438;
   double daughter_hl = 1.52655;
   double contam2_hl = 1920.;

   double parent_integral = 2984499.455;

   //spectrum is binned in s
   TF1 * cycle_fit_function = new TF1("Total Fit",FitFunction,t_bg,t_stop,9); //Changing 9 here means also changing it above in DrawComponent and IntegrateComponent
   cycle_fit_function->SetParName(0,"t_on");
   cycle_fit_function->SetParName(1,"t_off");
   cycle_fit_function->SetParName(2,"bg_intensity");
   cycle_fit_function->SetParName(3,"parent_hl");
   cycle_fit_function->SetParName(4,"parent_inten");
   cycle_fit_function->SetParName(5,"daughter_hl");
   cycle_fit_function->SetParName(6,"daughter_inten"); //Intensity of daughter as a CONTAMINANT
   cycle_fit_function->SetParName(7,"contam2_hl");
   cycle_fit_function->SetParName(8,"contam2_inten");

   //Parameter guesses
   cycle_fit_function->SetParameter(0,t_on);
   cycle_fit_function->SetParameter(1,t_off);
   cycle_fit_function->SetParameter(2,36.);
   cycle_fit_function->SetParameter(3,parent_hl);
   cycle_fit_function->SetParameter(4,2088000.);
   cycle_fit_function->SetParameter(5,daughter_hl);
   cycle_fit_function->SetParameter(6,0.);
   cycle_fit_function->SetParameter(7,contam2_hl);
   cycle_fit_function->SetParameter(8,0.);

   cycle_fit_function->FixParameter(0,t_on);
   cycle_fit_function->FixParameter(1,t_off);
//   cycle_fit_function->FixParameter(2,36.); //simulated bg rate*cycles*bin time
   cycle_fit_function->FixParameter(3,parent_hl);
   cycle_fit_function->FixParameter(5,daughter_hl);
   cycle_fit_function->FixParameter(7,contam2_hl);
   cycle_fit_function->FixParameter(8,0.);

   TFitResultPtr fit_result = hist_to_fit->Fit(cycle_fit_function,"N0LIRS");
   hist_to_fit->Draw();
   std::cout << "Chi2/NDF = " << fit_result->Chi2()/fit_result->Ndf() << std::endl;

   cycle_fit_function->SetLineColor(kRed);
   cycle_fit_function->SetNpx(5000);
   cycle_fit_function->Draw("SamE");
   TF1* parent_func = DrawComponent(cycle_fit_function, kMagenta, Parent_Component, "^{34}Ar");
   TF1* daughter_func = DrawComponent(cycle_fit_function, kBlack, Daughter_Component, "^{34}Cl");
   TF1* contam_func = DrawComponent(cycle_fit_function, kBlue, Contam1_Component, "^{34}Cl (contamination)");
   TF1* contam2_func = DrawComponent(cycle_fit_function, kViolet, Contam2_Component, "^{34m}Cl (contamination)");
   TF1* bg_func = DrawComponent(cycle_fit_function, kGreen, Background, "Background");
   TLegend * leg = gPad->BuildLegend();
   leg->SetBorderSize(1);

   //Now we want to integrate to get the total number of betas.
   Double_t background_integral_fit  = IntegrateComponent(cycle_fit_function,Background)/bin_width;
   Double_t parent_integral_fit      = IntegrateComponent(cycle_fit_function,Parent_Component)/bin_width;
   Double_t daughter_integral_fit    = IntegrateComponent(cycle_fit_function,Daughter_Component)/bin_width;
   Double_t daughter_contam_integral_fit    = IntegrateComponent(cycle_fit_function,Contam1_Component)/bin_width;

   Double_t background_integral_fit_err  = IntegrateComponentErr(cycle_fit_function,fit_result,1,1)/bin_width;
   Double_t parent_integral_fit_err      = IntegrateComponentErr(cycle_fit_function,fit_result,3,4)/bin_width;
   Double_t daughter_integral_fit_err    = IntegrateComponentErr(cycle_fit_function,fit_result,3,5)/bin_width;
   Double_t daughter_contam_integral_fit_err    = IntegrateComponentErr(cycle_fit_function,fit_result,5,6)/bin_width;

   std::cout << "Integrals" << std::endl;
   std::cout << "*************************************************" << std::endl;
   std::cout << "background: " << background_integral_fit << " +/- " << background_integral_fit_err << std::endl;
   std::cout << "parent: "     << parent_integral_fit     << " +/- " << parent_integral_fit_err     << std::endl;
   std::cout << "daughter: "   << daughter_integral_fit   << " +/- " << daughter_integral_fit_err   << std::endl;
   std::cout << "daughter_contam: "   << daughter_contam_integral_fit   << " +/- " << daughter_contam_integral_fit_err   << std::endl;
   std::cout << "daughter_total: "    << daughter_integral_fit + daughter_contam_integral_fit << " +/- " << TMath::Sqrt(daughter_integral_fit_err*daughter_integral_fit_err + daughter_contam_integral_fit_err*daughter_contam_integral_fit_err)         << std::endl;
   std::cout << "*************************************************" << std::endl;


   //std::cout << "background: "  << background_integral_fit  << " +/- " << background_integral_fit_err << std::endl;
   //MakeResidualPlot(hist_to_fit,fit_result);

  // new TCanvas;
  // cycle_fit_function->Draw();


}


