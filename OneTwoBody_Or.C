void OneTwoBody_Or(){

	const int N = 1000;         //Originally set to 10000 by Or. Currently running with 90000
	Double_t P1[N+10][3];       //Momentum vector for 1 nucleon
	Double_t P[N+10];
	Double_t Prel[3], Pcm[3];
	Double_t Pcom[3];           // Center of Mass Momentum
	Double_t tmp;
	Double_t tmp1;

	Double_t kfermi = 250./197.3; //MeV/(hbar*c)

	// Index in Pcom[i]:  0     ,     1     ,    2     ,     3    ,    4   ,     5
	Double_t Q[6] = {0.05*kfermi,0.25*kfermi,0.5*kfermi,1.0*kfermi,2*kfermi,3.5*kfermi};

	TH1F* hist_1body         = new TH1F("","",100,0,5);
	TH1F** hist_2body        = new TH1F*[7];
	TH1F** hist_2body_MF     = new TH1F*[7];
	TH1F** hist_2body_SRC    = new TH1F*[7];
	TH1F** hist_2body_SRC_MF = new TH1F*[7];
	TH1F** hist_2body_RealSRC= new TH1F*[7];
	TH1F** hist_Ratio        = new TH1F*[7];
	for(int i = 0 ; i < 7; i++){
		hist_2body[i]        = new TH1F("","",100,0,5);  hist_2body[i]   	->SetLineColor(1); // Red
		hist_2body_MF[i]     = new TH1F("","",100,0,5);  hist_2body_MF[i]    	->SetLineColor(2); // Red
		hist_2body_SRC[i]    = new TH1F("","",100,0,5);  hist_2body_SRC[i]   	->SetLineColor(4); // Blue
		hist_2body_SRC_MF[i] = new TH1F("","",100,0,5);  hist_2body_SRC_MF[i]	->SetLineColor(3); // Green
		hist_2body_RealSRC[i]= new TH1F("","",100,0,5);  hist_2body_RealSRC[i]	->SetLineColor(6); // Green
		hist_Ratio[i]        = new TH1F("","",100,0,5);
	}

	TH2F* h2_MF_MF    = new TH2F("MF_MF"   ,"MF_MF"   , 90, 0, 5, 90, 0, 7);
	TH2F* h2_SRC_MF   = new TH2F("SRC_MF"  ,"SRC_MF"  , 90, 0, 5, 90, 0, 7);
	TH2F* h2_SRC_SRC  = new TH2F("SRC_SRC" ,"SRC_SRC" , 90, 0, 5, 90, 0, 7);
	TH2F* h2_true_SRC = new TH2F("true_SRC","true_SRC", 90, 0, 5, 90, 0, 7);

	TF1* n_k     = new TF1(""," (x<(250/197.3))?0.8*3/((250/197.3)**3):0.5*0.2/(1./(250/197.3)-1./5)/(x**4)",0,5);
	TF1* n_k_k2  = new TF1("","((x<(250/197.3))?0.8*3/((250/197.3)**3):0.5*0.2/(1./(250/197.3)-1./5)/(x**4))*(x**2)",0,5);

	// ****************************************************************
	//           Initiate Monte Carlo sim to build database
	// ****************************************************************
	int counter1 = 0; int counter2=0;
	for(int i=0; i<N; i++){

		P[i] = n_k_k2 ->GetRandom();
		gRandom->Sphere(P1[i][0], P1[i][1], P1[i][2], P[i]);
		hist_1body->Fill( P[i] );

		if( P[i]>250/197.3 ){ //Above k fermi

			//Pcm[0] = gRandom->Gaus(0,140)/197.3;
			//Pcm[1] = gRandom->Gaus(0,140)/197.3;
			//Pcm[2] = gRandom->Gaus(0,140)/197.3;
			
			Pcm[0] = 0.;
                        Pcm[1] = 0.;
                        Pcm[2] = 0.;
			
			// This next piece constructs correlated pair
			P1[i+1][0] = Pcm[0]-P1[i][0];
			P1[i+1][1] = Pcm[1]-P1[i][1];
			P1[i+1][2] = Pcm[2]-P1[i][2];
			P[i+1] = sqrt(pow(P1[i+1][0],2)+pow(P1[i+1][1],2)+pow(P1[i+1][2],2));
			hist_1body->Fill( P[i+1] );


			Prel[0] = (P1[i][0]-P1[i+1][0])/2.;
			Prel[1] = (P1[i][1]-P1[i+1][1])/2.;
			Prel[2] = (P1[i][2]-P1[i+1][2])/2.;
			tmp = sqrt(pow(Prel[0],2)+pow(Prel[1],2)+pow(Prel[2],2));

			Pcom[0] = (P1[i][0]+P1[i+1][0]);
			Pcom[1] = (P1[i][1]+P1[i+1][1]);
			Pcom[2] = (P1[i][2]+P1[i+1][2]);
			tmp1 = sqrt(pow(Pcom[0],2)+pow(Pcom[1],2)+pow(Pcom[2],2));

			h2_true_SRC -> Fill(tmp , tmp1);

			for(int k=0; k<6; k++){
				if(tmp1 < Q[k]){
					hist_2body_RealSRC[k]->Fill( tmp );
				}
			}

			i++;
			counter1+=2;
		} //End if above k fermi

		else{counter2++;}

	} //End loop over i [0,N-1]
	hist_1body->Scale(0.5);


	cout << "Building BD Done\t(" << 100.*counter1/(counter1+counter2)<<"% SRC)" << endl;
	// ******************************************************

	for(int i=0; i<N; i++){
		for(int j=i+1; j<N; j++){
			Prel[0] = (P1[i][0]-P1[j][0])/2.;
			Prel[1] = (P1[i][1]-P1[j][1])/2.;
			Prel[2] = (P1[i][2]-P1[j][2])/2.;
			tmp = sqrt(pow(Prel[0],2)+pow(Prel[1],2)+pow(Prel[2],2));

			Pcom[0] = (P1[i][0]+P1[j][0]);
			Pcom[1] = (P1[i][1]+P1[j][1]);
			Pcom[2] = (P1[i][2]+P1[j][2]);
			tmp1 = sqrt(pow(Pcom[0],2)+pow(Pcom[1],2)+pow(Pcom[2],2));


			for(int k=0; k<6; k++){

				if(tmp1 < Q[k]){

					hist_2body[k]->Fill( tmp );
					if     (P[i]>250/197.3 && P[j]>250/197.3){hist_2body_SRC[k]   ->Fill( tmp );}
					else if(P[i]<250/197.3 && P[j]<250/197.3){hist_2body_MF[k]    ->Fill( tmp );}
					else                                     {hist_2body_SRC_MF[k]->Fill( tmp );}

				}
			}


			if     (P[i]>250/197.3 && P[j]>250/197.3){h2_SRC_SRC  	->Fill(tmp , tmp1);}
               		else if(P[i]<250/197.3 && P[j]<250/197.3){h2_MF_MF    	->Fill(tmp , tmp1);}
                      	else                                     {h2_SRC_MF	->Fill(tmp , tmp1);}




		}
		if(i%(N/100) == 0){cout << 100*i/N << endl;}
	}
	//    cout <<hist_2body_SRC->GetEntries()<<"\t"<<hist_2body_MF->GetEntries()<<"\t"<<hist_2body_SRC_MF->GetEntries()<<endl;

	// ******************************************************
	// Editing the histograms
	for(int k=0; k<6; k++){
		Pretty(hist_2body[k]         , 0, k);      
		Pretty(hist_2body_MF[k]      , 1, k);
		Pretty(hist_2body_SRC[k]     , 2, k);
		Pretty(hist_2body_SRC_MF[k]  , 3, k);
		Pretty(hist_2body_RealSRC[k] , 4, k);
		Pretty(hist_1body            , 5, k);
	}


	// ******************************************************
	// Drawing the histograms
	leg = new TLegend(0.6,0.6,0.89,0.89);
	leg->SetTextSize(0.05);
	leg->AddEntry(	hist_2body[0]		,"Total");
	leg->AddEntry(	hist_2body_MF[0]	,"MF - MF");
	leg->AddEntry(	hist_2body_SRC_MF[0] 	,"MF - SRC");
	leg->AddEntry(	hist_2body_SRC[0]	,"SRC - SRC");
	leg->AddEntry(	hist_2body_RealSRC[0]	,"True SRC");
	
	int ci = 925;
	color = new TColor(ci, 0, 0, 0, " ", 0);
	leg->SetLineColor(ci);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);

	TCanvas *can = new TCanvas("can","1",900,700);
	can->Divide(3,2);
	for(int k=0; k<6; k++){
		TPad *p1 = can->cd(k+1);
		p1->SetLogy();
		hist_2body[k]->Draw();
		hist_2body_MF[k]      -> Draw("same");
		hist_2body_SRC_MF[k]  -> Draw("same");
		hist_2body_SRC[k]     -> Draw("same");
		hist_2body_RealSRC[k] -> Draw("same");
		hist_1body            -> Draw("same");

		leg->Draw();
	}             

	TCanvas *can1 = new TCanvas("can1","1",1000,500);    
	can1->Divide(3,1);
	TPad *p2 = can1->cd(1);	p2->SetLogy();
	p2->SetBottomMargin(0.15);
	p2->SetRightMargin(0.03);
	hist_2body[0]->Draw();
	hist_2body_MF[0]      -> Draw("same");
	hist_2body_SRC_MF[0]  -> Draw("same");
	hist_2body_SRC[0]     -> Draw("same");
	hist_2body_RealSRC[0] -> Draw("same"); 


	TPad *p3 = can1->cd(2); p3->SetLogy();
	p3->SetBottomMargin(0.15);
	p3->SetRightMargin(0.03);
	hist_2body[3]->Draw();
	hist_2body_MF[3]      -> Draw("same");
	hist_2body_SRC_MF[3]  -> Draw("same");
	hist_2body_SRC[3]     -> Draw("same");
	hist_2body_RealSRC[3] -> Draw("same");


	leg->Draw();

	TPad *p4 = can1->cd(3); p4->SetLogy();
	p4->SetBottomMargin(0.15);
	p4->SetRightMargin(0.03);
	hist_2body[5]->Draw();
	hist_2body_MF[5]      -> Draw("same");
	hist_2body_SRC_MF[5]  -> Draw("same");
	hist_2body_SRC[5]     -> Draw("same");
	hist_2body_RealSRC[5] -> Draw("same");

	can1      -> Print("out.pdf" ,"pdf");

	
	Pretty2D(h2_MF_MF    , 1);
        Pretty2D(h2_SRC_SRC  , 2);
        Pretty2D(h2_SRC_MF   , 3);
	Pretty2D(h2_true_SRC , 0);

	TCanvas *can2 = new TCanvas("can2","3",1200,800);
	can2 -> Divide(2,2);
	can2->cd(1);	h2_SRC_SRC  -> Draw("colz");
	can2->cd(4);	h2_true_SRC -> Draw("colz");
	can2->cd(3);	h2_MF_MF    -> Draw("colz");
	can2->cd(2);	h2_SRC_MF   -> Draw("colz");
	/*
	leg1 = new TLegend(0.65,0.75,0.89,0.89);
        leg1->SetTextSize(0.04);
	leg1->SetLineColor(ci);
        leg1->SetLineStyle(1);
        leg1->SetLineWidth(1);
        leg1->AddEntry( h2_MF_MF     ,"MF - MF");
        leg1->AddEntry( h2_SRC_MF    ,"MF - SRC");
        leg1->AddEntry( h2_SRC_SRC   ,"SRC - SRC");
	leg1->AddEntry( h2_true_SRC  ,"true SRC");
	leg1->Draw();
	*/
}

// ****************************************************************************************
// Function to format plots
// ****************************************************************************************
void Pretty(TH1F *gP, int k, int opt){
	int color; int style; int linestyle;
	if( k == 0){ color = 1; style = 3004;}
	if( k == 1){ color = 3; style = 3004; linestyle = 2;}
	if( k == 2){ color = 2; style = 3005; linestyle = 7;}
	if( k == 3){ color = 4; style = 3007; linestyle = 9;}
	if( k == 4){ color = 7; style = 3006;}
	if( k == 5){ color = 6; style = 3005;}

	AddTitle(gP,opt);

	if( k == 0){gP -> SetFillColor(0); gP->SetLineWidth(2);}
	else{gP -> SetFillColor(color);}
	gStyle ->SetOptStat(0);
	gP -> SetFillStyle(style);
	//gP -> SetLineStyle(linestyle);
	gP -> SetLineColor(color);
	gP -> SetMarkerColor(color);
	gP -> GetXaxis()->SetTitle("q[fm^{-1}]");	
	gP -> GetXaxis()->SetRangeUser(0,5);
	gP -> GetXaxis()->SetTitleSize(0.06);
	gP -> GetYaxis()->SetLabelSize (0.06);
	gP -> GetXaxis()->SetLabelSize (0.07); 
	gP -> GetXaxis()->SetNdivisions(505);	
	gStyle->SetTitleSize(0.09,"t");
}

// ****************************************************************************************
// Function to add a title
// ****************************************************************************************
void AddTitle(TH1F *gP ,int opt){
	if(opt ==  0){gP -> SetTitle("Q = 0");}
	else if(opt == 1){gP -> SetTitle("Q < k_{f} /4");}
	else if(opt == 2){gP -> SetTitle("Q < k_{f} /2");}
	else if(opt == 3){gP -> SetTitle("Q < k_{f}"  );}
	else if(opt == 4){gP -> SetTitle("Q < 2k_{f}" );}
	else if(opt == 5){gP -> SetTitle("All Q"      );}
}

// ****************************************************************************************
// Function to format 2D plots
// ****************************************************************************************
void Pretty2D(TH2F *gP, int k){
        int color;
	
      	if( k == 0){ color = 1; gP ->SetTitle("True-SRC");} // Black 
        if( k == 1){ color = 3;	gP ->SetTitle("MF-MF"	);} // Green
        if( k == 2){ color = 2;	gP ->SetTitle("SRC-SRC"	);} // Red
        if( k == 3){ color = 4;	gP ->SetTitle("MF-SRC"	);} // Blue
        gP -> SetFillColor(color);
        gStyle ->SetOptStat(0);
	gStyle->SetPalette(54); 
	gP -> SetLineColor(color);
        gP -> SetMarkerColor(color);
        gP -> GetXaxis()->SetTitle("q[fm^{-1}]");
        gP -> GetYaxis()->SetTitle("Q[fm^{-1}]");
	gP -> GetYaxis()-> SetTitleOffset(1.4);
	gP -> GetYaxis()->SetTitleSize(0.055);
	gP -> GetXaxis()->SetTitleSize(0.055);
        gP -> GetYaxis()->SetLabelSize(0.055);
        gP -> GetXaxis()->SetLabelSize(0.055);
	gStyle->SetTitleSize(0.07,"t");
	gP -> GetXaxis()->SetTitleOffset(0.85);
	gP -> GetYaxis()->SetTitleOffset(0.60);
}














