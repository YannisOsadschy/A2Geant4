#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <limits>

struct SegPosition
{
    Double_t rMin;
    Double_t rMax;
    Double_t phiMin;
    Double_t phiMax;

};


//Global definitions
Double_t vdrift; //cm/ms //is now set automatically
Double_t qeslope = -4.26975e-05; //from output of Calibrate.C   //old value: -0.000208
Double_t qeintercept = 0.0013405; //from output of Calibrate.C	//old value: -0.89


void GetAnodeLayoutFromDataFile(std::string fileName, std::vector<Double_t>& radii, std::vector<Int_t>& nSegments)
{
	std::ifstream inputStream(fileName);

	if (!inputStream)
	{
		throw std::runtime_error(fileName+" does not exist");
	}
	std::string line;   
	Int_t ierr = 0;
	while (std::getline(inputStream, line))
	{
		if (line.empty() || line[0]=='#')
		{
			continue;
		}
		Int_t nSegment;
		Double_t radius;
		std::istringstream iss(line);
		std::string key;
		iss>>key;
		if (key=="Anode-Radii:")
		{
			radii.clear();
			while(iss>>radius)
			{
				radii.push_back(radius);
			}
			++ierr;
		}
		else if (key=="Anode-Segments:")
		{
			nSegments.clear();
			while(iss>>nSegment)
			{
				nSegments.push_back(nSegment);
			}
			++ierr;
		}	
	}
	if (ierr!=2)
	{
		throw std::runtime_error("format of anode layout data file is incorrect");
	}
}


SegPosition GetSegPosition(Int_t measuredSegID)
{
	std::vector<Int_t> nSegments;
    std::vector<Double_t> radii;
	//actual design used in the simulation
	GetAnodeLayoutFromDataFile("data/TPC.dat",radii,nSegments);
    Int_t segID=0;
    for (std::size_t i=0; i<radii.size()-1; ++i)
    {
        SegPosition segPostion;
        segPostion.rMin = radii[i];
        segPostion.rMax = radii[i+1];
        for (std::size_t j=0; j<nSegments[i]; ++j)
        {
            segPostion.phiMin = 360./nSegments[i]*j;
            segPostion.phiMax = 360./nSegments[i]*(j+1);
            ++segID;
            if (measuredSegID==segID)
            {
                return segPostion;
            }
        }
    }
    throw std::runtime_error("measuredSegID does not exist in Cathode design");
}

Double_t GetDistance(Double_t r1, Double_t phi1, Double_t r2, Double_t phi2 )
{
    return std::sqrt(r1*r1+r2*r2-2*r1*r2*cos((phi1-phi2)*3.1415926/180));
}

//Reconstruct recoil polar angle
Double_t GetTheta(Double_t delta_s,Double_t delta_t)
{
	Double_t theta = 0;
	if(delta_s!=0 && delta_t!=0)theta=90-atan((delta_t*vdrift*10)/delta_s)*(180/3.141);
	return theta;
}

//Reconstruct event position
Double_t GetZ(Double_t time){
	Double_t z_pos =-11.5 +vdrift *time;
	return z_pos;
}

//Reconstruct recoil energy
Double_t GetEnergy(Double_t charge){
	Double_t energy = qeintercept + (qeslope)*charge;
	return energy;
}

double GetMean(const std::vector<double>& list)
{
    double sum = 0;
    int n =list.size();
    for (std::size_t i=0; i<n; ++i)
    {
        sum += list[i]; 
    }
    return sum/n;
}

double GetStandardDeviation(const std::vector<double>& list)
{
    double sum = 0;
    int n =list.size();
    double mean = GetMean(list);
    for (std::size_t i=0; i<n; ++i)
    {
        sum += (list[i]-mean)*(list[i]-mean); 
    }
    return std::sqrt(sum/(n-1));
}

std::vector<std::size_t> GetArgSort(const std::vector<double>& list)
{
    int n =list.size();
    std::vector<std::size_t> index;
    for (std::size_t i=0; i<n; ++i)
    {
        index.push_back(i);
    }
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j){return list[i] < list[j];});
    return index;
}

void TransformCircleToCartresian(double& r, double& phi, double& x, double& y)
{
    x = r*std::cos(3.14159/180*phi);
    y = r*std::sin(3.14159/180*phi);
}



void GetDistanceAndTime1(Int_t ntpc, Float_t* ttpc, Float_t* qtpc, Int_t* itpc, Double_t& delta_s, Double_t& delta_t)
{
    Double_t timeFirstSeg = ttpc[0];
	SegPosition posAngularFirstSeg = GetSegPosition(itpc[0]);
	Double_t timeLastSeg = ttpc[ntpc-1];
	SegPosition posAngularLastSeg = GetSegPosition(itpc[ntpc-1]);

    Double_t rFirst = posAngularFirstSeg.rMin;
    Double_t phiFirst = (posAngularFirstSeg.phiMin+posAngularFirstSeg.phiMax)/2;
    Double_t rLast = posAngularLastSeg.rMin;
    Double_t phiLast = (posAngularLastSeg.phiMin+posAngularLastSeg.phiMax)/2;

    delta_s = GetDistance(rFirst, phiFirst, rLast, phiLast);
    delta_t = timeLastSeg-timeFirstSeg;
}


void GetDistanceAndTime3(Int_t ntpc, Float_t* ttpc, Float_t* qtpc, Int_t* itpc, Float_t charge, Double_t& delta_s, Double_t& delta_t, Bool_t& disregardLowChargePads, std::vector<Int_t>& disregardedPads, std::vector<Int_t>& selectedPads, Bool_t& useForAngular)
{
	//filter
	disregardLowChargePads = true;
	std::vector<Float_t> ttpc2, itpc2, qtpc2;
	Int_t ntpc2 = 0;
	Double_t maxCharge = 0;
   
	for (std::size_t i=0; i<ntpc; ++i)
	{
		if (std::abs(qtpc[i])>0.01*std::abs(charge))
		{
			ttpc2.push_back(ttpc[i]);
			qtpc2.push_back(qtpc[i]);
			itpc2.push_back(itpc[i]);
			++ntpc2;
		}
		else
		{
			disregardedPads.push_back(itpc[i]);
		}
		if (std::abs(qtpc[i]) > std::abs(maxCharge))
		{
			maxCharge = qtpc[i];
		}	
	}

	if (std::abs(maxCharge)>0.9*std::abs(charge) || ntpc2<2) //0.9
	{
		useForAngular = false;
	}
	if (useForAngular)
	{
	    Double_t timeFirstSeg = ttpc[0];
		SegPosition posAngularFirstSeg = GetSegPosition(itpc2[0]);
		selectedPads.push_back(itpc2[0]);
		Double_t timeLastSeg = ttpc[ntpc-1];
		SegPosition posAngularLastSeg = GetSegPosition(itpc2[ntpc2-1]);
		selectedPads.push_back(itpc2[ntpc2-1]);

	    Double_t rFirst = (posAngularFirstSeg.rMin+posAngularFirstSeg.rMax)/2;
	    Double_t phiFirst = (posAngularFirstSeg.phiMin+posAngularFirstSeg.phiMax)/2;
	    Double_t rLast = (posAngularLastSeg.rMin+posAngularLastSeg.rMax)/2;
	    Double_t phiLast = (posAngularLastSeg.phiMin+posAngularLastSeg.phiMax)/2;

	    delta_s = GetDistance(rFirst, phiFirst, rLast, phiLast);
	    delta_t = timeLastSeg-timeFirstSeg;
	}
	ttpc2.clear();itpc2.clear(); qtpc2.clear();
}

void GetDistanceAndTime4(Int_t ntpc, Float_t* qtpc, Int_t* itpc, Double_t* tMeantpc, Double_t* tSigmatpc, Float_t charge,
                         Double_t& delta_s, Double_t& delta_t, double& xStart, double& yStart, double& xEnd, double& yEnd, 
                         Bool_t& disregardLowChargePads, std::vector<Int_t>& disregardedPads, std::vector<Int_t>& selectedPads, 
                         Bool_t& useForAngular, std::size_t runIndex, double discardPadThreshhold = 0.003, double delta_tFac = 1.2 )
{
    disregardLowChargePads = true;
    /*
    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();

    for (const auto& row : tRawtpc) {
        for (const double x : row) {
            if (x < minVal) minVal = x;
            if (x > maxVal) maxVal = x;
        }
    }
    */
    //filter
	std::vector<Float_t> qtpc2;
    std::vector<Int_t> itpc2;
    std::vector<double> meanttpc2, sigmattpc2;
	Int_t ntpc2 = 0;
	for (std::size_t i=0; i<ntpc; ++i)
	{
		if (std::abs(qtpc[i])>discardPadThreshhold*std::abs(charge))
		{
			qtpc2.push_back(qtpc[i]);
			itpc2.push_back(itpc[i]);
            meanttpc2.push_back(tMeantpc[i]);
            sigmattpc2.push_back(tSigmatpc[i]);
            std::cout << tSigmatpc[i] << std::endl;
			++ntpc2;
		}
		else
		{
			disregardedPads.push_back(itpc[i]);
		}
    }
    std::vector<std::size_t> index = GetArgSort(meanttpc2);
    std::vector<Float_t> qtpc3;
    std::vector<Int_t> itpc3;
    std::vector<double> meanttpc3, sigmattpc3;
	Int_t ntpc3 = ntpc2;
    for (std::size_t i=0; i<ntpc2; ++i)
    {
        itpc3.push_back(itpc2[index[i]]);
        qtpc3.push_back(qtpc2[index[i]]);
        meanttpc3.push_back(meanttpc2[index[i]]);
        sigmattpc3.push_back(sigmattpc2[index[i]]);
    }
    if (ntpc2 < 2)
    {
        useForAngular = false;
        return;
    }
    bool smallTrack = false;
    double smallTrackFac = 1.;
    if (meanttpc3[0]+sigmattpc3[0]>meanttpc3[ntpc3-1]-sigmattpc3[ntpc3-1])
	{
        useForAngular = false;
	}
    if (useForAngular)
    {
        delta_t =  meanttpc3[ntpc3-1] + sigmattpc3[ntpc3-1] - (meanttpc3[0] - sigmattpc3[0]);

        std::vector<Int_t> itpcStart;
        Int_t ntpcStart = 0;
        std::vector<Int_t> itpcEnd;
        Int_t ntpcEnd = 0;
        
        for (std::size_t i=0; i<ntpc3; ++i)
        {
            //bool dominantPeakBool=2*sigmattpc3[i]>1*(maxVal-minVal);
            //bool dominantPeakBool=2*sigmattpc3[i]>0.45*(maxVal-minVal);
            //if (meanttpc3[i]<meanttpc3[0]+sigmattpc3[0] || dominantPeakBool)
            if (meanttpc3[i]<meanttpc3[0]+sigmattpc3[0])
            //if (meanttpc3[i]-sigmattpc3[i]<meanttpc3[0]+sigmattpc3[0] || dominantPeakBool)
            //if (meanttpc3[i]-0.3*sigmattpc3[i]<meanttpc3[0]+sigmattpc3[0] || dominantPeakBool)
            {
                itpcStart.push_back(itpc3[i]);
                ++ntpcStart;
            }
            //if (meanttpc3[i]>meanttpc3[ntpc3-1]-sigmattpc3[ntpc3-1] || dominantPeakBool)
            if (meanttpc3[i]+sigmattpc3[i]>meanttpc3[ntpc3-1]-sigmattpc3[ntpc3-1])
            {
                itpcEnd.push_back(itpc3[i]);
                ++ntpcEnd;
            }
        }
        for (std::size_t i=0; i<ntpcStart; ++i)
        {
            SegPosition seg = GetSegPosition(itpcStart[i]);
            double r =  (seg.rMin+seg.rMax)/2;
            double phi = (seg.phiMin+seg.phiMax)/2;
            double x,y;
            TransformCircleToCartresian(r,phi,x,y);
            xStart += x/ntpcStart;
            yStart += y/ntpcStart;
        }
        for (std::size_t i=0; i<ntpcEnd; ++i)
        {
            SegPosition seg = GetSegPosition(itpcEnd[i]);
            double r =  (seg.rMin+seg.rMax)/2;
            double phi = (seg.phiMin+seg.phiMax)/2;
            double x,y;
            TransformCircleToCartresian(r,phi,x,y);
            xEnd += x/ntpcEnd;
            yEnd += y/ntpcEnd;
        }
        delta_t = 1.2*delta_t;
        delta_s = std::sqrt((xStart-xEnd)*(xStart-xEnd)+(yStart-yEnd)*(yStart-yEnd));
    }
}



//Draw the reconstruction output
void Draw(TFile* output){
	//draws fractional error distributions
	//and reconstructed vs true graphs
	TTree *recon = (TTree*)output->Get("recon");
	TCanvas *c4 = new TCanvas("c4","Reconstruction Output");
	c4->Divide(3,2);
	c4->cd(1);
	recon->Draw("(theta_rec-theta_true)/theta_true","theta_rec!=0");
	TH1F *acc = (TH1F*)gPad->GetPrimitive("htemp");
	acc->SetName("accuracy");
	acc->SetTitle("Fractional Error: Angular Reconstruction");
	acc->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	acc->SetLineColor(kBlue+3);
	acc->SetLineWidth(3);
	acc->SetFillColor(kBlue-7);
	c4->cd(2);
	recon->Draw("(energy_rec-energy_true)/energy_true","energy_rec>0");
	TH1F *acce = (TH1F*)gPad->GetPrimitive("htemp");
	acce->SetName("accuracy");
	acce->SetTitle("Fractional Error: Energy Reconstruction");
	acce->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	acce->SetLineColor(kBlue+3);
	acce->SetLineWidth(3);
	acce->SetFillColor(kBlue-7);
	c4->cd(3);
	recon->Draw("(z_rec-z_true)/(z_true+11.5)","");
	TH1F *accz = (TH1F*)gPad->GetPrimitive("htemp");
	accz->SetName("accuracy");
	accz->SetTitle("Accuracy of Position Reconstruction");
	accz->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	accz->SetLineColor(kBlue+3);
	accz->SetLineWidth(3);
	accz->SetFillColor(kBlue-7);
	c4->cd(4);
	recon->Draw("theta_rec:theta_true","theta_rec!=0 && z_rec!=-11.5 && energy_rec>0","COLZ");
	TH2F *theta = (TH2F*)gPad->GetPrimitive("htemp");
	theta->SetTitle("Recoil Polar Angle (deg)");
	theta->GetXaxis()->SetTitle("True");
	theta->GetYaxis()->SetTitle("Reconstructed");
	c4->cd(5);
	recon->Draw("z_rec:z_true","std::abs(z_rec)<11.5 && energy_rec>0","");
	TH2F *zco = (TH2F*)gPad->GetPrimitive("htemp");
	zco->SetTitle("Event Z Position (cm)");
	zco->GetXaxis()->SetTitle("True");
	zco->GetYaxis()->SetTitle("Reconstructed");
	c4->cd(6);
	recon->Draw("energy_rec:energy_true","energy_rec>0 && z_rec!=-11.5","");
	TH2F *ke = (TH2F*)gPad->GetPrimitive("htemp");
	ke->SetTitle("Recoil Energy (MeV)");
	ke->GetXaxis()->SetTitle("True");
	ke->GetYaxis()->SetTitle("Reconstructed");
}


//Main function: read data, reconstruct, save to file, and draw
void Reconstruct(TString filename){
	//read in data file and tree
	TFile *in = new TFile(filename);
	TTree *h12 = (TTree*)in->Get("h12");
	//Get vdrift set in the simulation
    vdrift=((TParameter<Double_t>*)in->Get("TPCvDrift"))->GetVal()/10;
	//define variables for input tree information
	Int_t ntpc;							//true: number of electrons per event
	Int_t *itpc = new Int_t[100];     	//measured: id of pad
	Float_t *qtpc = new Float_t[100];  	//measured:something related to the charge 
	Float_t *ttpc = new Float_t[100];	//measured: time of impact of electrons
    Double_t *tMeantpc = new Double_t[100];
    Double_t *tSigmatpc = new Double_t[100];
	Float_t *vertex = new Float_t[3];	//true:  position of hadron
	Float_t *klab = new Float_t[3];		//true: energy of hadron  
	Float_t dircos[100][3];				//true: cosine of momentum direction in x,y,z values for each created particle during creation 
    std::vector<std::vector<double>>* tRawtpc =nullptr; //hit times
	//set branch addresses for input tree
	h12->SetBranchAddress("ntpc",&ntpc);
	h12->SetBranchAddress("itpc",itpc);
	h12->SetBranchAddress("qtpc",qtpc);
	h12->SetBranchAddress("ttpc",ttpc);
    h12->SetBranchAddress("tMeantpc",tMeantpc);
    h12->SetBranchAddress("tSigmatpc",tSigmatpc);
	h12->SetBranchAddress("klab",klab);
	h12->SetBranchAddress("vertex",vertex);
	h12->SetBranchAddress("dircos",&dircos);
    h12->SetBranchAddress("tRawtpc",&tRawtpc);

	//define variables to be used in computations
	Float_t charge;
	Int_t max = h12->GetEntries();
	//define variables for output tree information
	Float_t theta_true, z_true, energy_true;
	Float_t theta_rec, z_rec, energy_rec;
	Double_t hits[100][2];
	Bool_t disregardLowChargePads = false;
	std::vector<Int_t> disregardedPads;
	std::vector<Int_t> selectedPads;
    double xStart;
    double yStart;
    double xEnd ;
    double yEnd;
    Double_t delta_s, delta_t;
	//define output file and tree
	TFile *out = new TFile("Reconstructed.root","RECREATE","Reconstructed TPC Data");
	TTree *goat = new TTree("recon","Reconstructed TPC Data");
	//set branch addresses for output tree
	goat->Branch("theta_rec",&theta_rec,"theta_rec/F"); //reconstructed recoil theta
	goat->Branch("theta_true",&theta_true,"theta_true/F"); //simulated recoil theta
	goat->Branch("z_rec",&z_rec,"z_rec/F"); //reconstructed z coordinate
	goat->Branch("z_true",&z_true,"z_true/F"); //simulated z coordinate
	goat->Branch("energy_rec",&energy_rec,"energy_rec/F"); //reconstructed recoil energy
	goat->Branch("energy_true",&energy_true,"energy_true/F"); //simulated recoil energy
    goat->Branch("delta_s",&delta_s,"delta_s/D");
    goat->Branch("delta_t",&delta_t,"delta_t/D"); 
	goat->Branch("disregardLowChargePads", &disregardLowChargePads);
	goat->Branch("disregardedPads", &disregardedPads);
	goat->Branch("selectedPads", &selectedPads);
    goat->Branch("xStart",&xStart);
    goat->Branch("yStart",&yStart);
    goat->Branch("xEnd",&xEnd);
    goat->Branch("yEnd",&yEnd);

	Double_t mintime, maxtime, minx, maxx;
	//read input file, do computations, and write output file
	for (Int_t i=0;i<max;i++){
		h12->GetEntry(i);
		charge=0; //initialize
		//sum charge over all anode sections
		for (Int_t j=0;j<ntpc;j++){
			charge+=qtpc[j];
		}
		//get two data points for theta
		//Note index 0 = most recently saved by simulation = last to occur
		//index ntpc-1 = first saved by simulation = first to occur
		//so point 1 = "min" = tpc[ntpc-1], point 2 = "max" = tpc[0]


        Bool_t useForAngular = true;
        //std::vector<std::vector<double>>& raw = *tRawtpc;
        xStart = 0;
        yStart = 0;
        xEnd = 0;
        yEnd = 0;
		//GetDistanceAndTime1(ntpc, ttpc, qtpc, itpc, delta_s, delta_t);
		//GetDistanceAndTime3(ntpc, ttpc, qtpc, itpc, charge, delta_s, delta_t, 
        //disregardLowChargePads, disregardedPads, selectedPads, useForAngular);
        
        GetDistanceAndTime4(ntpc, qtpc, itpc, tMeantpc, tSigmatpc, charge, delta_s, delta_t, 
                            xStart, yStart, xEnd, yEnd, disregardLowChargePads, 
                            disregardedPads, selectedPads, useForAngular, i+1);
        
        //std::cout<<delta_s<< "  " << delta_t << std::endl;
		if (useForAngular)
		{
			theta_rec=GetTheta(delta_s, delta_t);
		}
		else
		{
			theta_rec=0;
		}
		z_rec=GetZ(ttpc[ntpc-1]); //central section: along beam line
		energy_rec=GetEnergy(charge);	
		//get true variables
		theta_true=acos(dircos[0][2])*180/3.14;  //changed 1 to 0, since the hadron is the first particle created in phase space mode
		z_true=vertex[2];
		energy_true=klab[0];			//changed 1 to 0, since the hadron is the first particle created in phase space mode
		//write to tree
		goat->Fill();
		disregardedPads.clear();
		selectedPads.clear();
	}
	out->Write(); //write the goat tree to a file and make the file

    std::cout << std::endl;
	Draw(out);
}