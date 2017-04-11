// This is a "simple" numerical fitter using TMinuit
// It assumes the data is in the form of four columns, you can
// of course change this by altering the "getData" function
// Be forwarned: without errors on the data (which were not provided),
// the error bars on the fit values are basically meaningless

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "TMinuit.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TAxis.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include "CalcH.hh"
#include "Galaxy.hh"
#include "GalaxyDB.hh"




using namespace std;

// global pointer to the calcH object for minuit's use
CalcH* g_calcH;


// load the input files for a galaxy
bool LoadGalaxyFiles(const GalaxyDB::GalaxyFiles& files, Galaxy& galaxy)
{
	for (GalaxyDB::GalaxyFiles::const_iterator fileIt = files.begin(); fileIt != files.end(); ++fileIt)
	{
		if (fileIt->second.empty())
		{
			cerr << "empty filename for galaxy file with label " << fileIt->first << endl;
			return false;
		}
		galaxy.readfile(fileIt->second, fileIt->first);
	}
	return true;
}


// chi-squared calculation, for fitting
void fcn(int& npar, double* /*deriv*/, double& f, double par[], int /*flag*/)
{
    //cout << "in fcn: " << npar << "  " << par[0] << "  " << par[1] << "  " << par[2] << "  " << par[3];
    int ncomp = g_calcH->get_other()->get_n_components();
    for (int i=0; i<ncomp; i++)
    {
    	g_calcH->get_other()->setScaling(i, par[i]);
    }
    g_calcH->get_other()->calcLumScaled();
    g_calcH->get_other()->calcPhiB();
    g_calcH->get_other()->calcDeltaC();
    g_calcH->calculate_vLCM();

    Double_t alpha = par[ncomp];
    
    Int_t nEntries = g_calcH->size();
    //Compute chi2, you can do this as a separate function if you want
    Double_t chi2 = 0;
    Double_t vSum, diff;
    for (int i=0; i<nEntries; i++)
    {
        int galaxy_i = g_calcH->galaxy_index(i);
        int other_i = g_calcH->get_galIndices()[galaxy_i].other;
        vSum = g_calcH->get_other()->vlum(other_i) * g_calcH->get_other()->vlum(other_i) + alpha * g_calcH->vLCM(i);
        diff = sqrt(vSum) - g_calcH->data_v(i);
        chi2 += diff * diff / (g_calcH->data_vError(i)  * g_calcH->data_vError(i));
        //cout <<  "  "<< g_calcH->vLCM(i);
    }
    f = chi2;
    //cout << "  " << f << endl;
}                         // end of fcn


int KDFit3(string galaxyLabel, string galaxyDBFilename, string graphfile)
{
    
    //To use this, simply load the macro into ROOT using .L KDFit.C
    //Then run KDFit("data filename", "filename to store graphs");
    //The filename to store graphs must end in .ps, .pdf, .png, .jpg
    //or a couple of others (check ROOT's documentation)

	// Load the galaxy DB
	GalaxyDB allGalaxies(galaxyDBFilename, "Label");
	if (! allGalaxies.DBIsReady())
	{
		cerr << "Error: galaxy DB was not loaded" << endl;
		return -1;
	}

	// Ensure that the Milky Way is in the DB
	string mwLabel("MilkyWay");
	if (! allGalaxies.HasGalaxy(mwLabel))
	{
		cerr << "Error: Milky Way was not in the galaxy DB" << endl;
		return -1;
	}

	// Ensure that the other galaxy is in the DB
	if (! allGalaxies.HasGalaxy(galaxyLabel))
	{
		cerr << "Error: unknown galaxy: " << galaxyLabel << endl;
		return -1;
	}


    // Load the other galaxy
    Galaxy other;
	GalaxyDB::GalaxyFiles otherFiles = allGalaxies.GetFiles(galaxyLabel);
    if (! LoadGalaxyFiles(otherFiles, other))
    {
    	cerr << "Error: files did not load for galaxy <" << galaxyLabel << ">" << endl;;
        return -1;
    }
    other.initialCalcLum();
    other.calcLumScaled();
    other.calcPhiB();
    other.calcDeltaC();

    // Load the Milky Way
    Galaxy milky;
	GalaxyDB::GalaxyFiles milkyFiles = allGalaxies.GetFiles(mwLabel);
    if (! LoadGalaxyFiles(milkyFiles, milky))
    {
    	cerr << "Error: files did not load for galaxy <" << mwLabel << ">" << endl;
        return -1;
    }
    milky.initialCalcLum();
    milky.calcLumScaled();
    milky.calcPhiB();
    milky.calcDeltaC();

	// Load dataFile, dataErrorFile, and outFile and make sure they're present
	string dataFile = allGalaxies.GetDataFile(galaxyLabel, "Data");
	string dataErrorFile = allGalaxies.GetDataFile(galaxyLabel, "DataError");

	if (dataFile.empty() || dataErrorFile.empty())
	{
		cerr << "One of the required files is missing:\n" <<
				"\tData: " << dataFile << '\n' <<
				"\tDataError: " << dataErrorFile << '\n' << endl;
		return -1;
	}

	CalcH calcH( &milky, &other );
	calcH.readData(dataFile, dataErrorFile);
	calcH.prepareDataAndGalaxies();
    //exit(-1);

	g_calcH = &calcH;

    
    //Change this to use errors
    bool error = true;
    //bool error = false;
    
    // Initialize minuit, set initial values etc. of parameters.
    
    int ncomp = g_calcH->get_other()->get_n_components();
    int n_other_par = 1;
    int npar = ncomp + n_other_par;
    TMinuit minuit(npar);
    minuit.SetFCN(fcn);
    
    for (int i = 0; i < ncomp; i++)
    {
        if (calcH.get_other()->get_component_label(i) == "Gas")
        {
            minuit.DefineParameter(i, calcH.get_other()->get_component_label//(i).c_str(), 1., 0.01, 1.0, 1.0);
                 (i).c_str(), 1., 0.01, 0.75, 1.5);                  
        }
        else
        {
            minuit.DefineParameter(i, calcH.get_other()->get_component_label(i).c_str(), 1., 0.01, 0.000000000001, 5.0);
        }
    }
 minuit.DefineParameter(ncomp, "alpha", 1.00, 0.0001, 0.0000, 0.0000);
 //  minuit.DefineParameter(ncomp, "alpha", 0.08, 0.0001, 0.7100, 0.9000);
    //minuit.FixParameter(ncomp);
    // Fix parameter beta
   
    
    // Do the minimization!
    minuit.SetMaxIterations(2);
    minuit.Migrad();       // Minuit's best minimization algorithm
    //minuit.mnmnos();

    // Get results out for plots, etc
    
    calcH.Print(&CalcH::kappa);
    //calcH.Print(&CalcH::DeltaC);
    

    double* outpar = new double[npar];
    double* outparerr = new double[npar];
    for (int i=0; i<npar; i++)
    {
        minuit.GetParameter(i,outpar[i],outparerr[i]);
        //minuit.mnerrs(i, ePlus, eMinus, eParab, gcc);
        //cout << parName[i] << " + " << ePlus << " - " << eMinus << endl;
    }

    // alpha (the first parameter after the scaling components)
    double alpha = outpar[ncomp];
    double alphaerr = outparerr[ncomp];

    //Compute chi2 for best fit point
    for (int i=0; i<ncomp; i++)
    {
    	calcH.get_other()->setScaling(i, outpar[i]);
    }
    calcH.get_other()->calcLumScaled();
    calcH.get_other()->calcPhiB();
    calcH.get_other()->calcDeltaC();
    calcH.calculate_vLCM();

    Int_t nEntries = calcH.size();

    //Compute chi2
    Double_t chi2 = 0;
    Double_t ldm, vSum, diff;
    for (int i=0; i<nEntries; i++)
    {
        int galaxy_i = g_calcH->galaxy_index(i);
        int other_i = g_calcH->get_galIndices()[galaxy_i].other;
        vSum = g_calcH->get_other()->vlum(other_i) * g_calcH->get_other()->vlum(other_i) + alpha * g_calcH->vLCM(i);
        ldm = alpha * g_calcH->vLCM(i);
        diff = sqrt(vSum) - g_calcH->data_v(i);
        chi2 += diff * diff / (calcH.data_vError(i)  *calcH.data_vError(i));
    }
    Double_t dof = nEntries - minuit.GetNumFreePars();
    Double_t c = chi2 ;
    
    //Print results to screen
    for (int i=0; i < ncomp; i++)
    {
		cout << "scaling " << i << " = " << outpar[i] << " +- " << outparerr[i] << endl;
    }
    cout << "alpha = " << alpha << " +- " << alphaerr << endl;
    cout << "dof = " << dof << endl;
    cout << "chi2/dof for best fit = " << c << "/" << dof << endl;
    
    delete [] outpar;
    delete [] outparerr;
    
    //cout << "phi_0 other is: " << calcH.get_other()->get_data().rbegin()->phi << endl;
    // make a plot!
    TCanvas* canvas = new TCanvas();

    // Set a bunch of parameters to make the plot look nice

    canvas->SetFillColor(0);
    canvas->UseCurrentStyle();
    canvas->SetBorderMode(0);       // still leaves red frame bottom and right
    canvas->SetFrameBorderMode(0);   // need this to turn off red hist frame!
    gROOT->SetStyle("Plain");
    canvas->UseCurrentStyle();
    gROOT->ForceStyle();

    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleSize(0.04);
    gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
    gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
    gROOT->ForceStyle();

    // load graphs
    TGraphErrors* dataGraph = new TGraphErrors(nEntries);
    TGraph* modelGraph = new TGraph(nEntries);
    TGraph* vlumGraph = new TGraph(nEntries);
    TGraph* ldmGraph = new TGraph(nEntries);
    for (int i=0; i<nEntries; i++)
    {
        int galaxy_i = g_calcH->galaxy_index(i);
        int other_i = g_calcH->get_galIndices()[galaxy_i].other;
        vSum = g_calcH->get_other()->vlum(other_i) * g_calcH->get_other()->vlum(other_i) + alpha * g_calcH->vLCM(i);
        ldm =  alpha * g_calcH->vLCM(i);
        ldmGraph->SetPoint(i, calcH.r(i), sqrt(ldm));
        modelGraph->SetPoint(i, calcH.r(i), sqrt(vSum));
        vlumGraph->SetPoint(i, calcH.r(i), calcH.get_other()->vlum(other_i));
        dataGraph->SetPoint(i, calcH.data_r(i), calcH.data_v(i));
        if (error)
        {
        	dataGraph->SetPointError(i, 0., calcH.data_vError(i));
        }
    }
    /*ldmGraph->SetLineColor(kBlue);
    ldmGraph->SetLineStyle(4);
    ldmGraph->SetLineWidth(2);
    ldmGraph->SetFillColor(kWhite);
    
    
    modelGraph->SetLineColor(kRed);
    modelGraph->SetLineStyle(5);
    modelGraph->SetLineWidth(3);
    modelGraph->SetFillColor(kWhite);

    vlumGraph->SetLineColor(kBlack);
    vlumGraph->SetLineStyle(3);
    vlumGraph->SetLineWidth(3);
    

    dataGraph->SetLineColor(kBlack);
    dataGraph->SetLineStyle(9);
    dataGraph->SetMarkerColor(kBlack);
    dataGraph->SetMarkerStyle(24);
    dataGraph->SetMarkerSize(0.602);
    dataGraph->SetFillColor(kWhite);*/
     //This is for the Black and white pix for journal's print copy
    ldmGraph->SetLineColor(kBlack);
    ldmGraph->SetLineStyle(4);
    ldmGraph->SetLineWidth(2);
    ldmGraph->SetFillColor(kWhite);
    
    
    modelGraph->SetLineColor(kBlack);
    modelGraph->SetLineStyle(5);
    modelGraph->SetLineWidth(3);
    modelGraph->SetFillColor(kWhite);
    
    vlumGraph->SetLineColor(kBlack);
    vlumGraph->SetLineStyle(3);
    vlumGraph->SetLineWidth(3);
    
    
    dataGraph->SetLineColor(kBlack);
    dataGraph->SetLineStyle(9);
    dataGraph->SetMarkerColor(kBlack);
    dataGraph->SetMarkerStyle(24);
    dataGraph->SetMarkerSize(0.602);
    dataGraph->SetFillColor(kWhite);
    
    Double_t xMin=999999., xMax=0.;
    Double_t yMin=0.0, yMax=0.;
    
    ldmGraph->Draw("APL");
    canvas->Update();
    xMin = TMath::Min(canvas->GetUxmin(), xMin);
    xMax = TMath::Max(canvas->GetUxmax(), xMax);
    yMax = TMath::Max(canvas->GetUymax(), yMax);
    
    modelGraph->Draw("APL");
    canvas->Update();
    xMin = TMath::Min(canvas->GetUxmin(), xMin);
    xMax = TMath::Max(canvas->GetUxmax(), xMax);
    yMax = TMath::Max(canvas->GetUymax(), yMax);
    
    vlumGraph->Draw("APL");
    canvas->Update();
    xMin = TMath::Min(canvas->GetUxmin(), xMin);
    xMax = TMath::Max(canvas->GetUxmax(), xMax);
    yMax = TMath::Max(canvas->GetUymax(), yMax);
    
    dataGraph->Draw("APL");
    canvas->Update();
    xMin = TMath::Min(canvas->GetUxmin(), xMin);
    xMax = TMath::Max(canvas->GetUxmax(), xMax);
    yMax = TMath::Max(canvas->GetUymax(), yMax);
    
    
    
    // final drawing
    
    TH1I* hAxes = new TH1I("hAxes", "", 1, xMin, xMax);
    hAxes->GetYaxis()->SetRangeUser(yMin, yMax);
    hAxes->SetXTitle(" R (kpc)");
    hAxes->SetYTitle(" V (km/s)");
    gStyle->SetOptTitle(0);
    hAxes->Draw();
    
    ldmGraph->Draw("PL");
    modelGraph->Draw("PL");
    vlumGraph->Draw("PL");
    dataGraph->Draw("zP");
    
    
    //TLegend *mylegend = new TLegend(0.6,0.90,0.8,0.999);
    //mylegend->AddEntry(&dataGraph, "Measured");
    //mylegend->AddEntry(&modelGraph, "LCM");
    //mylegend->AddEntry(&Lorenz_unscaledGraph,"Luminous "); //"NFW");
    //mylegend->SetFillColor(kWhite);
    //mylegend->Draw();
    
    canvas->SaveAs(graphfile.c_str());
    
    canvas->Close();
    
    return 0;
    
}


int main(int argc, char** argv)
{
    if (argc < 4)
	{
	    cout << "Usage: KDFit3e [galaxy label] [galaxy db file] [plot file]" << endl;
	    return 0;
	}

    string galaxyLabel(argv[1]);
    string galaxyDBFile(argv[2]);
    string plotFile(argv[3]);

    return KDFit3(galaxyLabel, galaxyDBFile, plotFile);
}
