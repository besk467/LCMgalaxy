/*
 This file will generate the .txt file needed for Chi Square
 Outputs data with 4 columns: r, v_lum, v_dm, h

 Created by: Noah Oblath & David Chester, MIT 2012-14
 Date:6/19/2012
 */

#include "Galaxy.hh"
#include "Constants.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>


using namespace std;

Galaxy::Galaxy() :
		n_components(0),
		allRawData(),
		data(),
		scalings()
{
}

Galaxy::~Galaxy()
{}

void Galaxy::readfile (const string& filename, const string& label)
{
	ifstream indata; //ifstream reads, ofstream writes, fstream does both.
	indata.open(filename.c_str());
	if(!indata){
		cerr << "Error: file could not be opened: " << filename << endl;
		exit(1);
	}
	cout << "Opened " << filename << endl;

	GalaxyRawData newRawData;

	newRawData.label = label;

	double r_temp, vlum_temp;
	int count = 0;

	string line;
	getline(indata,line);
	if (isdigit(line[0]))
	{
		cout << "Didn't remove 1st line." << endl;
		stringstream linestream;
		linestream << line;
		linestream >> r_temp;
		linestream >> vlum_temp;
		newRawData.r.push_back(r_temp);
		newRawData.vlum.push_back(vlum_temp);
		count++;
	}

	while (! indata.eof())
	{
        r_temp = -1.;
        vlum_temp = -1.;
        indata >> r_temp;
        indata >> vlum_temp;
        if (r_temp < 0 || vlum_temp < 0)
        {
            continue;
        }
		newRawData.r.push_back(r_temp);
		newRawData.vlum.push_back(vlum_temp);
       // cout << newRawData.r.back() << "  " << newRawData.vlum.back()<< endl;
		count++;
	}


	/*
     for (int i=0; i<num; i++){
     cout << "x = " << x[i] << ", y = " << y[i] << endl;
     }
	 */

	indata.close();

	allRawData.push_back(newRawData);

	scalings.push_back(1.);

	++n_components;

	return;
}


void Galaxy::initialCalcLum()
{
	if (allRawData.empty())
	{
		cerr << "No raw data has been loaded" << endl;
		exit(-1);
	}

	unsigned lumSize = allRawData[0].r.size();
	for (unsigned iFile = 1; iFile < allRawData.size(); ++iFile)
	{
		if (allRawData[iFile].r.size() != lumSize)
		{
			cerr << "Raw data sizes are unequal" << endl;
			exit(-1);
		}
	}

	data.resize(lumSize);

	for (unsigned i=0; i<lumSize; i++)
	{
		unsigned count = 0;
		double rSum = 0.;
		double rSqSum = 0.;
		double vlumSum = 0.;
        data[i].vlum_components.resize(n_components);
		for (unsigned iFile = 0; iFile < n_components; ++iFile)
		{
			++count;
			rSum += allRawData[iFile].r[i];
			rSqSum += allRawData[iFile].r[i] * allRawData[iFile].r[i];
			data[i].vlum_components[iFile] = allRawData[iFile].vlum[i];
			vlumSum += data[i].vlum_components[iFile] * data[i].vlum_components[iFile];
		}
		r(i) = rSum / (double)count;
		dr(i) = sqrt(rSqSum / (double)count - r(i) );

		vlum_unscaled(i) = sqrt(vlumSum);
	}

	cout << "Successfully calculated v_lum_unscaled in quadrature" << endl << endl;

	return;
}

void Galaxy::calcLumScaled()
{
	for (int i=0; i<(int)data.size(); i++)
	{
		double vlumSum = 0.;
		for (unsigned iFile = 0; iFile < n_components; ++iFile)
		{
			vlumSum +=  data[i].vlum_components[iFile] * data[i].vlum_components[iFile] * scalings[iFile];
		}

        vlum(i) = sqrt(vlumSum);
	}
}

void Galaxy::calcPhiB()
{
	const double oneOverG = 1. / Constants::G;
	const double c = Constants::c;
	mass(0) = vlum(0) * vlum(0) * r(0) * oneOverG;
	for (int i=1; i<(int)data.size(); i++)
	{ //void Galaxy::printMass() const
        //{cout<<"mass: "<< endl; Print(mass);}
		mass(i) = vlum(i)*vlum(i)*r(i)*oneOverG;
        //cout<<"mass: "<< endl; Print(mass); 
	}
	for (int i=0; i<(int)data.size(); i++)
	{
		force(i) = 2.0*vlum(i)*vlum(i)/(r(i)*c*c) ;
		angmom(i) = 3.0*vlum(i)*vlum(i)*vlum(i)*r(i);
	}
	integrate(&Galaxy::force, &Galaxy::phi);
	for (int i=0; i<(int)data.size(); i++) 
	{//cout<<"test: " <<i<< "  " << force(i)<< "  " <<vlum(i)<< "  "  <<r(i)<< "  "  << phi(i) << "  "   << endl;
		angmom(i) = 3.0*phi(i)*c*c*vlum(i)*r(i);
	}
	integrate(&Galaxy::angmom, &Galaxy::J);
	for (int i=0; i<(int)data.size(); i++){
		b(i) = J(i)/(r(i)*r(i)*phi(i)*c*c*c);
	}
}

void Galaxy::integrate(const_data_access ptr_in, data_access ptr_out)
{
	int len = data.size();
	(this->*ptr_out)(0) = 0.5 * r(0) * (this->*ptr_in)(0);
	for (int i=1; i<len; i++){
		(this->*ptr_out)(i) = 0.5 * ((this->*ptr_in)(i) + (this->*ptr_in)(i-1)) * (r(i) - r(i-1)) + (this->*ptr_out)(i-1);
	}
	return;
}

void Galaxy::calcDeltaC()
{
	const double c = Constants::c;
	for(int i=0; i<(int)data.size(); i++){
		gtp(i) = -1.0*r(i)*phi(i)*b(i);
		gtt(i) = -1.0*(1.0-phi(i));
		gpp(i) = r(i)*r(i)*(1.0+b(i)*b(i)*(1.0+phi(i)));
		n(i)=1/sqrt(-gtt(i));
		deltac(i) = c*(1- sqrt(-gtt(i)));
        //cout <<  i <<" "<< deltac(i) << endl;
	}
}

void Galaxy::calcError()
{
	const double c = Constants::c;
	int size = data.size();
	deque< double > DrO(size);
	for (int i=0; i<(int)data.size(); i++)
	{
		if (i==0){dvlum(i) = (vlum(i+1)-vlum(i))/(r(i+1)-r(i))*dr(i);}
		else {dvlum(i) = (vlum(i)-vlum(i-1))/(r(i)-r(i-1))*dr(i);}
		dforce(i) = sqrt(pow(4.*vlum(i)*dvlum(i)/(c*c*r(i)),2)
				+ pow(4.*vlum(i)*vlum(i)*dr(i)/(c*c*r(i)*r(i)),2));
		if (i==0){dphi(i) = pow(r(i+1)*dforce(i)/2.0,2)
			+ pow(force(i+1)*dr(i)/2.0,2);}
		else if (i==((int)data.size()-1)){dphi(i) = dphi(i-1) + pow(r(i-1)*dforce(i)/2.0,2) + pow((force(i-1)+force(i))*dr(i)/2.0,2);}
		else {dphi(i) = dphi(i-1) + pow((r(i+1)+r(i-1))*dforce(i)/2.0,2)
			+ pow((force(i-1)-force(i+1))*dr(i)/2.0,2);}
		dangmom(i) = sqrt(pow(9.0*vlum(i)*vlum(i)*r(i)*dvlum(i),2)
				+ pow(3.0*vlum(i)*vlum(i)*vlum(i)*dr(i),2));
		if (i==0){dJ(i) = pow(r(i+1)*dangmom(i)/2.0,2)
			+ pow(angmom(i+1)*dr(i)/2.0,2);}
		else if (i==((int)data.size()-1)){dJ(i) = dJ(i-1) + pow(r(i-1)*dangmom(i)/2.0,2) + pow((angmom(i-1)+angmom(i))*dr(i)/2.0,2);}
		else {dJ(i) = dJ(i-1) + pow((r(i+1)+r(i-1))*dangmom(i)/2.0,2)
			+ pow((angmom(i-1)-angmom(i+1))*dr(i)/2.0,2);}
		db(i) = sqrt(pow(dJ(i)/(c*pow(r(i)*vlum(i),2)),2)
				+ pow(2.0*J(i)*dr(i)/(pow(r(i)*vlum(i),2)*r(i)*c),2)
				+ pow(2.0*J(i)*dvlum(i)/(pow(r(i)*vlum(i),2)*vlum(i)*c),2));
		dgtp(i) = sqrt(pow(phi(i)*b(i)*dr(i),2) + pow(r(i)*b(i)*dphi(i),2)
				+ pow(phi(i)*r(i)*db(i),2));
		dgtt(i) = dphi(i);
		dgpp(i) = sqrt(pow(2.0*r(i)*(1.0+b(i)*b(i)*(1.0+phi(i)))*dr(i),2)
				+ pow(2.0*r(i)*r(i)*b(i)*(1.0+phi(i))*db(i),2)
				+ pow(pow(r(i)*b(i),2)*dphi(i),2));
		dn(i) = sqrt(pow((1.0+gtp(i)/sqrt(gtp(i)*gtp(i)-gtt(i)*gpp(i)))*dgtp(i)/(gtt(i)*sqrt(gpp(i))),2)
				+ pow((1.0/(2.0*sqrt(gpp(i)*(gtp(i)*gtp(i)-gtt(i)*gpp(i))))+(gtp(i)+sqrt(gtp(i)*gtp(i)-gtt(i)*gpp(i)))/(2.0*gtt(i)*pow(gpp(i),1.5)))*dgpp(i),2)
				+ pow((sqrt(gpp(i)/(gtp(i)*gtp(i)-gtt(i)*gpp(i)))/(2.0*gtt(i))+(gtp(i)+sqrt(gtp(i)*gtp(i)-gtt(i)*gpp(i)))/(sqrt(gpp(i))*gtt(i)*gtt(i)))*dgtt(i),2));
		DrO[i] = sqrt(pow(dn(i)/(n(i)*n(i)),2)/(1.0+b(i)*b(i)*(1+phi(i)))
				+ pow(b(i)*(1.0+phi(i))*db(i)/n(i),2)/pow(1.0+b(i)*b(i)*(1.0+phi(i)),3)
				+ pow(b(i)*b(i)*dphi(i)/(2.0*n(i)),2)/pow(1.0+b(i)*b(i)*(1.0+phi(i)),3));
		ddeltac(i) = c*DrO[i];
	}
	// Still need to square root Dphi and DJ
	for (int i=0; i<(int)data.size(); i++)
	{
		dphi(i) = sqrt(dphi(i));
		dJ(i) = sqrt(dJ(i));
	}
}

void Galaxy::setScaling(int component, double scaling)
{
	scalings[component] = scaling;
	return;
}

void Galaxy::printError()
{
	for (int i=0; i<(int)data.size(); i=i+10){cout << "r["<<i<<"] = " <<r(i)<<" +/- "<<dr(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "v["<<i<<"] = " <<vlum(i)<<" +/- "<<dvlum(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "F["<<i<<"] = " <<force(i)<<" +/- "<<dforce(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "phi["<<i<<"] = " <<phi(i)<<" +/- "<<dphi(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "ang["<<i<<"] = " <<angmom(i)<<" +/- "<<dangmom(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "J["<<i<<"] = " <<J(i)<<" +/- "<<dJ(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "b["<<i<<"] = " <<r(i)<<" +/- "<<db(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "gtp["<<i<<"] = " <<gtp(i)<<" +/- "<<dgtp(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "gtt["<<i<<"] = " <<gtt(i)<<" +/- "<<dgtt(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "gpp["<<i<<"] = " <<gpp(i)<<" +/- "<<dgpp(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "n["<<i<<"] = " <<n(i)<<" +/- "<<dn(i)<< endl;}
	for (int i=0; i<(int)data.size(); i=i+10){cout << "dc["<<i<<"] = " <<deltac(i)<<" +/- "<<ddeltac(i)<< endl;}
}


void Galaxy::Print(const_data_access print_ptr) const
{
    int length = data.size();
    for (int i=0; i<length; i++)
    {
        cout << r(i) <<" "<< (this->*print_ptr)(i) << endl;
    }
}
