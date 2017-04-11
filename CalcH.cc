/*
 This file will generate the .txt file needed for Chi Square
 Outputs data with 4 columns: r, v_lum, v_dm, h

 Created by:NOAH OBLATH y  David Chester
 Date:6/19/2012
 */

#include "CalcH.hh"
#include "Galaxy.hh"
#include "Constants.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

CalcH::CalcH() :
		milky(),
		other(),
        galIndices(),
		data(),
		results(),
		rTolerance(0.05)
//h()
{}

CalcH::CalcH(Galaxy *milky_galaxy, Galaxy *other_galaxy) :
		milky(milky_galaxy),
		other(other_galaxy),
        galIndices(),
		data(),
		results(),
		rTolerance(0.05)
{
	calcGalaxyInfo();
}
CalcH::~CalcH()
{}

void CalcH::setGalaxies(Galaxy *milky_galaxy, Galaxy *other_galaxy)
{
	milky = milky_galaxy;
	other = other_galaxy;

	calcGalaxyInfo();
	return;
}

void CalcH::calcGalaxyInfo()
{
	deque< GalaxyData >& data_mw = milky->get_data();
	deque< GalaxyData >& data_o = other->get_data();
	int lenmw = data_mw.size();
	int leno = data_o.size();


	//This checks for which galaxy has the smaller max radius and sets that as the max range for both data sets.
	double max_radius = min(data_mw.rbegin()->r, data_o.rbegin()->r);


	bool repeatedRadiiMilky = false;
    bool repeatedRadiiOther = false;
	for (int j=0; j<leno && other->r(j)<=max_radius; j++)
	{
		double diff = fabs(other->r(j) - milky->r(0));
		int i;
		for (i=1; i<lenmw && milky->r(i) <= max_radius; i++)
		{
			double newdiff = fabs(other->r(j) - milky->r(i));
			if (newdiff < diff)
			{
				diff = newdiff;
			}
			else
			{
                // this (j, i) pair are larger than the pair (j, i-1), so (j, i-1) are the pair we want,
                // as long as the difference is less than the tolerance
                if (diff > rTolerance)
                    break;
                galIndices.push_back(galaxyIndices(j, i-1));
				if (galIndices.size() > 1)
				{
					int indicesLen = galIndices.size();
					if (galIndices[indicesLen-1].other == galIndices[indicesLen-2].other)
						repeatedRadiiOther = true;
					if (galIndices[indicesLen-1].mw == galIndices[indicesLen-2].mw)
						repeatedRadiiMilky = true;
				}
				break;
			}
		}
		//once i reaches lenm or rm[i] is beyond max_radius, add the last data point to im and end the loop over j
		if (i >= lenmw || milky->r(i) > max_radius)
		{
            galIndices.push_back(galaxyIndices(j, i-1));
            if (galIndices.size() > 1)
            {
                int indicesLen = galIndices.size();
                if (galIndices[indicesLen-1].other == galIndices[indicesLen-2].other)
                    repeatedRadiiOther = true;
                if (galIndices[indicesLen-1].mw == galIndices[indicesLen-2].mw)
                    repeatedRadiiMilky = true;
            }
			break;
		}
	}

	if (repeatedRadiiMilky)
	{
		cerr << "Warning: There are repeated radii in the milky galaxy." << endl;
		for (int i=0; i<(int)galIndices.size(); i++)
		{
			cerr << i << "  " << other->r(galIndices[i].other) << "  " << milky->r(galIndices[i].mw) << endl;
		}
	}
	if (repeatedRadiiOther)
	{
		cerr << "Error: There are repeated radii in the other galaxy." << endl;
		for (int i=0; i<(int)galIndices.size(); i++)
		{
			cerr << i << "  " << other->r(galIndices[i].other) << "  " << milky->r(galIndices[i].mw) << endl;
		}
        exit(-1);
	}



	//The following commands truncates all the data sets to match the radius files, then saves it to the above 8 files.
	int lenfinal = galIndices.size();
                                 
    // remove points from the galaxies after the last matching radius
    // this way, when we ask for any final points, it will be the last matching point
    int i_prev = lenfinal - 1;
    // other galaxy
    // remove indices after lenfinal-1
    for (int i_to_remove = data_o.size()-1; i_to_remove > galIndices[i_prev].other; i_to_remove--)
    {
        deque< GalaxyData >::iterator it_to_remove = data_o.begin() + i_to_remove;
        data_o.erase(it_to_remove);
    }
                                 
    // milky galaxy
    // remove indices after lenfinal-1
    for (int i_to_remove = data_mw.size()-1; i_to_remove > galIndices[i_prev].mw; i_to_remove--)
    {
        deque< GalaxyData >::iterator it_to_remove = data_mw.begin() + i_to_remove;
        data_mw.erase(it_to_remove);
    }



	results.resize(lenfinal);
    cout << "end of galaxyInfo" << endl;
	for (int i=0; i<lenfinal; i++)
	{
        galaxy_index(i) = -1;
		r(i) = other->r(galIndices[i].other);
		kappa(i) = other->deltac(galIndices[i].other) / milky->deltac(galIndices[i].mw);
//        cout << r(i) <<"  " <<kappa(i)<< "  " << milky->vlum(galIndices[i].mw) << "  " << other->vlum(galIndices[i].other) << endl;
	}
}

deque< deque<double> > CalcH::readfile (const string& filename)
{
	ifstream indata; //ifstream reads, ofstream writes, fstream does both.
	indata.open(filename.c_str());
	if(!indata)
	{
		cerr << "Error: file could not be opened: " << filename << endl;
		exit(1);
	}
	cout << "Opened " << filename << endl;

	deque<double> x, y;
	double xtemp, ytemp;
	int count = 0;

	string line;
	getline(indata,line);
	if (isdigit(line[0]))
	{
		cout << "Didn't remove 1st line." << endl;
		stringstream linestream;
		linestream << line;
		linestream >> xtemp;
		linestream >> ytemp;
		x.push_back(xtemp);
		y.push_back(ytemp);
		count++;
	}

    indata >> xtemp;
    indata >> ytemp;
	while (! indata.eof())
	{
		x.push_back(xtemp);
		y.push_back(ytemp);
		count++;
		indata >> xtemp;
		indata >> ytemp;
	}

	/*
     for (int i=0; i<num; i++){
     cout << "x = " << x[i] << ", y = " << y[i] << endl;
     }
	 */

	indata.close();

	deque< deque<double> > answer(2);
	answer[0] = x;
	answer[1] = y;
	return answer;
}

void CalcH::readData(const string& data_filename, const string& dataerror_filename)
{
	// Read the Data file
	deque< deque<double> > readData = readfile(data_filename);

	//data.r = Data[0];
	//data.v = Data[1];
	//cout << "Successfully uploaded  data" << endl << endl;

	// Read the DMError file
	// Assumes error is upper error bar, so err must be subtracted
	deque< deque<double> > readDataError = readfile(dataerror_filename);

	if (readData[0].size() != readDataError[0].size())
	{
		cout << "Error: Data Error file is different length than Data" << endl;
		exit(-1);
	}

	data.resize(readData[0].size()); // Make data the right length.
	for (int i = 0; i<(int)data.size(); i++)
	{
		data_r(i) = readData[0][i];
		data_v(i) = readData[1][i];
		data_vError(i) = readDataError[1][i] - data_v(i);
        //cout << "data(" << i << ") = " << data_r(i) << "  " << data_v(i) << "  " << data_vError(i) << endl;
	}
    
	cout << "Successfully uploaded measured data." << endl << endl;
}

void CalcH::prepareDataAndGalaxies()
{
    // this next procedure assumes that the results vector matches the galaxy vectors
    int last_results_i = -1;
    int remove_data_after = -1; // inclusive
    int remove_data_before = -1; // inclusive
    for (int data_i = 0; data_i < data.size(); data_i++)
    {
        // find radius in results vector that best matches data_r(data_i)
        int minDiffIndex = last_results_i + 1;
        if (minDiffIndex == size())
        {
            remove_data_after = data_i;
            break;
        }
        double minDiff = fabs(data_r(data_i) - r(minDiffIndex));
        for (int results_i = minDiffIndex + 1; results_i < size(); results_i++)
        {
            double diff = fabs(data_r(data_i) - r(results_i));
            if (diff < minDiff)
            {
                minDiff = diff;
                minDiffIndex = results_i;
            }
        }
        if (minDiff > rTolerance)
        {
            remove_data_before = data_i;
            continue;
        }
        last_results_i = minDiffIndex;
        galaxy_index(minDiffIndex) = minDiffIndex;
        cout << "matching radii: data_r(" << data_i << ")=" << data_r(data_i) << " ~ " << "r(" << minDiffIndex << ")=" << r(minDiffIndex) << endl;
    }
    
    // check if we need to remove data that extends beyond the galaxies
    if (remove_data_after >= 0)
    {
        for (int data_i = data.size(); data_i > remove_data_after; data_i--)
        {
            cout << "removing data position <" << data_i << "> from the data vector" << endl;
            deque< MeasData >::iterator it_to_remove = data.begin() + data_i;
            data.erase(it_to_remove);
        }
    }
    
    // check if we need to remove data that extends before the galaxies
    if (remove_data_before >= 0)
    {
        for (int data_i = remove_data_before; data_i >= 0; data_i--)
        {
            cout << "removing data position <" << data_i << "> from the data vector" << endl;
            deque< MeasData >::iterator it_to_remove = data.begin() + data_i;
            data.erase(it_to_remove);
        }
    }
    
    // remove components of the results vector for which results.galaxyIndex = -1
    for (int results_i = size()-1; results_i >= 0; results_i--)
    {
        if (galaxy_index(results_i) == -1)
        {
            cout << "removing position <" << results_i << "> from the results vector" << endl;
            deque< CalcHData >::iterator it_to_remove = results.begin() + results_i;
            results.erase(it_to_remove);
        }
    }
    cout << "size of results vector is: " << size() << endl;
    cout << "size of the data vector is: " << data.size() << endl;
    
    cout << "testing . . ." << endl;
    for (int i=0; i<(int)data.size(); i++)
	{
        int galaxy_i = galaxy_index(i);
		//cout <<   r(i) <<   "  MW " << (   milky->vlum(galaxy_i) ) 
        //<<endl;
    }
    
    // trim results array
    
    /*
	// Make the fronts of the vectors line up in r
	if (r(0) - data_r(0) > rTolerance)
	{
		// data starts before milky
		int i;
		for (i = 1; i < (int)data.size(); i++)
		{
			if (fabs(r(0) - data_r(i)) < rTolerance) break;
		}
		if (i == (int)data.size())
		{
			cerr << "Error: no overlap between data and galaxy models (data.r < milky.r)" << endl;
			exit(-1);
		}
		// trim off the front of data
		for (int j = 0; j < i; j++)
		{
			data.pop_front();
		}
	}
	else if (data_r(0) - r(0) > rTolerance)
	{
		// milky starts before data
		int i;
		for (i = 1; i < (int)results.size(); i++)
		{
			if (fabs(data_r(0) - r(i)) < rTolerance) break;
		}
		if (i == (int)results.size())
		{
			cerr << "Error: no overlap between data and galaxy models (milky.r < data.r)" << endl;
			exit(-1);
		}
		// trim off the front of milky and other
		for (int j = 0; j < i; j++)
		{
			milky->get_data().pop_front();
			other->get_data().pop_front();
			results.pop_front();
		}
	}

	// Make the backs of the vectors line up in r
	if (data.size() > results.size())
	{
		// data is longer than milky
		int sizediff = data.size() - size();
		for (int i = 0; i < sizediff; i++)
		{
			data.pop_back();
		}
	}
	else if (results.size() > data.size())
	{
		// milky is longer than data
		int sizediff = results.size() - data.size();
		for (int i = 0; i < sizediff; i++)
		{
			milky->get_data().pop_back();
			other->get_data().pop_back();
			results.pop_back();
		}
	}
     */

	cout << "Successfully calculated everything needed for Chi Square" << endl << endl;
}

//makeFile shouldn't be used once code runs from Fitting Program, though is still referenced for past functions and construction. Removed from compiling on 1/31/14.

/*void CalcH::makeFile(const string& file)
{
    const double c = Constants::c;
    const double overc = 1. / c;

    trimVectors();

    if (milky.r.size() != data.r.size())
    {
        cerr << "Error: for some reason the galaxy vector size is not equal to the data vector size: " << endl;
 */
//other->get_data().rbegin()->phi



void CalcH::calculate_vLCM()
{
	const double c = Constants::c;
	const double overc = 1. / c;

	if (results.size() != data.size())
	{
		cerr << "Error: for some reason the galaxy vector size is not equal to the data vector size: " << endl;

		cerr << "       results.size = " << results.size() << endl;
		for (int i=0; i<(int)results.size(); i++)
		{
			cerr << i << "  " << r(i) << endl;
		}

		cerr << "       data.size = " << data.size() << endl;
		for (int i=0; i<(int)data.size(); i++)
		{
			cerr << i << "  " << data_r(i) << endl;
		}
		exit(-1);
	}
    
    //cout << "the final value for other vlum is: " << other->get_data().rbegin()->vlum << endl;

    // first update kappa (because the other galaxy's scalings may have changed
	for (int i=0; i<(int)data.size(); i++)
	{
        int galaxy_i = galaxy_index(i);
        //DOES NOT GIVE RIGHT VALUES:kappa(i) ne other->deltac(i) / milky->deltac(i);
        kappa(i) = other->deltac(galIndices[galaxy_i].other) / milky->deltac(galIndices[galaxy_i].mw);
       cout <<milky->r(i)<< "   "<<milky->vlum(i)<<"   "<<kappa_f()<< endl;
       //kappa(i)<<"   "<<kappa_f()<< endl;
         //cout << "deltaC(MW, other)" << i<<"  " <<c*(1.-1./milky->n(i))<<  "  " <<c*(1.-1./other->n(i))<<endl;
    }
    
    //cout << "R,milky, other.phi: "<< r(0) <<  "  "  << milky->phi(0) << "  "  << other->phi(0) << endl;
   // cout << "kappa: " << kappa(8) << endl;
    // then calculate vLCM
	 for (int i=0; i<(int)data.size(); i++)
	{
        int mw_i = galIndices[galaxy_index(i)].mw;
        int other_i = galIndices[galaxy_index(i)].other;
        //LCM error func (v_gal/gal)
	//This is the current version as of 2017/2/11
	//vLCM(i) =c*c*((1.-1./other->n(other_i))/(1.-1./milky->n(mw_i)))*((1.-1./other->n(other_i))/(1.-1./milky->n(mw_i)))*(((sqrt((1.+other->vlum(other_i)*overc)/(1.-other->vlum(other_i)*overc)))+(other->n(other_i)/milky->n(mw_i)))/((sqrt((1.+other->vlum(other_i)*overc)/(1.-other->vlum(other_i)*overc)))-(other->n(other_i)/milky->n(mw_i))))*((2/(other->n(other_i)/milky->n(mw_i)+milky->n(mw_i)/other->n(other_i)))-1);
	//rott - 2017/2/11 - refactoring to be more readable and stable
	//first, define w = n - 1   for computational stability
	//Note that the final value can differ from the vLCM(i) above by
	//about one part in 10000, this is because the new f_i is more
	//numerically stable
	double w_gal = other->n(other_i) - 1;
	double w_mw = milky->n(mw_i) - 1;
	double kappa_i = (w_gal/w_mw) * (w_mw + 1) / (w_gal + 1);
	double a_i = kappa_i * kappa_i;
	//d_i may still have problems, need to investigate more
	double d_i = (((sqrt((1.+other->vlum(other_i)*overc)/(1.-other->vlum(other_i)*overc)))+(other->n(other_i)/milky->n(mw_i)))/((sqrt((1.+other->vlum(other_i)*overc)/(1.-other->vlum(other_i)*overc)))-(other->n(other_i)/milky->n(mw_i))));
	double f_i = -(w_gal - w_mw)*(w_gal - w_mw) / ((w_gal+1)*(w_gal+1) + (w_mw+1)*(w_mw+1));
	double vLCM_new = c*c*a_i*d_i*f_i;

        //remember: these guys must be calling the correct galaxy index or they mismatch.
 // (((other->n(other_i)/milky->n(mw_i)) - 1./(other->n(other_i)/milky->n(mw_i)))/((other->n(other_i)/milky->n(mw_i)) + 1./(other->n(other_i)/milky->n(mw_i))));
     cout <<other->r(other_i)<< "   "<<other->vlum(other_i)<< endl;
	//cout <<other->r(other_i) <<  "  "  <<milky->vlum(mw_i) <<endl;
  //((1.-1./other->n(other_i))/(1.-1./milky->n(mw_i)))<<  "  "  << kappa(i)<< "  "  
//	cout <<other->deltac(other_i) / milky->deltac(mw_i)<< "  "  <<kappa_f();
  //(other->phi(other_i)/milky->phi(mw_i)) <<endl;

	}
    //Print(&CalcH::vLCM);
    //Print(&Galaxy::vlum);
	//cout << "Successfully calculated for vLCM" << endl;

}

void CalcH::setTolerance(double tol)
{
	rTolerance = tol;
}

void CalcH::Print(const_data_access print_ptr) const
{
    int length = data.size();
    for (int i=0; i<length; i++)
    {
        cout << r(i) <<" "<< (this->*print_ptr)(i) << endl;
    }
}
