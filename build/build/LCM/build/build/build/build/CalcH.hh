/*
 This file will generate the .txt file needed for Chi Square
 Outputs data with 4 columns: r, v_lum, v_dm, h

 Created by: David Chester
 Date:6/19/2012
 */


#include <deque>
#include <string>

class Galaxy;

struct galaxyIndices
{
    int other;
    int mw;
    galaxyIndices(int oI, int mwI)
    {
        other = oI;
        mw = mwI;
    }
};

struct MeasData
{
	double r, v, vError;
};

struct CalcHData
{
    int galaxyIndex;
	double r, kappa, vLCM;
};



class CalcH{
public:
	typedef double (CalcH::*const_data_access)(int) const;
	typedef double& (CalcH::*data_access)(int);

public:
	CalcH();
	CalcH(Galaxy *milky_galaxy, Galaxy *_other_galaxy);
	~CalcH();

	void setGalaxies(Galaxy *milky_galaxy, Galaxy *other_galaxy);

	void readData(const std::string& data_filename, const std::string& dataerror_filename);

	//void makeFile(const std::string& file); Old code, read description in CalcH.cc
	void prepareDataAndGalaxies();

	void calculate_vLCM();

	int size() {return results.size();}

	Galaxy* get_milky() {return milky;}
	Galaxy* get_other() {return other;}
    
    const std::deque< galaxyIndices >& get_galIndices() const {return galIndices;}
    
    double r_f() const {return results[results.size()-1].r;}
    double kappa_f() const {return results[results.size()-1].kappa;}
    double vLCM_f() const {return results[results.size()-1].vLCM;}

    int galaxy_index(int i) const {return results[i].galaxyIndex;}
	double r(int i) const {return results[i].r;}
	double kappa(int i) const {return results[i].kappa;}
	double vLCM(int i) const {return results[i].vLCM;}

    int& galaxy_index(int i) {return results[i].galaxyIndex;}
	double& r(int i) {return results[i].r;}
	double& kappa(int i) {return results[i].kappa;}
	double& vLCM(int i) {return results[i].vLCM;}

	double data_r(int i) const {return data[i].r;}
	double data_v(int i) const {return data[i].v;}
	double data_vError(int i) const {return data[i].vError;}

	double& data_r(int i) {return data[i].r;}
	double& data_v(int i) {return data[i].v;}
	double& data_vError(int i) {return data[i].vError;}

	void setTolerance(double tol);

	// prints r and whichever variable is specified
	// usage: e.g. for kappa, Print(&CalcH::kappa);
	void Print(const_data_access print_ptr) const;

	//void printKappa() const;

private:
	void calcGalaxyInfo();

	std::deque< std::deque<double> > readfile (const std::string& filename);

	// Galaxy model information
	Galaxy* milky;
	Galaxy* other;
    
    std::deque< galaxyIndices > galIndices;

	// Actual rotation curve of the galaxy in question
	std::deque< MeasData > data;

	std::deque< CalcHData > results;

	// Tolerance for equating radial positions
	double rTolerance;
};
