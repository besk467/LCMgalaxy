/*
 This file will generate the .txt file needed for Chi Square
 Outputs data with 4 columns: r, v_lum, v_dm, h

 Created by: Noah Oblath & David Chester 2012-14, MIT
 Date:6/19/2012
 */


#include <string>
#include <deque>

struct GalaxyRawData
{
	std::deque<double> r;
	std::deque<double> vlum;
	std::string label;
    GalaxyRawData() {};
    ~GalaxyRawData() {};
};

struct GalaxyData
{
	std::deque<double> vlum_components;
	double r, vlum, vlum_unscaled, force, phi, angmom, J, b, gtp, gtt, gpp, n, deltac, mass;
	double dr, dvlum, dforce, dphi, dangmom, dJ, db, dgtp, dgtt, dgpp, dn, ddeltac;
    GalaxyData() :
        vlum_components()
    {};
    ~GalaxyData() {};
};



class Galaxy
{
public:
	typedef double (Galaxy::*const_data_access)(int) const;
	typedef double& (Galaxy::*data_access)(int);

public:
	Galaxy();
	~Galaxy();

	void readfile (const std::string& filename, const std::string& label);

	void initialCalcLum();

	void calcLumScaled();

	void calcPhiB();

	void calcDeltaC();

	void calcError();

	void setScaling(int component, double scaling);

	int get_n_components() const {return n_components;}

	const std::string& get_component_label(int i) {return allRawData[i].label;}

	int size() const {return data.size();}

	const std::deque< GalaxyData >& get_data() const {return data;}
	std::deque< GalaxyData >& get_data() {return data;}

    double r_f() const {return data[data.size()-1].r;}
    double vlum_f() const {return data[data.size()-1].vlum;}
    double phi_f() const {return data[data.size()-1].phi;}
    
	double r(int i) const {return data[i].r;}
	double vlum(int i) const {return data[i].vlum;}
	double vlum_unscaled(int i) const {return data[i].vlum_unscaled;}
	double force(int i) const {return data[i].force;}
	double phi(int i) const {return data[i].phi;}
	double angmom(int i) const {return data[i].angmom;}
	double J(int i) const {return data[i].J;}
	double b(int i) const {return data[i].b;}
	double gtp(int i) const {return data[i].gtp;}
	double gtt(int i) const {return data[i].gtt;}
	double gpp(int i) const {return data[i].gpp;}
	double n(int i) const {return data[i].n;}
	double deltac(int i) const {return data[i].deltac;}
	double mass(int i) const {return data[i].mass;}

	double& r(int i) {return data[i].r;}
	double& vlum(int i) {return data[i].vlum;}
	double& vlum_unscaled(int i) {return data[i].vlum_unscaled;}
	double& force(int i) {return data[i].force;}
	double& phi(int i) {return data[i].phi;}
	double& angmom(int i) {return data[i].angmom;}
	double& J(int i) {return data[i].J;}
	double& b(int i) {return data[i].b;}
	double& gtp(int i) {return data[i].gtp;}
	double& gtt(int i) {return data[i].gtt;}
	double& gpp(int i) {return data[i].gpp;}
	double& n(int i) {return data[i].n;}
	double& deltac(int i) {return data[i].deltac;}
	double& mass(int i) {return data[i].mass;}

	double dr(int i) const {return data[i].dr;}
	double dvlum(int i) const {return data[i].dvlum;}
	double dforce(int i) const {return data[i].dforce;}
	double dphi(int i) const {return data[i].dphi;}
	double dangmom(int i) const {return data[i].dangmom;}
	double dJ(int i) const {return data[i].dJ;}
	double db(int i) const {return data[i].db;}
	double dgtp(int i) const {return data[i].dgtp;}
	double dgtt(int i) const {return data[i].dgtt;}
	double dgpp(int i) const {return data[i].dgpp;}
	double dn(int i) const {return data[i].dn;}
	double ddeltac(int i) const {return data[i].ddeltac;}

	double& dr(int i) {return data[i].dr;}
	double& dvlum(int i) {return data[i].dvlum;}
	double& dforce(int i) {return data[i].dforce;}
	double& dphi(int i) {return data[i].dphi;}
	double& dangmom(int i) {return data[i].dangmom;}
	double& dJ(int i) {return data[i].dJ;}
	double& db(int i) {return data[i].db;}
	double& dgtp(int i) {return data[i].dgtp;}
	double& dgtt(int i) {return data[i].dgtt;}
	double& dgpp(int i) {return data[i].dgpp;}
	double& dn(int i) {return data[i].dn;}
	double& ddeltac(int i) {return data[i].ddeltac;}

	void printError();

	// prints r and whichever variable is specified
	// usage: e.g. for phi, Print(&Galaxy::phi);
	void Print(const_data_access print_ptr) const;

private:

	void integrate(const_data_access ptr_in, data_access ptr_out);

	unsigned n_components;

	// Used to read files
	std::deque< GalaxyRawData > allRawData;

	// Used in calculations
	std::deque< GalaxyData > data;

	std::deque< double > scalings;
};
