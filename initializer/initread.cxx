#include "initread.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <limits>

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right

inline std::string trim(std::string s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

using namespace std;

Parameters::Parameters(string pName) :
	m_omega_cdm(-1.0),
	m_deut(-1.0),
	m_omega_nu(-1.0),
	m_hubble(-1.0),
	m_ss8(-1.0),
	m_ns(-1.0),
	m_w_de(-1.0),
	m_wa_de(-1.0),
	m_Tcmb(-1.0),
	m_Zdec(-1.0),
	m_neff_massless(-1.0),
	m_neff_massive(-1.0),
	m_cdm_mass(-1.0),
	zpr_size(0),
	m_zin(-1.0),
	m_zfin(-1.0),
	m_new_ics(1),
	m_debug(0),
	m_trans(0),
	m_iseed(-1),
	m_alpha(1.0),
	m_ng(-1),
	m_rL(-1.0),
	m_nsteps(-1),
    m_omega_baryon(-1.0),
    m_omega_cb(-1.0),
    m_omega_matter(-1.0),
	m_omega_radiation(-1.0),
    m_f_nu_massless(-1.0),
    m_f_nu_massive(-1.0),
	m_pknbin1d(0),
	m_pklogbins(false),
	m_pknmode3d(0)
	{
		readNewParams(pName);
	}

void Parameters::readNewParams(string &fn) {
	m_params.clear();

	std::stringstream myfile;
	if (!getRank0Stream(fn.c_str(), myfile)) {
		std::ostringstream ost;
		ost << "initread cannot open '"
			<< fn
			<< "'";
		throw std::runtime_error( ost.str() );
	}

	while (myfile.good()){
		string line;
		getline(myfile, line);
		size_t pos = line.find_first_not_of(" \t");
		if(pos < line.size()){
			line = line.substr(pos);
		}
		
		if(line.empty()) {
			continue;
		} else if (line[0]== '#'){
			continue;
		}

		pos = line.find_first_of(" \t");
		string key = line.substr(0, pos);

		if(pos < line.size()){
			line = line.substr(pos);
		} else {
			line = "";
		}

		pos = line.find_first_not_of(" \t");
		if(pos < line.size()){
			line = line.substr(pos);
		} else {
			line = "";
		}

		string value = line;

		m_params.insert(std::make_pair(key, value));
	}
	
	m_omega_cdm = atof(m_params["OMEGA_CDM"].c_str());
	m_deut = atof(m_params["DEUT"].c_str());
	m_omega_nu = atof(m_params["OMEGA_NU"].c_str());
	m_hubble = atof(m_params["HUBBLE"].c_str());
	m_ss8 = atof(m_params["SS8"].c_str());
	m_ns = atof(m_params["NS"].c_str());
	m_w_de = atof(m_params["W_DE"].c_str());
	m_wa_de = atof(m_params["WA_DE"].c_str());
	m_Tcmb = atof(m_params["T_CMB"].c_str());
	m_Zdec = atof(m_params["Z_DEC"].c_str());
	m_neff_massless = atof(m_params["N_EFF_MASSLESS"].c_str());
	m_neff_massive = atof(m_params["N_EFF_MASSIVE"].c_str());
	
	std::istringstream ss(m_params["CDM_MASS"].c_str());
	ss >> m_cdm_mass; //in scientific notation

	ss.clear();

	ss.str(m_params["SCATTER"].c_str());
	ss >> m_scatter;
	
	ss.clear();

	m_zin = atof(m_params["Z_IN"].c_str());
	m_new_ics = atoi(m_params["NEW_ICS"].c_str());
	m_debug = atoi(m_params["DEBUG"].c_str());
	if(m_params["TRANS"] == "CMB"){
		m_trans = 0;
	} else if (m_params["TRANS"] == "KH"){
		m_trans = 1;
	} else if (m_params["TRANS"] == "AXCAMB"){
		m_trans = 2;
	}
	m_iseed = atoi(m_params["ISEED"].c_str());
	
	if(m_params.count("ALPHA")){
		m_alpha = atof(m_params["ALPHA"].c_str());
	}
	if(m_params.count("PK_NBIN_1D")){
		m_pknbin1d = atoi(m_params["PK_NBIN_1D"].c_str());
	}
	if(m_params.count("PK_NMODE_3D")){
		m_pknbin1d = atoi(m_params["PK_NMODE_3D"].c_str());
	}

	istringstream(m_params["PK_LOGBINS"].c_str()) >> std::boolalpha >> m_pklogbins;
	
	m_ng = atoi(m_params["NG"].c_str());
	m_rL = atof(m_params["RL"].c_str());
	m_zfin = atof(m_params["Z_FIN"].c_str());
	m_nsteps = atoi(m_params["NSTEPS"].c_str());

	std::istringstream line(m_params["KMODE"].c_str());
	for(int i=0; i<3; i++){
		int mode = 0;
		line >> mode;
		m_kmode[i] =mode;
	}

	std::istringstream instream1(m_params["Z_PR"].c_str());
	while(!instream1.eof()){
		REAL z;
		if(instream1 >> z){
			zpr_size++;
		}
	}
	m_zpr = (REAL *)malloc(sizeof(REAL)*zpr_size);
	m_apr = (REAL *)malloc(sizeof(REAL)*zpr_size);
	std::istringstream instream(m_params["Z_PR"].c_str());
	for(int i=0; i<zpr_size; i++){
		REAL z;
		instream >> z;
		m_zpr[i]=z;
	}

	for(int i=0; i<zpr_size; i++){
		m_apr[i] = 1.0/(1.0 + m_zpr[i]);
	}

	m_inTransfer = m_params["TRANSFER0"];
	m_initTrans = m_params["TRANSFERIN"];

	m_qho_freq = atof(m_params["FREQUENCY"].c_str());
	m_width = atof(m_params["GAUSSIAN_WIDTH"].c_str());
	m_shift = atof(m_params["INITIAL_SHIFT"].c_str());
	m_duration = atof(m_params["DURATION"].c_str());

	ss.str(m_params["OUTBASE"].c_str());
	m_outBase  = trim(ss.str()); 
	ss.clear();

	m_omega_baryon = m_deut/ m_hubble / m_hubble;
    m_omega_cb = m_omega_cdm + m_omega_baryon;
    m_omega_matter = m_omega_cb + m_omega_nu;
	m_omega_radiation = 2.471e-5*pow(m_Tcmb/2.725f,4.0f)/pow(m_hubble,2.0f);
    m_f_nu_massless = m_neff_massless*7.0/8.0*pow(4.0/11.0,4.0/3.0);
    m_f_nu_massive = m_neff_massive*7.0/8.0*pow(4.0/11.0,4.0/3.0);

	//TimeStepper variables
	m_ain = 1.0 / (1 + m_zin);
	m_afin = 1.0 / (1 + m_zfin);
	m_pp = pow(m_ain, m_alpha);
	m_pfin = pow(m_afin, m_alpha);

}

string Parameters::getParamsString(){
	std::stringstream ss;
	ss << "TRANS ";
	switch(m_trans) {
		case 0: 
			ss <<"CMB\n";
			break;
		case 1: 
			ss <<"KH\n";
			break;
		case 2: 
			ss <<"AxionCAMB\n";
			break;
		default:
			ss << "CMB\n";
			break;
	}
    ss << "TRANSFER0 " << m_inTransfer << "\n";
    ss << "TRANSFER IN " << m_initTrans << "\n";
    ss << "OUTBASE " << m_outBase << "\n";
	ss << "\n";

    ss << "NG " << m_ng << "\n";
    ss << "NS " << m_ns << "\n";
    ss << "I_SEED " << m_iseed << "\n";
    ss << "Z_IN " << m_zin << "\n";
    ss << "Z_FIN " << m_zfin << "\n";
    ss << "N_STEPS " << m_nsteps << "\n";
    ss << "CDM_MASS " << m_cdm_mass << "\n";
    ss << "SCATTER " << m_scatter << "\n";
	ss << "HUBBLE " << m_hubble << "\n";
    ss << "OMEGA_CDM " << m_omega_cdm << "\n";
    ss << "DEUT " << m_deut << "\n";
    ss << "RL " << m_rL << "\n";
    ss << "SS8 " << m_ss8 << "\n";
    ss << "W_DE " << m_w_de << "\n";
    ss << "WA_DE " << m_wa_de << "\n";
    ss << "OMEGA_NU " << m_omega_nu << "\n";
    ss << "T_CMB " << m_Tcmb << "\n";
    ss << "N_EFF_MASSLESS " << m_neff_massless << "\n";
    ss << "N_EFF_MASSIVE " << m_neff_massive << "\n";
	ss << "KMODE ("<<m_kmode[0]<<","<<m_kmode[1]<<","<<m_kmode[2]<<")"<< "\n";
	ss << "Z_PR ";
	for(int i=0; i<zpr_size; i++){
		ss<<m_zpr[i];
		ss<<" ";
	}
	ss<<"\n";

    return ss.str();
}
