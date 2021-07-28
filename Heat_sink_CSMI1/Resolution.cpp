#include<iostream>
#include<math.h>
#include<iomanip>
#include<fstream>
#include<vector>
#include <map>
#include<iomanip>
#include"randomnumber.hpp"

class Geometry
{
public:
	//constructeur par défaut
	Geometry()
	:
	M_Lx(0.04),M_Ly(0.004),M_Lz(0.05){}

	//constructeur à 3 arguments
	Geometry(double Lx,double Ly,double Lz)
	:
	M_Lx(Lx),M_Ly(Ly),M_Lz(Lz){}
	
	//destructeur
	virtual ~Geometry() {}

	//accesseur
	double getLx() const {return M_Lx;}
	double getLy() const {return M_Ly;}
	double getLz() const {return M_Lz;}
	
	//mutateur
	double setLx(double Lx){return M_Lx = Lx;}
	double setLy(double Ly){return M_Ly = Ly;}
	double setLz(double Lz){return M_Lz = Lz;}

	//p : perimetre transversale
	double p() const
	{
		return (M_Ly +M_Lz)/2;
	}

	//S : aire transversale
	double S() const
	{
		return (M_Ly * M_Lz);
	}
	
	// retourne l'objet lui meme
	Geometry getG() const {return *this;}

protected:
	//attributs
	double M_Lx,M_Ly,M_Lz;
};
//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Geometry const& c)
{
	return o << " Lx "<< c.getLx() <<  " Ly "<< c.getLy() <<  " Lz "<< c.getLz() ;
}

//surcharge entrée

std::istream& operator>> (std::istream& i, Geometry & c)
{
	char t[4]; double Lx, Ly, Lz;
	i >> t >> Lx >> t >> Ly>> t >> Lz;// Lx f Ly f Lz f
	c.setLx(Lx);
	c.setLy(Ly);
	c.setLz(Lz);
	return i;
}

//heritage de la source de production de la classe Geometry
class Convection : public Geometry
{
public:
	//constructeur par defaut
	Convection()
	:
	Geometry(),M_hc(200),M_Te(293)
	{}
	
	//constructeur à 3 arguments
	Convection(Geometry const& G, double const& hc,double const& Te )
	:
	Geometry(G),M_hc(hc),M_Te(Te)
	{}

	//constructeur à 2 arguments 
	//car la temperature ambiante est toujours constante
	Convection(Geometry const& G, double const& hc)
	:
	Geometry(G),M_hc(hc),M_Te(293)
	{}

	//destructeur
	virtual ~Convection() {}

	//accesseur
	double gethc() const {return M_hc;}
	double getTe() const {return M_Te;}
	//Geometry getG() const {return M_G;}
	
	//mutateur
	double sethc(double hc){return M_hc = hc;}
	double setTe(double Te){return M_Te = Te;}
	

protected:
	double M_hc,M_Te;

};
//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Convection const& c)
{
	return o << " hc "<< c.gethc() << " Te "<< c.getTe();
}

class Conduction
{
public:
	//constructeur par défaut
	Conduction()
	:
	M_K(164){}

	//constructeur à 1 arguments
	Conduction(double K)
	:
	M_K(K){}

	//accesseur
	double getK() const {return M_K;}

	//mutateur
	double setK(double K){return M_K = K;}
	
	//destructeur
	virtual ~Conduction() {}


protected:
	//attributs
	double M_K;
};
//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Conduction const& c)
{
	return o << " K = "<< c.getK();
}

class Inertie
{
public:
	//constructeur par défaut
	Inertie()
	:
	M_rho(2700),M_Cp(940),M_N(600),M_tfinal(300)
	{
		M_tstep = M_tfinal/M_N;
		M_coeff = M_rho*M_Cp/M_tstep;
	}

	//constructeur à 2 arguments
	Inertie(double rho, double Cp)
	:
	M_rho(rho),M_Cp(Cp)
	{
		M_coeff = M_rho*M_Cp/M_tstep;
		M_tstep = M_tfinal/M_N;
	}

	//accesseur
	double getrho() const {return M_rho;}
	double getCp() const {return M_Cp;}
	double gettstep() const {return M_tstep;}
	double gettfinal() const {return M_tfinal;}
	double getcoeff() const {return M_coeff;}
	int getN() const {return M_N;}

	//mutateur
	double setrho(double rho){return M_rho = rho;}
	double setCp(double Cp){return M_Cp = Cp;}
	
	//destructeur
	virtual ~Inertie() {}

protected:
	//attributs
	double M_rho,M_Cp,M_tfinal,M_tstep,M_coeff;
	int M_N;
};
//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Inertie const& c)
{
	return o << " rho = "<< c.getrho()<< " Cp = "<< c.getCp();
}

class Parametres : public Inertie, public Conduction, public Convection
{
public:
	//constructeur par defaut
	Parametres() 
	: 
	M_M(100000),M_stat(1){M_h =M_Lx/M_M;}

	//construceteur à 2 argument
	Parametres(int const&  M,bool const&  stat) 
	: 
	M_M(M),M_h(M_Lx/M),M_stat(stat){}
	
	//construceteur à 5 argument
	Parametres(Inertie const& I,Conduction const& C,Convection const& V, int const& M,bool const& stat)
	: 
	Inertie(I),Conduction(C),Convection (V),M_M(M),M_h(M_Lx/M),M_stat(stat){}

	//destructeur
	virtual ~Parametres() {}

	//accesseur
	int getM() const {return M_M;}
	double geth() const {return M_h;}
	bool getstat() const {return M_stat;}

	//mutateur
	int setM(int M)
	{
		M_h = M_Lx/M;
		return M_M = M;
	}

	bool setstat(bool stat)
	{
		return M_stat = stat;
	}

	typedef std::vector<double> vec_type;

	//construire le vecteur a de M elements
	vec_type vec_a() const
	{
		vec_type v;

		// remplir a[0]...a[M-2]
		for(int i=0;i<M_M-1;i++) v.push_back(-M_K/(M_h*M_h));
		// remplir a[M-1]
		v.push_back(M_K/M_h);
		return v;	
	}
	
	//construire le vecteur b de M+1 elements
	vec_type vec_b() const
	{
		vec_type v;
		int option = (M_stat+1)%2;
		// remplir b[0]
		v.push_back(M_K/M_h);
		// remplir b[1]...b[M-1]
		for(int i=1;i<M_M;i++)
		{ 
			v.push_back(2*M_K/(M_h*M_h)+this->p()*M_hc/this->S() + option * this->M_coeff);
		}
		// remplir b[M]
		v.push_back(-M_K/M_h);
		return v;	
	}
	
	//construire le vecteur c de M elements
	vec_type vec_c() const
	{
		vec_type v;
		// remplir c[0]
		v.push_back(-M_K/M_h);
		// remplir c[1]...c[M-1]
		for(int i=1;i<M_M;i++) v.push_back(-M_K/(M_h*M_h));

		return v;	
	}
	
	
protected:
	int M_M;
	double M_h;
	bool M_stat;
};

//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Parametres const& c)
{
	return o << " M = "<< c.getM() << " h= " << c.geth();
}

class Matrice : virtual public Parametres
{
public:
	//constructeur par défaut
	Matrice()
	:
	Parametres(){}

	//constructeur 
	Matrice(Parametres const& P)
	:
	Parametres(P){}

	//Methode LU
	typedef std::vector<double> vec_type;
	typedef std::pair<vec_type,vec_type> paire_vec;
	//construire le vecteur b* et c* de M+1 et M elements
	paire_vec bc_LU() const
	{
		vec_type L,U;
		vec_type a = this->vec_a();
		vec_type b = this->vec_b();
		vec_type c= this->vec_c();
		// remplir b*[0] et c*[0]
		L.push_back(b[0]);
		U.push_back(c[0]/L[0]);
		// remplir b*[1]...b*[M-1] et c*[1]...c*[M-1]
		for(int i=1;i<M_M;i++)
		{ 
			L.push_back(b[i]- a[i-1]*U[i-1]);
			U.push_back(c[i]/L[i]);
		}
		// remplir b*[M]
		L.push_back(b[M_M]- a[M_M-1]*U[M_M-1]);

		//retourner un paire (b*,c*)
		paire_vec p = std::make_pair(L, U);
		return p;	
	}
	
};	

class Source : virtual public Parametres
{
public:
	//constructeur par défaut
	Source()
	:
	Parametres(),M_phi(1.25*std::pow(10,5)){}

	//constructeur à 1 argument
	Source(Parametres const& P)
	:
	Parametres(P),M_phi(1.25*std::pow(10,5)){}

	//constructeur à 2 arguments
	Source(Parametres const& P,double const& phi)
	:
	Parametres(P),M_phi(phi){}

	//accesseur
	double getphi() const {return M_phi;}

	//mutateur
	double setphi(double const& phi ){return M_phi = phi;}
	
	typedef std::vector<double> vec_type;

	//construire le vecteur F de M+1 elements
	vec_type F() const
	{
		vec_type v;
		// remplir F[0]
		v.push_back(M_phi);
		// remplir F[1]...F[M-1]
		for(int i=1;i<M_M;i++) v.push_back((M_hc*this->p()*M_Te)/this->S());
		// remplir F[M]
		v.push_back(0);
		return v;	
	}

protected:
	//attributs
	double M_phi;
};

class Resolution : public Source, public Matrice
{
public:
	//constructeur par defaut
	Resolution()
	:
	Source(),Matrice(){M_now = this->M_tfinal;}

	//constructeur à 2 argument
	Resolution(Source const& S,Matrice const& M)
	:
	Source(S),Matrice(M){M_now = this->M_tfinal;}

	//accesseur
	double getnow() const {return M_now;}

	//mutateur
	double setnow(double const& now )
	{
		if(now <= M_tfinal)
			M_now = now;
		return M_now;
	}

	typedef std::vector<double> vec_type;

	//to celicius
	vec_type To_Celcius(vec_type X) const
	{
		for(int i=0;i<this->M_M+1;i++) 
		{
			X[i] = X[i] - 273;
		}
		return X;
	}
	//Resolution
	vec_type Solve_T() const
	{
		//initialisation
		vec_type Y;
		vec_type X(this->M_M+1);
		vec_type  F_i = this->F();
		vec_type a = this->vec_a();
		vec_type b_etoile = this->bc_LU().first;
		vec_type c_etoile = this->bc_LU().second;
		
		if(M_stat == 1) //cas static
		{
			
			//resolution LY = F
			// remplir Y[0]
			Y.push_back(F_i[0]/b_etoile[0]);
			// remplir Y[1]...Y[M]
			for(int i=1;i<this->M_M+1;i++) 
			{
				Y.push_back((F_i[i]- a[i-1]*Y[i-1])/b_etoile[i]);
			}

			//resolution UX = Y
			// remplir X[M]
			X[this->M_M]=Y[this->M_M];
			// remplir X[1]...X[M]
			for(int i=this->M_M-1;i>-1;i--) 
			{
				X[i]=Y[i]- c_etoile[i]*X[i+1];
			}
			X = this->To_Celcius(X); //conversion °C
			return X; 	
		}
		else //cas dynamic
		{
		
			//intialiser T⁰
			vec_type Tinit(M_M +1,M_Te);
			Tinit[0]=0;
			Tinit[M_M]=0;
			vec_type F_dynamic(M_M +1);
			F_dynamic[0] = F_i[0];
			F_dynamic[M_M] = F_i[M_M];
			std::cout << "---------------------Parametres dynamic-----------------" << std::endl;
			std::cout << " Le coefficient de l'inértie =" << this->M_coeff << std::endl;
			std::cout << " Discritisation temporelle = " << this->M_N << std::endl;
			std::cout << " L'instant de l'evaluation final =" << this->M_tfinal << std::endl;
			std::cout << "----------------------------------------------------" << std::endl;
			//boucle sur le temps
			int n = 0;
			double t=0;
			while(t <= M_now )
			{
				//resolution LY = F_dynamic
				//définir F' = F + coeff* T^n
				
				for(int i =1;i<M_M;i++)
				{
					F_dynamic[i] = F_i[i]+this->M_coeff*Tinit[i];
				}
					
				// remplir Y[0]
				Y.push_back(F_dynamic[0]/b_etoile[0]);
				// remplir Y[1]...Y[M]
				for(int i=1;i<this->M_M+1;i++) 
				{
					Y.push_back((F_dynamic[i]- a[i-1]*Y[i-1])/b_etoile[i]);
				}

				//resolution UX = Y
				// remplir X[M]
				X[this->M_M]=Y[this->M_M];
				// remplir X[1]...X[M]
				for(int i=this->M_M-1;i>-1;i--) 
				{
					//std::cout << " i " << i << std::endl;
					X[i]=Y[i]- c_etoile[i]*X[i+1];
				}
				//T^n <-- T^n+1
				for(int i=1;i<M_M;i++) 
				{
					Tinit[i] = X[i];
				}
					
				//l'instant suivant
				n++;
				t = n*M_tstep;
				//std::cout << " --- instant en (s) ----" << t << std::endl;
			}
			std::cout << "--------------Résultat à l'instant : " << M_now << " sec---------------" << std::endl;
			X = this->To_Celcius(X); //conversion °C
			return X; 	
		}
		
	}

	
	//Solution exacte
	paire_vec Exacte() const
	{
		vec_type X,Texact;
		double a = sqrt(M_hc*this->p()/(M_K*this->S()));
		for(int i=0;i<M_M+1;i++) 
		{
			X.push_back(i*M_h);
			Texact.push_back(M_Te + (M_phi*cosh(a*(M_Lx-X[i])))/(M_K*a*sinh(a*M_Lx)));
		}
		Texact = this->To_Celcius(Texact); //conversion °C	

		//retourner un paire (x,Texact)
		paire_vec p = std::make_pair(X, Texact);
		return p;
	}

	//export xi,Tsolve, Texact
	void Export()
	{
		int nb_discrit = this->M_M+1;//nombre de points de discritisation
      		std::ofstream ofile("Statique_output.csv");
      		//ofile << nb_discrit<< std::endl; 
		ofile << " xi " <<","<< " Tsolution " <<","<<  " Texacte " << "\n";
		paire_vec p = this-> Exacte();
		vec_type T = this-> Solve_T();
		for (int i=0;i<nb_discrit;++i)
		{
			ofile  << std::scientific<< p.first[i]<<","<<  T[i] <<","<< p.second[i] << "\n";
		}
		std::cout << " ------------------------------------------------- " << std::endl;
		std::cout << " ----------Résultat stationnaire exporté---------- " << std::endl;
		std::cout << " ------------------------------------------------- " << std::endl;
		ofile.close();
	}
private:
	double M_now;
		
};
//surcharge sortie
	
std::ostream& operator<< (std::ostream& o, Resolution const& c)
{
	return o <<" Geometry: "<< " Lx "<< c.getLx() <<  " Ly "<< c.getLy() <<  " Lz "<< c.getLz() << "\n" << " Aire: " << c.S() << " Perimetres: " << c.p() << "\n Inertie : rho "<<c.getrho() << " Cp " << c.getCp() << "\n" << " Conduction est : K= " << c.getK() << "\n" << " Convection : Te " << c.getTe() << " hc " << c.gethc() << "\n" << " Flux source :" << c.getphi() << "\n le nombre de point de discritisation : " << c.getM();
}

//surcharge entrée

std::istream& operator>> (std::istream& i, Resolution & R)
{
	char t[11]; 
	bool stat;
	int M;
	double Lx, Ly, Lz,phi,hc,Te;
	
	i >> t >> Lx >> t >> Ly >> t >> Lz;  
	i >> t >> M >> t >> phi;
	i >> t >> hc >> t >> Te;
	i >> t >> stat;  
	R.setLx(Lx);
	R.setLy(Ly);
	R.setLz(Lz);
	R.setM(M);
	R.setphi(phi);
	R.sethc(hc);
	R.setTe(Te);
	R.setstat(stat);
	return i;
}



int main()
{
						//\\
					       //!!\\
					      //----\\
///------les tests unitaires éfféctués lors de la verification sont commentée en bas du code-----//

//--------------------------------------//Validation//--------------------------------------//

	//lire les parametres de configuration
	Resolution Ailette;
	std::ifstream rfile("file.cfg",std::ios::in);
	if(rfile)
	{
		rfile >> Ailette; 
	}
	rfile.close();
	typedef std::vector<double> vec_type;
	typedef std::pair<vec_type,vec_type> paire_vec;
	vec_type X = Ailette.Solve_T();
	paire_vec p;
	p = Ailette.Exacte(); 
	int M = Ailette.getM();
	if(Ailette.getstat() == 0)
	{
		std::cout << " -------Resolution Dynamique--------- " << std::endl; 
		std::cout << "\n Voila les parametres de la simulation: \n\n" << Ailette << std::endl;

		Ailette.setnow(15);//instant de l'evaluation 	
		std::cout << " ------------------------------------------------- " << std::endl;
		std::cout << " La solution T(0):" << X[0]<< " °C"<<std::endl;
		std::cout << " La solution T(M/2):" << X[M/2]<< " °C"<< std::endl;
		std::cout << " La solution T(M):" << X[M] << " °C"<< std::endl;
		std::cout << " ------------------------------------------------- " << std::endl;
	}
	else
	{
		std::cout << " ---------------Resolution statique----------- " << std::endl;
		std::cout << "\n Voila les parametres de la simulation: \n\n" << Ailette << std::endl;
	
		std::cout << " ------------------------------------------------- " << std::endl;
		std::cout << " La solution T(0):" << X[0]<< " °C"<<std::endl;
		std::cout << " La solution T(M/2):" << X[M/2]<< " °C"<< std::endl;
		std::cout << " La solution T(M):" << X[M] << " °C"<< std::endl;
		std::cout << " ------------------------------------------------- " << std::endl;

		vec_type xi = p.first;
		vec_type Ti = p.second;
		std::cout << " -------------------Solution exacte-------------------- " << std::endl;
		std::cout << " La solution exacte T(0):" << Ti[0]<< " °C"<<std::endl;
		std::cout << " La solution exacte T(M/2):" << Ti[M/2]<< " °C"<< std::endl;
		std::cout << " La solution exacte T(M):" << Ti[M] << " °C"<< std::endl;
		std::cout << " -------------------Exportation-------------------- " << std::endl;
		Ailette.Export();
		std::cout << " ------------------------------------------------- " << std::endl;
		
	}
/*
//-----------------------------//tests unitaires//----------------------------------//	
	
	std::cout << " ------------------------- " << std::endl;
	Geometry G(1,1,1);
	std::cout << " ailette test est :" << G << std::endl;

	std::cout << " ------------------------- " << std::endl;
	//lire les parametres de configuration
	std::ifstream rfile("file.cfg",std::ios::in);

	Geometry L;
	if(rfile)
	{
		rfile >> L; 
	}
	rfile.close();
	Geometry K(G);//par copie
	std::cout << " ailette étudiée est :" << K << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " Le perimetre est :" << K.p() << std::endl;
	std::cout << " L'aire est :" << L.S() << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Convection V;
	std::cout << " ventilateur par défaut est :" << V<< std::endl;
	std::cout << " la geometrie est :" << V.getG()<< std::endl;
	std::cout << " lx est :" << V.getLx()<< std::endl;
	std::cout << " Le perimetre est :" << V.p() << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Convection V1(G,202);
	std::cout << " Le perimetre test est :" << V1.p() << std::endl;
	V1.setLy(10);
	std::cout << " Le perimetre test est :" << V1.p() << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Conduction C;
	std::cout << " La conduction:" << C << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Inertie I;
	std::cout << " L'inertie:" << I << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------Parametres----------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Parametres P;
	P.setM(110);
	P.setLx(0.14);//modifier les parametres
	std::cout << " Les Parametres sont:" << P << std::endl;
	std::cout << " ------------------------- " << std::endl;
	int M = 100;	
	Parametres newP(I,C,V1,M,0);
	std::cout << " Les Parametres sont:" << newP << std::endl;
	std::cout << " Le vecteur a(0):" << newP.vec_a()[0]<< std::endl;
	std::cout << " Le vecteur a(M-1):" << newP.vec_a()[M-1]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " Le vecteur b(0):" << newP.vec_b()[0]<< std::endl;
	std::cout << " Le vecteur b(M-1):" << newP.vec_b()[M-1]<< std::endl;
	std::cout << " Le vecteur b(M):" << newP.vec_b()[M]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " Le vecteur c(0):" << newP.vec_c()[0]<< std::endl;
	std::cout << " Le vecteur c(M-1):" << newP.vec_c()[M-1]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------Matrice Stationnaire---------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Matrice Mat(newP);
	std::cout << " ------------------------- " << std::endl;
	std::cout << " Le vecteur b*(0):" << Mat.bc_LU().first[0]<< std::endl;
	std::cout << " Le vecteur b*(M-1):" << Mat.bc_LU().first[Mat.getM()-1]<< std::endl;
	std::cout << " Le vecteur b*(M):" << Mat.bc_LU().first[Mat.getM()]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " Le vecteur c*(0):" << Mat.bc_LU().second[0]<< std::endl;
	std::cout << " Le vecteur c*(M-1):" << Mat.bc_LU().second[Mat.getM()-1]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------Source---------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Source S(newP,1.25*std::pow(10,5));
	std::cout << " Le flux source est:" << S.getphi()<< std::endl;
	std::cout << " Le vecteur F(0):" << S.F()[0]<< std::endl;
	std::cout << " Le vecteur F(M-1):" << S.F()[S.getM()-1]<< std::endl;
	std::cout << " Le vecteur F(M):" << S.F()[S.getM()]<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ------------Inertie---------- " << std::endl;
	std::cout << " ------------------------- " << std::endl;
	Inertie Id;
	std::cout << " Le pas de temps:" << Id.gettstep()<< std::endl;
	std::cout << " Le coefficient de l'inertie:" << Id.getcoeff()<< std::endl;
	Id.settstep(1);
	std::cout << " Le pas de temps:" << Id.gettstep()<< std::endl;
	std::cout << " Le coefficient de l'inertie aprés modif du pas:" << Id.getcoeff()<< std::endl;
	std::cout << " ------------------------- " << std::endl;
	std::cout << " ----------Parametres: stat-------- " << std::endl;
	bool stat = 0;
	int coef = (stat+1)%2;
	std::cout << " stat = " << stat << std::endl;
	std::cout << " stat coef = " << coef << std::endl;
	std::cout << " ------------------------- " << std::endl;
*/

//--------------------------------------//verification//--------------------------------------//

/*	
	std::cout << " ------------------------------------------------- " << std::endl;
	std::cout << " -----------Resolution cas Stationnaire----------- " << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	Resolution R;
	int M= 10000;
	R.setM(M);
	//R.sethc(10);
	typedef std::vector<double> vec_type;
	vec_type X = R.Solve_T();
	//parametres
	std::cout << "\n Voila les parametres de la simulation: \n\n" << R << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	std::cout << " La solution T(0):" << X[0]<< " °C"<<std::endl;
	std::cout << " La solution T(M-1):" << X[M-1]<< " °C"<< std::endl;
	std::cout << " La solution T(M):" << X[M] << " °C"<< std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	std::cout << " -----------Verification cas Stationnaire----------- " << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	typedef std::pair<vec_type,vec_type> paire_vec;
	paire_vec p;
	p = R.Exacte();
	vec_type x = p.first;
	vec_type T = p.second;
	std::cout << " La solution T(0):" << T[0]<< " °C"<<std::endl;
	std::cout << " La solution T(M-1):" << T[M-1]<< " °C"<< std::endl;
	std::cout << " La solution T(M):" << T[M] << " °C"<< std::endl;
	std::cout << " Pour un x allant de " << x[0] << " à " << x[M] << " avec un pas de "<< R.geth() <<std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	R.Export();//export
	std::cout << " ------------------------------------------------- " << std::endl;
	std::cout << " -----------Verification cas intationnaire----------- " << std::endl;
	std::cout << " ------------------------------------------------- " << std::endl;
	Resolution Rd;
	int Md= 50000;
	Rd.setM(Md);
	Rd.setstat(0);
	//R.sethc(10);
	vec_type Xd = Rd.Solve_T();
	//parametres
	std::cout << "\n Voila les parametres de la simulation: \n\n" << R << std::endl;
	std::cout << " -----------------------t=300s--------------------------: " << Rd.gettfinal()<< std::endl;
	std::cout << " La solution T(0):" << Xd[0]<< " °C"<<std::endl;
	std::cout << " La solution T(M-1):" << Xd[M-1]<< " °C"<< std::endl;
	std::cout << " La solution T(M):" << Xd[M] << " °C"<< std::endl;
	std::cout << " -----------------------t=15s-------------------------- " << std::endl;
	Rd.setnow(5);
	Xd = Rd.Solve_T();
	std::cout << " l'instant du résultat est "<< Rd.getnow() << std::endl;
	std::cout << " 0: dynamic et 1: static ? "<< Rd.getstat() << std::endl;
	std::cout << " La solution T(0):" << Xd[0]<< " °C"<<std::endl;
	std::cout << " La solution T(M-1):" << Xd[M-1]<< " °C"<< std::endl;
	std::cout << " La solution T(M):" << Xd[M] << " °C"<< std::endl;
	std::cout << " ------------------------------------------------------ " << std::endl;
*/
 
	return 0;
}









