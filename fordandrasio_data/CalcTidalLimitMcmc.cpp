#include<cmath>
#include <string.h> 
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector>
#include<unistd.h>

double rand01()
{
  return rand()/static_cast<double>(RAND_MAX);
};

/*
class RvPlanet
{
public:
  double mP, mK, me;
  double mMStar, mSini, mA, mmu;
  double P() const { return mP; };
  double K() const { return mK; };
  double e() const { return me; };
  void SetMstarSinI(const double M, const double sini) 
  { 
    mMStar = M;
    mSini = sini;
    double p = (mP/365.25);
    mmu = 0.;
    mA = pow(mMStar*(1.+mmu)*p*p,1./3.);
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/sini;
    mA = pow(mMStar*(1.+mmu)*p*p,1./3.);
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/sini;
    mA = pow(mMStar*(1.+mmu)*p*p,1./3.);
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/sini;
    mA = pow(mMStar*(1.+mmu)*p*p,1./3.);
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/sini;
  };
  double Mstar() const { return mMStar; };
  double sini() const { return mSini; };
  double a() const { return mA; };
  double mu() const { return mmu; };
};

class TransitPlanet
{
public:
  double mP, mDepth, mDur, mK, me;
  double mMStar, mRStar, mRPlanet, mSini, mA, mmu;
  double P() const { return mP; };
  double Depth() const { return mDepth; };
  double Duration() const { return mDur; };
  double K() const { return mK; };
  double e() const { return me; };
  TransitPlanet(const double p, const double dep, const double dur, const double k, const double e=0.)
    : mP(p), mDepth(dep), mDur(dur), mK(k), me(e) {};
  void SetMRstar(const double M, const double R)
  { 
    mMStar = M;
    mRStar = R;
    mRPlanet = mRStar * sqrt(mDepth);
    double p = (mP/365.25);
    mSini = 1.;
    mmu = 0.;
    mA = pow(mStar*(1.+mmu)*p*p,1./3.);
    double sqrtarg = 1.-R*R/(mA*mA)+M_PI*M_PI*mDur*mDur;
    mSini = (sqrtarg<1.) ? sqrt(sqrtarg) : 1.;
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/mSini;

    mA = pow(mStar*(1.+mmu)*p*p,1./3.);
    sqrtarg = 1.-R*R/(mA*mA)+M_PI*M_PI*mDur*mDur;
    mSini = (sqrtarg<1.) ? sqrt(sqrtarg) : 1.;
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/mSini;

    mA = pow(mStar*(1.+mmu)*p*p,1./3.);
    sqrtarg = 1.-R*R/(mA*mA)+M_PI*M_PI*mDur*mDur;
    mSini = (sqrtarg<1.) ? sqrt(sqrtarg) : 1.;
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/mSini;

    mA = pow(mStar*(1.+mmu)*p*p,1./3.);
    sqrtarg = 1.-R*R/(mA*mA)+M_PI*M_PI*mDur*mDur;
    mSini = (sqrtarg<1.) ? sqrt(sqrtarg) : 1.;
    mmu = mK/(2.*M_PI*mA)*p*sqrt(1.-me*me)/mSini;

  };
  double Mstar() const { return mMStar; };
  double Rstar() const { return mRStar; };
  double RPlanet() const { return mRPlanet; };
  double sini() const { return mSini; };
  double a() const { return mA; };
  double mu() const { return mmu; };
};
*/


class ModelParametersClass
{
public:
  double normalization;
  double alpha, beta, gamma; // Period exponent, Mass exponent, Inner limit exponent
  double xmax, ymax;  // Maximum period, maximum mass
  double xl; // Inner limit constant
  double Kd; // Detectable amplitude (m/s)
  double mo, Po; // Normalizing mass (in Msol) and period (in years)
  double uu, vv; // y_min,a for nominal mass & period & radius; y_min,K for nominal mass & period

  static const double RadiusJup = 0.000477279108;  // in AUs
  static const double MassJup = 0.001; // in M_sol

  ModelParametersClass(const double n, const double a, const double b, const double g, const double maxP, const double maxM, const double pxl = 2., const double kd = 28.4, const double pmo = MassJup, const double pPo = 4./365.25) :
    normalization(n), alpha(a), beta(b), gamma(g), xmax(log(maxP/pPo)), ymax(log(maxM/pmo)), Kd(kd), xl(pxl), mo(pmo), Po(pPo)
  {
    uu = 2./(3.*gamma)*log(Po)-log(mo)-log(RadiusJup/0.462)/gamma-log(xl)/gamma;
    vv = log(Kd/28.4) - log(mo/MassJup)+log(Po)/3.;
    //        std::cerr << "# uu = " << uu << "  vv = " << vv << "\n";
  };
  
  double ProbPlanet (const double P, const double q, const double Rpl, const double MassStar)
  { // q = MassPlanet/MassStar, 
    double x = log(P/Po), y = log(q/mo);
    double w = log(Rpl/RadiusJup), z = log(MassStar);
    double u = uu + z/(3.*gamma) - w/gamma;
    double v = vv - z/3.;
    double epsq =0; // log(1.+muo*exp(ymina));
    double ymina = 2/(3*gamma)*x+u+0*epsq/(3*gamma);
    //    double xmin = 3*gamma/2*(ymax-u-0*epsq);
    double yminK = x/3+v;
    //    double xt = 3*gamma*(u-v+0*epsq)/(gamma-2);

    if(x>xmax)  return 0.;
    if(y>ymax)  return 0.;
    if(y<ymina) return 0.;
    if(y<yminK) return 0.;

    double prob = normalization*exp(-alpha*x-beta*y);
    if(prob<0.)
      prob = 0.;
    if(prob>1.)
      prob = 1.;
    return prob;
  };
  double ProbNotPlanet(const double Rpl, const double MassStar)
  {  
    double w = log(Rpl/RadiusJup), z = log(MassStar);
    double u = uu + z/(3.*gamma) - w/gamma;
    double v = vv - z/3.;

  double cg = -(0.6e1 * exp(-beta * ymax - 0.3e1 / 0.2e1 * gamma * alpha * ymax + 0.3e1 / 0.2e1 * gamma * alpha * u) * beta * alpha + 0.2e1 * exp(-beta * ymax - 0.3e1 / 0.2e1 * gamma * alpha * ymax + 0.3e1 / 0.2e1 * gamma * alpha * u) * beta * beta + 0.3e1 * exp(-(0.3e1 * gamma * alpha * u - 0.3e1 * gamma * alpha * v + beta * u * gamma - 0.2e1 * beta * v) / (gamma - 0.2e1)) * alpha * gamma * beta - 0.6e1 * exp(-(0.3e1 * gamma * alpha * u - 0.3e1 * gamma * alpha * v + beta * u * gamma - 0.2e1 * beta * v) / (gamma - 0.2e1)) * alpha * beta + 0.9e1 * exp(-alpha * xmax - beta * xmax / 0.3e1 - beta * v) * alpha * alpha * gamma + 0.6e1 * exp(-alpha * xmax - beta * xmax / 0.3e1 - beta * v) * alpha * beta - 0.9e1 * exp(-alpha * xmax - beta * ymax) * alpha * alpha * gamma - 0.6e1 * exp(-alpha * xmax - beta * ymax) * alpha * beta - 0.3e1 * exp(-alpha * xmax - beta * ymax) * beta * alpha * gamma - 0.2e1 * exp(-alpha * xmax - beta * ymax) * beta * beta) / alpha / (0.3e1 * alpha * gamma + 0.2e1 * beta) / beta / (0.3e1 * alpha + beta);

    double prob = normalization*cg;
    if(prob<0.)
      {
	//	std::cerr << "Null Integral = " << alpha << ' ' << beta << ' ' << gamma << ' ' << normalization << ' ' << u << ' ' << v << ' ' << xmax << ' ' << ymax << ' ' << prob << '\n';
      prob = 0.;
      }
    if(prob>1.)
      prob = 1.;
    return prob;
  };

};


int main(int argc, char **argv)
{

  const double RJup = 0.000477279108;
  const double mJup = 0.001;
  const double RSun = 0.005;
  pid_t pid = getpid();
  srand(time(NULL)^(pid+pid<<15));
  double RP = 1.2;
  double Kmin = 30.;

  if(argc>1) { sscanf(argv[1],"%lg",&RP); if((RP<0.5)||(RP>2.0)) std::cerr << "# WARNING: Rp = " << RP << " R_J\n"; }
  RP *= RJup;

  
  std::vector<double> msini, P, a, e, K, Mstar, sini;
  std::vector<std::string> name;
  std::ifstream RvData("RvPlanets.txt");
  while(!RvData.eof())
    {
      double Pmin = 0.13, Pmax = 365.25*1;
      double emax = 0.8;
//      double Pmin = 0.13, Pmax = 365.25*1;
//      double emax = 0.8;
      double mmin = 0.1;

      double tmpId, tmpMsini, tmpP, tmpA, tmpE, tmpK, tmpMstar, tmpsini, junk;
      std::string tmpName;
      char c = RvData.peek();
      if(c=='#')
	{
	  RvData.ignore(256,'\n');
	  continue;
	}
      RvData >> junk >> tmpName >> tmpMsini >> tmpP >> tmpA >> tmpE >> tmpK >> tmpMstar;
      RvData.ignore(256,'\n');
      if(tmpMsini<mmin) continue;
      if(tmpE>emax) continue;
      if(tmpK<Kmin) continue;
      if(tmpP>Pmax) continue;
      if(tmpMstar<=0.) 
	{ 
	  std::cerr << "# WARNING: Star " << tmpName << " mstar = " << tmpMstar << "\n"; 
	  continue;
	}
      if(RvData.eof()) break;
      double thissini = -1.;
      if(strcmp(tmpName.c_str(),"hd209458")==0) thissini = 1.;

      tmpP /= 365.25;
      msini.push_back(tmpMsini);
      P.push_back(tmpP);
      a.push_back(tmpA);
      e.push_back(tmpE);
      K.push_back(tmpK);
      Mstar.push_back(tmpMstar);
      sini.push_back(thissini);
      name.push_back(tmpName); 
    }
  RvData.close();
  std::cerr << "# Read " << name.size() << " planets from RvPlanets.txt\n";

  // Information about posterior distribution
  double Pc = 60./365.25, Mc = mJup;
  double MaxP = 365.25, MaxM = 10.*mJup;
  double xmaxP = log(MaxP/Pc), ymax = log(MaxM/Mc);

  int nAlphas = 25, nBetas = 25, nCs = 25, nGammas = 25, nHs = 25;
  int nPlanets = name.size(), nSinIs = 8, nStars = 2000;
  std::vector<double> ProbAlpha(nAlphas,0.), ProbBeta(nBetas,0.), ProbC(nCs,0.);
  std::vector<double> ProbGamma(nGammas,0.), ProbXl(nHs,0.);
  double atanalphamin = -M_PI/2., atanalphamax = M_PI/2., atanbetamin = -M_PI/2., atanbetamax = M_PI/2., cmin = 0.00, cmax = 1.00;
  double atangammamin = -M_PI/2., atangammamax = M_PI/2.;
  //  double gammamin = -1./6., gammamax = -1.;
    //  double gammamin = -1./3., gammamax = -1./3.;
  double xlmin = 0.;
  double xlmax = 0.036/((RP/0.462)*pow(Mc,tan(atangammamin+0.5*(atangammamax-atangammamin)/nGammas)));


  double ProbData = 0.;
  for(int a=0;a<nAlphas;++a)
    {
      //      double alpha = alphamin+(alphamax-alphamin)*(a+0.5)/nAlphas;      
      double alpha = tan(atanalphamin+(atanalphamax-atanalphamin)*(a+0.5)/nAlphas);      
      for(int b=0;b<nBetas;++b)
	{
	  //double beta = betamin+(betamax-betamin)*(b+0.5)/nBetas;
	  double beta = tan(atanbetamin+(atanbetamax-atanbetamin)*(b+0.5)/nBetas);
	  //	  double xl = 2.;
	  for(int g=0;g<nGammas;++g)
	    {
	      //	      double gamma = -1./3.;
	      //	      double gamma = gammamin+(gammamax-gammamin)*(g+0.5)/nGammas;
	      double gamma = tan(atangammamin+(atangammamax-atangammamin)*(g+0.5)/nGammas);

	      for(int h=0;h<nHs;++h)
		{
		  double xl = xlmin+(xlmax-xlmin)*(h+0.5)/nHs;
		  
	  for(int c=0;c<nCs;++c)
	    {
	      double cc = cmin + (cmax-cmin)*(c+0.5)/nCs;

	      ModelParametersClass Model(cc, alpha, beta, gamma, MaxP, MaxM, xl, Kmin, Mc, Pc);
	      double ProbDgM = 1.;
	      double ProbNotOtherPlanets;

	      { // Null detections
		double MstarNull = 1.0;
		ProbNotOtherPlanets = Model.ProbNotPlanet(RP,MstarNull);
		ProbNotOtherPlanets = pow(ProbNotOtherPlanets,nStars);
		ProbDgM *= ProbNotOtherPlanets;
	      } // Null detections
	      if(ProbDgM<=0.) continue;

	      bool UsePlanet;
	      for(int p=0;p<nPlanets;++p)
		{
		  double ProbDpgM = 0., ProbNotOtherPlanets = 0.;

		  UsePlanet = false;
		  for(int j=0;j<nSinIs;++j)
		    {
		      double thissini = sini[p];
		      if(thissini<0.)
			{
			  // double cosi = rand()/static_cast<double>(RAND_MAX);
			  double cosi = 0.+1.*(j+0.5)/nSinIs; 
			  thissini = sqrt(1.-cosi*cosi);
			  if(nSinIs==1) thissini=1.;
			}
		      double q = msini[p]*mJup / (Mstar[p]*thissini);
		      //		  double aR = RP / (0.462 * pow(q*mJup,1./3.));
		      
		      double ProbDpgMI = Model.ProbPlanet(P[p], q, RP, Mstar[p]); 
		      if((ProbDpgMI>=0.)&&(ProbDpgMI<=1.))
			{
			  UsePlanet = true;
			  ProbDpgM +=  ProbDpgMI / nSinIs;
			}

		      /*
		      double x = log(P[p]/Pc);
		      double y = log(q/Mc);

		      std::cerr << P[p]*365.25 << ' ' << q << "  " << x << ' ' << y << "  " << xmaxP << ' ' << xmaxK << ' ' << xminA << ' ' << ymax << "   " << ProbDpgM*nSinIs << ' ' << cg << '\n';
		      */
		      
		    } // sini
		  if(UsePlanet)
		    {
		      ProbDgM *= ProbDpgM * ProbNotOtherPlanets;
		      //		  std::cerr << p << ' ' << ProbDpgM << ' ' << ProbNotOtherPlanets << "  " << ProbDgM << '\n';
		    }
		} // planet
	      if(ProbDgM>=0.)
		{
		  //		  std::cout << " " << alpha << ' ' << beta << ' ' << cc << "  " << gamma << " " << xl << "  " << ProbDgM << ' ' << ProbDgM/ProbNotOtherPlanets << ' ' << ProbNotOtherPlanets << '\n';

		  double pp = ProbDgM / (nAlphas*nBetas*nCs*nGammas*nHs);
		  ProbAlpha[a] += pp;
		  ProbBeta[b] += pp;
		  ProbC[c] += pp;
		  ProbGamma[g] += pp;
		  ProbXl[h] += pp;
		  ProbData += ProbDgM;
		  if(pp>0.)
		    std::cout << alpha << ' ' << beta << ' ' << gamma << ' ' << xl << ' ' << cc << ' ' << pp << '\n';
		}
	    }
		}
	    } // c's
	} // beta
    } // alpha

  

  std::cerr << "# Prob Data = " << ProbData << "\n";
  // Renormalize marginal probability distributions
  for(int i=0;i<nAlphas;++i)
    {
      //      double alpha = alphamin+(i+0.5)*(alphamax-alphamin)/nAlphas;
      ProbAlpha[i] /= ProbData;
    }
  for(int i=0;i<nBetas;++i)
    {
      //      double beta = betamin+(i+0.5)*(betamax-betamin)/nBetas;
      ProbBeta[i] /= ProbData;
    }

  for(int i=0;i<nCs;++i)
    {
      ///      double c = cmin+(i+0.5)*(cmax-cmin)/nCs;
      ProbC[i] /= ProbData;
    }

  for(int i=0;i<nGammas;++i)
    {
      //      double c = cmin+(i+0.5)*(cmax-cmin)/nCs;
      ProbGamma[i] /= ProbData;
    }

  for(int i=0;i<nHs;++i)
    {
      //      double c = cmin+(i+0.5)*(cmax-cmin)/nCs;
      ProbXl[i] /= ProbData;
    }

  std::ofstream alphaout("alpha.out");
  std::cerr << "# Posterior Marginal Distribution for Alpha\n";
  double cum = 0.;
  for(int i=0;i<nAlphas;++i)
    {
      //      double alpha = alphamin+(i+0.5)*(alphamax-alphamin)/nAlphas;
      double alpha = tan(atanalphamin+(i+0.5)*(atanalphamax-atanalphamin)/nAlphas);
      cum += ProbAlpha[i];
      if(ProbAlpha[i]!=0.)
	//	alphaout << alpha << ' ' << ProbAlpha[i] / (alphamax-alphamin)*nAlphas << ' ' << cum << '\n';
	alphaout << alpha << ' ' << ProbAlpha[i] / (atanalphamax-atanalphamin)*nAlphas << ' ' << cum << '\n';
    }
  alphaout.close();

  std::ofstream betaout("beta.out");
  std::cerr << "# Posterior Marginal Distribution for Beta\n";
  cum = 0.;
  for(int i=0;i<nBetas;++i)
    {
      //      double beta = betamin+(i+0.5)*(betamax-betamin)/nBetas;
      double beta = tan(atanbetamin+(i+0.5)*(atanbetamax-atanbetamin)/nBetas);
      cum += ProbBeta[i];
      if(ProbBeta[i]!=0.)
	//	betaout << beta << ' ' << ProbBeta[i] / (betamax-betamin)*nBetas << ' ' << cum << '\n';
	betaout << beta << ' ' << ProbBeta[i] / (atanbetamax-atanbetamin)*nBetas << ' ' << cum << '\n';
    }
  betaout.close();


  std::ofstream constout("const.out");
  std::cerr << "# Posterior Marginal Distribution for C\n";
  cum = 0.;
  for(int i=0;i<nCs;++i)
    {
      double c = cmin+(i+0.5)*(cmax-cmin)/nCs;
      cum += ProbC[i];
      if(ProbC[i]!=0.)
        constout << c << ' ' << ProbC[i] / (cmax-cmin)*nCs << ' ' << cum << '\n'
;
    }
  constout.close();

  std::ofstream gammaout("gamma.out");
  std::cerr << "# Posterior Marginal Distribution for gamma\n";
  cum = 0.;
  for(int i=0;i<nGammas;++i)
    {
      //      double gamma = gammamin+(i+0.5)*(gammamax-gammamin)/nGammas;
      double gamma = atangammamin+(i+0.5)*(atangammamax-atangammamin)/nGammas;
      cum += ProbGamma[i];
      if(ProbGamma[i]!=0.)
	//	gammaout << gamma << ' ' << ProbGamma[i] / fabs(gammamax-gammamin)*nGammas << ' ' << cum << '\n';
	gammaout << gamma << ' ' << ProbGamma[i] / fabs(atangammamax-atangammamin)*nGammas << ' ' << cum << '\n';
    }
  gammaout.close();

  std::ofstream xlout("xl.out");
  std::cerr << "# Posterior Marginal Distribution for xl\n";
  cum = 0.;
  for(int i=0;i<nHs;++i)
    {
      double xl = xlmin+(i+0.5)*(xlmax-xlmin)/nHs;
      cum += ProbXl[i];
      if(ProbXl[i]!=0.)
	xlout << xl << ' ' << ProbXl[i] / (xlmax-xlmin)*nHs << ' ' << cum << '\n';
    }
  xlout.close();

  return 0;
}
