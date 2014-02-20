///////////////////////////////////////////////////////////////////
//SOFT BREHM / ELASTIC KINEMATICS
//Denner & Pozzorini, 1999
//Electromagnetic Section of their paper, includes Z bosons though
///////////////////////////////////////////////////////////////////



std::complex<double> se_c (se,0); //Need this for some calculations




///////////////////////////////////////////////////////////////////
//Below, the 'I' are relations for box diagrams
///////////////////////////////////////////////////////////////////

double epsilon = 0.;
std::complex<double> ieps (0,epsilon);

std::complex<double> I5ggST(double x)
    {
    return alpha/twopi*( se/(2.*(se+te(x))) * log(te(x)/(se+ieps)) -
                         se*(se+2.*te(x))/(4.*pow(se+te(x),2.)) * (pow(log(te(x)/(se+ieps)),2.)+pi*pi));
    }

std::complex<double> I5ggSU(double x)
    {
    return alpha/twopi*(se/(2.*(se+ue(x)))*log(ue(x)/(se+ieps)) -
                        se*(se+2.*ue(x))/(4.*pow(se+ue(x),2.))*(pow(log(ue(x)/(se+ieps)),2.)+pi*pi));
    }

std::complex<double> I5ggTU(double x)
    {
    return alpha/twopi*(te(x)/(2.*(te(x)+ue(x)))*log(ue(x)/(te(x)+ieps)) -
                        te(x)*(te(x)+2.*ue(x))/(4.*pow(te(x)+ue(x),2.))*(pow(log(ue(x)/(te(x)+ieps)),2.)+pi*pi));
    }

std::complex<double> I5ggTS(double x)
    {
    return alpha/twopi*(te(x)/(2.*(te(x)+se))*log(se/(te(x)+ieps)) -
                        te(x)*(te(x)+2.*se)/(4.*pow(te(x)+se,2.))*(pow(log(se/(te(x)+ieps)),2.)+pi*pi));
    }

std::complex<double> I5ggUT(double x)
    {
    return alpha/twopi*(ue(x)/(2.*(ue(x)+te(x)))*log(te(x)/(ue(x)+ieps)) -
                        ue(x)*(ue(x)+2.*te(x))/(4.*pow(ue(x)+te(x),2.))*(pow(log(te(x)/(ue(x)+ieps)),2.)+pi*pi));
    }

std::complex<double> I5ggUS(double x)
    {
    return alpha/twopi*(ue(x)/(2.*(ue(x)+se))*log(se/(ue(x)+ieps)) -
                        ue(x)*(ue(x)+2.*se)/(4.*pow(ue(x)+se,2.))*(pow(log(se/(ue(x)+ieps)),2.)+pi*pi));
    }



std::complex<double> I5gzTS(double x)
    {
    return alpha/twopi*(te(x)-mz2)/(2.*(te(x)+se))*(
               std::log(se_c/(te(x)-mz2)) - mz2/te(x) * std::log(mz2/(mz2-te(x)))
               +(te(x)+2.*se+mz2)/(te(x)+se) * (TMath::DiLog(te(x)/mz2)-
                       TMath::DiLog(-se/mz2)+std::log(-se_c/mz2)*std::log((mz2-te(x))/(mz2+se))));
    }


std::complex<double> I5gzTU(double x)
    {
    return alpha/twopi*(te(x)-mz2)/(2.*(te(x)+ue(x)))*(
               std::log(ue(x)/(te(x)-mz2))-mz2/te(x)*std::log(mz2/(mz2-te(x)))
               +(te(x)+2.*ue(x)+mz2)/(te(x)+ue(x))*(TMath::DiLog(te(x)/mz2)-
                       TMath::DiLog(-ue(x)/mz2)+std::log(-ue(x)/mz2)*std::log((mz2-te(x))/(mz2+ue(x)))));
    }


std::complex<double> I5gzUS(double x)
    {
    return alpha/twopi*(ue(x)-mz2)/(2.*(ue(x)+se))*(
               std::log(se_c/(ue(x)-mz2))-mz2/ue(x)*std::log(mz2/(mz2-ue(x)))
               +(ue(x)+2.*se+mz2)/(ue(x)+se)*(TMath::DiLog(ue(x)/mz2)-
                       TMath::DiLog(-se/mz2)+std::log(-se_c/mz2)*std::log((mz2-ue(x))/(mz2+se))));
    }


std::complex<double> I5gzUT(double x)
    {
    return alpha/twopi*(ue(x)-mz2)/(2.*(ue(x)+te(x)))*(
               std::log(te(x)/(ue(x)-mz2))-mz2/ue(x)*std::log(mz2/(mz2-ue(x)))
               +(ue(x)+2.*te(x)+mz2)/(ue(x)+te(x))*(TMath::DiLog(ue(x)/mz2)-
                       TMath::DiLog(-te(x)/mz2)+std::log(-te(x)/mz2)*std::log((mz2-ue(x))/(mz2+te(x)))));
    }
std::complex<double> I5gzST(double x)
    {
    return alpha/twopi*(se-mz2)/(2.*(se+te(x)))*(
               std::log(te(x)/(se_c-mz2))-mz2/se*std::log(mz2/(mz2-se))
               +(se+2.*te(x)+mz2)/(se+te(x))*(TMath::DiLog(se/mz2)-
                       TMath::DiLog(-te(x)/mz2)+std::log(-te(x)/mz2)*std::log((mz2-se)/(mz2+te(x)))));
    }



std::complex<double> I5gzSU(double x)
    {
    return alpha/twopi*(se-mz2)/(2.*(se+ue(x)))*(
               std::log(ue(x)/(se_c-mz2))-mz2/se*std::log(mz2/(mz2-se_c))
               +(se+2.*ue(x)+mz2)/(se+ue(x))*(TMath::DiLog(se/mz2)-
                       TMath::DiLog(-ue(x)/mz2)+std::log(-ue(x)/mz2)*std::log((mz2-se)/(mz2+ue(x)))));
    }

///////////////////////////////////////////////////////////////////
//More Definitions for D&P paper
///////////////////////////////////////////////////////////////////

double Z = 3.*log(se/(me*me))+2.*pi*pi/3.-4.;

double Xp (double x)
    {
    return pow(log(ue(x)/te(x)),2.)+pi*pi/3.;
    }

double Yt (double x)
    {
    return -pow(log(-te(x)/se),2.)-2.*log(-te(x)/se)*log((te(x)+se)/se)+3.*log(-te(x)/se)-pi*pi;
    }

double Yu (double x)
    {
    return -pow(log(-ue(x)/se),2.)-2.*log(-ue(x)/se)*log((ue(x)+se)/se)+3.*log(-ue(x)/se)-pi*pi;
    }


double gam (double x, double dE)
    {
    return 4.*alpha/pi*log(2.*dE/sqrt(se))*(log(ue(x)*te(x)/(se*me*me))-1.);
    }

double dut (double x)
    {
    return -log(-te(x)/se)*(log((mz2-ue(x))/(-ue(x)))+
                            log((mz2-ue(x))/mz2))+TMath::DiLog((mz2+te(x))/mz2)-TMath::DiLog((mz2+se)/se);
    }
double dtu (double x)
    {
    return -log(-ue(x)/se)*(log((mz2-te(x))/(-te(x)))+
                            log((mz2-te(x))/mz2))+TMath::DiLog((mz2+ue(x))/mz2)-TMath::DiLog((mz2+se)/se);
    }


double gzm = (2.*sw2-1.)/(2.*sqrt(sw2*cw2));
double gzp = sqrt(sw2/cw2);


///////////////////////////////////////////////////////////////////
//Matrix Elements, D&P
///////////////////////////////////////////////////////////////////

double M1u (double x)
    {
    return 2.*se/ue(x);
    }

double M2u (double x)
    {
    return 2.*te(x)/ue(x);
    }

double M1t (double x)
    {
    return 2.*se/te(x);
    }

double M3t (double x)
    {
    return 2.*ue(x)/te(x);
    }


double M1uzp (double x)
    {
    return gzp*gzp*2.*se/(ue(x)-mz2);
    }

double M1uzm (double x)
    {
    return gzm*gzm*2.*se/(ue(x)-mz2);
    }

double M2uz (double x)
    {
    return gzp*gzm*2.*te(x)/(ue(x)-mz2);
    }

double M1tzp (double x)
    {
    return gzp*gzp*2.*se/(te(x)-mz2);
    }

double M1tzm (double x)
    {
    return gzm*gzm*2.*se/(te(x)-mz2);
    }

double M3tz (double x)
    {
    return gzp*gzm*2.*ue(x)/(te(x)-mz2);
    }



///////////////////////////////////////////////////////////////////
//Correction Factors, D&P
///////////////////////////////////////////////////////////////////

std::complex<double> delu1g (double x)
    {
    return alpha/twopi*(Z+Yu(x)+Xp(x))+2.*I5ggUT(x);
    }
std::complex<double> delu2g (double x)
    {
    return alpha/twopi*(Z+Yu(x)+Xp(x))-2.*I5ggUS(x);
    }

std::complex<double> delt1g (double x)
    {
    return alpha/twopi*(Z+Yt(x)+Xp(x))+2.*I5ggTU(x);
    }

std::complex<double> delt3g (double x)
    {
    return alpha/twopi*(Z+Yt(x)+Xp(x))-2.*I5ggTS(x);
    }

std::complex<double> delu1z (double x)
    {
    return alpha/twopi*(Z+Yu(x)+Xp(x)+2.*dut(x))+4.*I5gzUT(x);
    }

std::complex<double> delu2z (double x)
    {
    return alpha/twopi*(Z+Yu(x)+Xp(x)+2.*dut(x))-4.*I5gzUS(x);

    }

std::complex<double> delt1z (double x)
    {
    return alpha/twopi*(Z+Yt(x)+Xp(x)+2.*dtu(x))+4.*I5gzTU(x);

    }

std::complex<double> delt3z (double x)
    {
    return alpha/twopi*(Z+Yt(x)+Xp(x)+2.*dtu(x))-4.*I5gzTS(x);

    }

///////////////////////////////////////////////////////////////////
//Soft Correction = Delta (CS), D&P
///////////////////////////////////////////////////////////////////

double corr_soft(double x, double dE)
    {
    return 0.25*alpha*alpha/(4.*se)*(

               //M1, lambda +, Z
               std::real(M1u(x)*M1u(x) * (delu1g(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1uzp(x) * (delu1g(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1tzp(x) * (delu1g(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1t(x) * (delu1g(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               std::real(M1t(x)*M1uzp(x) * (delt1g(x)+std::conj(delu1z(x))+gam(x,dE)))+
               std::real(M1t(x)*M1u(x) * (delt1g(x)+std::conj(delu1g(x))+gam(x,dE)))+
               std::real(M1t(x)*M1t(x) * (delt1g(x)+std::conj(delt1g(x))+gam(x,dE)))+
               std::real(M1t(x)*M1tzp(x) * (delt1g(x)+std::conj(delt1z(x))+gam(x,dE)))+

               std::real(M1uzp(x)*M1u(x) * (delu1z(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1uzp(x)*M1uzp(x) * (delu1z(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1uzp(x)*M1tzp(x) * (delu1z(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1uzp(x)*M1t(x) * (delu1z(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               std::real(M1tzp(x)*M1u(x) * (delt1z(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1tzp(x)*M1uzp(x) * (delt1z(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1tzp(x)*M1tzp(x) * (delt1z(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1tzp(x)*M1t(x) * (delt1z(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               //M1, lambda -, Z
               std::real(M1u(x)*M1u(x) * (delu1g(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1uzm(x) * (delu1g(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1tzm(x) * (delu1g(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1u(x)*M1t(x) * (delu1g(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               std::real(M1t(x)*M1uzm(x) * (delt1g(x)+std::conj(delu1z(x))+gam(x,dE)))+
               std::real(M1t(x)*M1u(x) * (delt1g(x)+std::conj(delu1g(x))+gam(x,dE)))+
               std::real(M1t(x)*M1t(x) * (delt1g(x)+std::conj(delt1g(x))+gam(x,dE)))+
               std::real(M1t(x)*M1tzm(x) * (delt1g(x)+std::conj(delt1z(x))+gam(x,dE)))+

               std::real(M1uzm(x)*M1u(x) * (delu1z(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1uzm(x)*M1uzm(x) * (delu1z(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1uzm(x)*M1tzm(x) * (delu1z(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1uzm(x)*M1t(x) * (delu1z(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               std::real(M1tzm(x)*M1u(x) * (delt1z(x) + std::conj(delu1g(x)) + gam(x,dE)) )+
               std::real(M1tzm(x)*M1uzm(x) * (delt1z(x) + std::conj(delu1z(x)) + gam(x,dE)) )+
               std::real(M1tzm(x)*M1tzm(x) * (delt1z(x) + std::conj(delt1z(x)) + gam(x,dE)) )+
               std::real(M1tzm(x)*M1t(x) * (delt1z(x) + std::conj(delt1g(x)) + gam(x,dE)) )+

               //M2, Z
               std::real(M2u(x)*M2u(x) * (delu2g(x) + std::conj(delu2g(x)) + gam(x,dE)) )+
               std::real(M2u(x)*M2uz(x) * (delu2g(x) + std::conj(delu2z(x)) + gam(x,dE)) )+
               std::real(M2uz(x)*M2u(x) * (delu2z(x) + std::conj(delu2g(x)) + gam(x,dE)) )+
               std::real(M2uz(x)*M2uz(x) * (delu2z(x) + std::conj(delu2z(x)) + gam(x,dE)) )+

               //M3, Z

               std::real(M3t(x)*M3t(x) * (delt3g(x) + std::conj(delt3g(x)) + gam(x,dE)) )+
               std::real(M3t(x)*M3tz(x) * (delt3g(x) + std::conj(delt3z(x)) + gam(x,dE)) )+
               std::real(M3tz(x)*M3t(x) * (delt3z(x) + std::conj(delt3g(x)) + gam(x,dE)) )+
               std::real(M3tz(x)*M3tz(x) * (delt3z(x) + std::conj(delt3z(x)) + gam(x,dE)) )+

               //M2, Z, again - repeated since the unpolarized CS is (1/4)*Sum_lambda(polarized CS)
               std::real(M2u(x)*M2u(x) * (delu2g(x) + std::conj(delu2g(x)) + gam(x,dE)) )+
               std::real(M2u(x)*M2uz(x) * (delu2g(x) + std::conj(delu2z(x)) + gam(x,dE)) )+
               std::real(M2uz(x)*M2u(x) * (delu2z(x) + std::conj(delu2g(x)) + gam(x,dE)) )+
               std::real(M2uz(x)*M2uz(x) * (delu2z(x) + std::conj(delu2z(x)) + gam(x,dE)) )+

               //M3, Z, again

               std::real(M3t(x)*M3t(x) * (delt3g(x) + std::conj(delt3g(x)) + gam(x,dE)) )+
               std::real(M3t(x)*M3tz(x) * (delt3g(x) + std::conj(delt3z(x)) + gam(x,dE)) )+
               std::real(M3tz(x)*M3t(x) * (delt3z(x) + std::conj(delt3g(x)) + gam(x,dE)) )+
               std::real(M3tz(x)*M3tz(x) * (delt3z(x) + std::conj(delt3z(x)) + gam(x,dE)) )


           );}
