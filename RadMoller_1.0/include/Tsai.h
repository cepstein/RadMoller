//Tsai's Soft Corrections

double corr_soft_tsai(double x, double dE){
	return (-4.*(1./137.)/pi*(23./18.-11./12.*log(-te(x)/pow(me,2.))+
		log(Ecmp/dE)*(log((te(x)*ue(x))/(se*me*me))-1.)));}
