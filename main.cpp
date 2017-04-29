#include <iostream>
#include<cmath>
#include<cstdlib>

using namespace std;
const int n=4;

double maxi(double a, double b){
return (a>b)? a:b;
}
//implementation of normal cdf from :http://www.johndcook.com/cpp_phi.html
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}
class BestOfPut{
private:
	double r; double T; double K; double S[n]; double vol[n]; double corr[n][n];  double L[n][n];
	void Choleski();
	double simulation();
public:
	BestOfPut(){
	int i,j;
	double c;

	//entering option parameters
	cout << "Enter the risk free rate " << endl;
	cin>>r;
	cout << "Enter time to maturity of the option " << endl;
	cin>>T;
	cout << "Enter the strike of the option " << endl;
	cin>>K;

	//entering stocks parameters
	for(i=0; i<n; i++){
	cout << "Enter the price of stock " << i+1 << endl;
	cin>>S[i];
	cout << "Enter the volatility of stock " << i+1 << endl;
	cin>>vol[i];
	if(i>0){
		for (j=0; j<i; j++){
			c=2;
			do{
			cout << "Enter the correlation of stock " << i+1 << " with stock " << j+1 <<endl;
			cin>>c;
			//restriction on correlation matrix ensures it is positive-semidefinite
			if (c<-1||c>1) cout << "Correlation must be a real number between -1 and 1"<<endl;
			}while(c<-1||c>1);
			//symmetry of correlation matrix
			corr[i][j]=c;
			corr[j][i]=c;
		}
	}	
	corr[i][i]=1;
	}
	Choleski();

}
	double montecarlo(int N);

};

void  BestOfPut::Choleski(){
	
	//Choleski decomposition of corr matrix into LL'
	int i,j,k;
	double sum;

	//set L matrix equal to 0
	for(i=0; i<n; i++) {
		for (j=0; j<n; j++)
			L[i][j]=0;
	}

	// compute 1st column
	for(i=0; i<n; i++) 
		L[i][0]=corr[i][0];

	// compute jth column
	for (j=1; j<n; j++){
		
		//compute diagonal term
		sum=0;
		for(k=0; k<j; k++)
			sum +=L[j][k]*L[j][k];
		L[j][j]=sqrt(corr[j][j]-sum);

		//compute other terms
		for(i=j+1; i<n; i++){
			sum=0;
			for(k=0; k<j; k++)
				sum +=L[i][k]*L[j][k];
			L[i][j]=(corr[i][j]-sum)/L[j][j];
		}
	
	}

}

void Marsagliapolar(double& z1, double& z2){

      double u1, u2, v1, v2,s, sq;
      
	  do{
      u1=(double) rand()/RAND_MAX;
      u2=(double) rand()/RAND_MAX;
      v1=2*u1-1;
      v2=2*u2-1;
      s=v1*v1+v2*v2;
      }while(s>1 || s==0);

      sq=sqrt(-2*log(s)/s);
      z1=v1*sq;
      z2=v2*sq;
}
double BestOfPut::simulation(){
	int i,j;
	double Z[n], X[n];
	double sum, sumput=0, max1=0, max2=0, ST1, ST2;

	//simulating independent standard normal rvs
	for (i=0; i<n; i+=2)
		Marsagliapolar(Z[i],Z[i+1]);
	
	//introducing correlation
	for (i=0; i<n; i++){
		sum=0;
		for(j=0; j<n; j++)
			sum +=Z[j]*L[i][j];
		X[i]=sum;
	}

	//calculating stock prices at maturity
	for (i=0; i<n; i++){
		ST1=S[i]*exp((r-0.5*vol[i]*vol[i])*T+vol[i]*sqrt(T)*X[i]);
		if(ST1>max1) max1=ST1;
		sumput +=maxi(0,K-ST1); // control variate - average price of n put options

		// variance reduction - antithetic rv
		ST2=S[i]*exp((r-0.5*vol[i]*vol[i])*T-vol[i]*sqrt(T)*X[i]);
		if(ST2>max2) max2=ST2;
		sumput +=maxi(0,K-ST2);
	}

return (maxi(0,K-max1)+maxi(0,K-max2)-sumput/n)/2;
}

double BestOfPut::montecarlo(int N){
	int i;
	double sum=0, d1, d2, put, sumput=0;

	for (i=1; i<N; i++)
		sum +=simulation();

	//control variate - average put analytical price
	for (i=0; i<n; i++){
		d1=(log(S[i]/K)+(r+0.5*vol[i]*vol[i])*T)/(vol[i]*sqrt(T));
		d2=d1-vol[i]*sqrt(T);
		put=K*exp(-r*T)*phi(-d2)-S[i]*phi(-d1);
		sumput +=put;
	}

	return sum/N*exp(-r*T)+sumput/n;
}

int main(){
	srand(5);
	BestOfPut bop;
	cout << "The price of Best-of-put option is: "<< bop.montecarlo(1000) <<endl;;
	
	return 0;
}
