#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>
#include <utility>
#define PI 3.14159265

using namespace std;


//
// note: should be run with the -std=c++11 option.
//



const int n = 20;
double U[n][n];
double S = n;
double H = -sqrt(n)*S;
double thetaRange = 0.24367; // this is just the initial guess. It will be dynamically adjusted in the first block.
bool tune_thetaRange;
double beta;
double betaRn;
double theta;
bool rowColumn;
int index1;
int index2;
ofstream theta_log;


// set up random number generator. Use the clock as a seed.
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
mt19937 generator (seed);  // generate random numbers using the Mersenne Twister
uniform_real_distribution<double> uniform(0.0,1.0);
uniform_int_distribution<> coinFlip(0,1);
uniform_int_distribution<> random_index(0, n-1);


void initializeU() {
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      if (i==j)
	U[i][j] = 1;
      else
	U[i][j] = 0;
    }
  }
}


bool accept(double dS) {
  if (dS >= 0)
    return 1;
  else {
    double p = exp(betaRn*dS);
    if (tune_thetaRange == 1) { // this option is only turned on in the first block and is used to optimize thetaRange
      //thetaRange *= exp(p - 0.5);
      if (p>0.5){
	thetaRange *= 1.0001;
	if (thetaRange > 2*PI)
	  thetaRange = 2*PI;
      }
      else
	thetaRange *= 1.0/1.0001;
      theta_log << thetaRange << endl;
    } // end tuning-specific part of code
    double number = uniform(generator);
    if (number <= p)
      return 1;
    else
      return 0;
  }
}

bool givensRotation() {
  // "inputs" are the global parameters index1, index2, theta, and rowColumn
  // rowColum = 0 if we want to rotate rows; 1 if rotate columns
  // indices are the indices of the rows or columns we are rotating
  // theta is the angle of rotation in radians
  // implements a givens rotation on array U and computes the corresponding change in H
  double c = cos(theta);
  double s = sin(theta);
  double Unew[n][2];
  double dS = 0;
  // First compute the rotation and find the change in energy. Then decide whether to accept the rotation and implement if appropriate.
  if (rowColumn == 0){
    for(int i=0; i<n; i++){
      Unew[i][0] = c*U[index1][i] + s*U[index2][i];
      Unew[i][1] = c*U[index2][i] - s*U[index1][i];
      dS += abs(Unew[i][0]) - abs(U[index1][i])+ abs(Unew[i][1]) - abs(U[index2][i]);
    }
  }
  else {
    for(int i=0; i<n; i++){
      Unew[i][0] = c*U[i][index1] + s*U[i][index2];
      Unew[i][1] = c*U[i][index2] - s*U[i][index1];
      dS += abs(Unew[i][0]) - abs(U[i][index1]) + abs(Unew[i][1]) - abs(U[i][index2]);
    }
  }
  // cout << "dS = " << dS << endl;
  if (accept(dS)) { // implement the rotation if it is accepted
    S += dS;
    if (rowColumn == 0) {
      for(int i=0; i<n; i++){
	U[index1][i] = Unew[i][0];
	U[index2][i] = Unew[i][1];
      }
    }
    else {
      for(int i=0; i<n; i++){
	U[i][index1] = Unew[i][0];
	U[i][index2] = Unew[i][1];
      }
    }
    return 0; // rotation accepted
  }
  else
    return 1; // rotation not accepted
}

void randomParams() {
  theta = thetaRange * (uniform(generator) - 0.5);
  rowColumn = coinFlip(generator);
  index1 = random_index(generator);
  index2 = random_index(generator);
  while (index1 == index2)
    index2 = random_index(generator);
}


void updateU() {
  // perform givens rotations on a pair of rows or columns
  randomParams();
  while (givensRotation())
    randomParams();
  H = -sqrt(n)*S;
}





pair<double, double> runningAverage(int blockSize) {
  //returns the average energy density (first output) and heat capacity (second output), both per unit volume
  double runningE = 0;
  double runningE2 = 0;
  for(int i = 0; i< blockSize; i++){
    updateU();
    runningE += H;
    runningE2 += H*H;
  }
  double avg_E = runningE/blockSize;
  double avg_e = avg_E/(n*n);
  double avg_E2 = runningE2/blockSize;
  double avg_e2 = avg_E2/(n*n*n*n);
  return make_pair(avg_e,avg_e2);
}


void displayU() {
  // displays the current matrix U.
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
	cout << U[i][j] << " ";
    }
    cout<<endl;
  }
}


int main() {
  initializeU();
  
  int block_size = n*(n-1)*1E4;
  int number_of_blocks = 1;
  beta = 6;
  betaRn = beta*sqrt(n);

  // don't tune thetaRange
  tune_thetaRange = 0; // stop tuning the thetaRange.

  // generate sample
  for(int i = 0; i<number_of_blocks; i++){
    pair<double, double> e_and_e2 = runningAverage(block_size);
    double e = e_and_e2.first;
    double e2 = e_and_e2.second;
  }

  return 0;
}
