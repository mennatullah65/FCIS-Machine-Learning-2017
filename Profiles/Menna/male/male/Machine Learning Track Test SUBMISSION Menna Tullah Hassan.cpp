#include<iostream>
#include<algorithm>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#define Rand() ((double)rand()/RAND_MAX)

using namespace std;
double **M, **U, **V,**Vinverse,**Mdash;
int r, c;
//matrix calculations
void printMat(double **M,int r, int c){
	for (int i = 0; i < r; i++){
		cout << ((i == 0) ? "[[" : " [");
		for (int j = 0; j < c; j++){
			cout << M[i][j];
			cout << ((j == c - 1) ? "]" : " ");
		}
		cout << ((i == r - 1) ? "]\n" : "\n");
	}
}
void initializeMat(double ** &M, int r, int c){
	M = new double*[r];
	for (int i = 0; i < r; i++){
		M[i] = new double[c];
	}
}
void enterMat(double **M, int r, int c){
	for (int i = 0; i < r; i++){
		for (int j = 0; j < c; j++){
			cin >> M[i][j];
		}
	}
}
void transposeMat(double **M,double **Minverse ,int r, int c){
	for (int j = 0; j < c; j++){
		for (int i = 0; i < r; i++){
			Minverse[i][j] = M[j][i];
		}
	}
}
void multiplyMat(int n, int k ,int m,double**Mdash,double**U,double**V){
	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<k; j++)
		{
			for (int l = 0; l<m; l++)
			{
				Mdash[i][j] += U[i][l] * V[l][j];
			}
		}
	}
}

int main(){
	r = c = 5;
	//Declare Original Matrix
	initializeMat(M, r, c);
	cout << "M = " << endl;
	enterMat(M, r, c);
	//inizialize Mats
	initializeMat(U, r, 2);
	initializeMat(V, r, 2);
	initializeMat(Vinverse, 2, r);
	initializeMat(Mdash, r, c);
	
	//indicate U & V
	/*Generate random variables of U and V till the error estimation error would be the square root of the Mean Squared Error*/
	int dim = 2; // problem dimensions
	int numParticles = 5;
	int maxEpochs = 1000;
	double exitError = 0.0; // exit early if reach this error
	double minX = -10.0; // problem-dependent
	double maxX = 10.0;
	vector<double> res(dim);
	res = Solve(dim, numParticles, minX, maxX, maxEpochs, exitError);
	

	//Print Mats
	transposeMat(V, Vinverse, r, c);
	multiplyMat(5, 2, 5,Mdash,U,V);
	
	cout << "U = " << endl;
	printMat(U, r, c);
	cout << "V = " << endl;
	printMat(U, r, c);
	cout << "M' = " << endl;
	printMat(U, r, c);
	
	return 0;
}


double Error(vector<double> x)
    {
      // 0.42888194248035300000 when x0 = -sqrt(2), x1 = 0
      double trueMin = -0.42888194; // true min for z = x * exp(-(x^2 + y^2))
      double z = x[0] * exp( -((x[0]*x[0]) + (x[1]*x[1])) );
      return (z - trueMin) * (z - trueMin); // squared diff
    }

class Particle
{
public:
	vector<double> position;
	double error;
	vector<double> velocity;
	vector<double> bestPosition;
	double bestError;
	Particle();
	Particle(vector<double>pos, double err, vector<double>vel, vector<double>bestPos, double bestErr)
	{
		this->position.reserve(pos.size());
		this->position = pos;

		this->error = err;

		this->velocity.reserve(vel.size());
		this->velocity = vel;

		this->bestPosition.reserve(bestPos.size());
		this->bestPosition = bestPos;
		
		this->bestError = bestErr;
	}
	void setterParticle(vector<double>pos, double err, vector<double>vel, vector<double>bestPos, double bestErr)
	{
		this->position.reserve(pos.size());
		//this->position = pos;
		copy(pos.begin(),pos.end(),this->position.begin());

		this->error = err;

		this->velocity.reserve(vel.size());
		//this->velocity = vel;
		copy(vel.begin(), vel.end(), this->velocity.begin());

		this->bestPosition.reserve(bestPos.size());
		//this->bestPosition = bestPos;
		copy(bestPos.begin(), bestPos.end(), this->bestPosition.begin());

		this->bestError = bestErr;
	}

}; // Particle

vector<double> Solve(int dim, int numParticles, double minX, double maxX, int maxEpochs, double exitError)
{
	// assumes existence of an accessible Error function and a Particle class
	//Random rnd = new Random(0);

	Particle *swarm = new Particle[numParticles];
	vector<double>bestGlobalPosition(dim); // best solution found by any particle in the swarm
	double bestGlobalError = 1000; // smaller values better

	// swarm initialization
	for (int i = 0; i < numParticles; ++i)
	{
		vector<double>randomPosition(dim);
		for (int j = 0; j < dim; ++j)
			randomPosition[j] = (maxX - minX) * Rand() + minX; // 

		double error = Error(randomPosition);
		vector<double>randomVelocity(dim);

		for (int j = 0; j < dim; ++j)
		{
			double lo = minX * 0.1;
			double hi = maxX * 0.1;
			randomVelocity[j] = (hi - lo) * Rand() + lo;
		}
		swarm[i].setterParticle(randomPosition, error, randomVelocity, randomPosition, error);

		// does current Particle have global best position/solution?
		if (swarm[i].error < bestGlobalError)
		{
			bestGlobalError = swarm[i].error;
			copy(bestGlobalPosition.begin(),bestGlobalPosition.end(), swarm[i].position.begin());
		}
	} // initialization

	// prepare
	double w = 0.729; // inertia weight.
	double c1 = 1.49445; // cognitive/local weight
	double c2 = 1.49445; // social/global weight
	double r1, r2; // cognitive and social randomizations
	double probDeath = 0.01;
	int epoch = 0;

	vector<double> newVelocity (dim);
	vector<double> newPosition (dim);
	double newError;

	// main loop
	while (epoch < maxEpochs)
	{
		for (int i = 0; i < numParticles; ++i) // each Particle
		{
			Particle currP = swarm[i]; // for clarity

			// new velocity
			for (int j = 0; j < currP.velocity.size(); ++j) // each component of the velocity
			{
				r1 = Rand();
				r2 = Rand();

				newVelocity[j] = (w * currP.velocity[j]) +
					(c1 * r1 * (currP.bestPosition[j] - currP.position[j])) +
					(c2 * r2 * (bestGlobalPosition[j] - currP.position[j]));
			}
			copy(newVelocity.begin(), newVelocity.end(), currP.velocity.begin());

			// new position
			for (int j = 0; j < currP.position.size(); ++j)
			{
				newPosition[j] = currP.position[j] + newVelocity[j];
				if (newPosition[j] < minX)
					newPosition[j] = minX;
				else if (newPosition[j] > maxX)
					newPosition[j] = maxX;
			}
			copy(newPosition.begin(), newPosition.end(), currP.position.begin());
			newError = Error(newPosition);
			currP.error = newError;

			if (newError < currP.bestError)
			{
				copy(newPosition.begin(), newPosition.end(), currP.bestPosition.begin());
				currP.bestError = newError;
			}

			if (newError < bestGlobalError)
			{
				copy(newPosition.begin(), newPosition.end(), bestGlobalPosition.begin());
				bestGlobalError = newError;
			}

			// death?
			double die = Rand();
			if (die < probDeath)
			{
				// new position, leave velocity, update error
				for (int j = 0; j < currP.position.size(); ++j)
					currP.position[j] = (maxX - minX) * Rand() + minX;
				currP.error = Error(currP.position);
				copy(currP.position.begin(), currP.position.end(), currP.bestPosition.begin());
				currP.bestError = currP.error;
				if (currP.error < bestGlobalError) // global best by chance?
				{
					bestGlobalError = currP.error;
					copy(currP.position.begin(), currP.position.end(), bestGlobalPosition.begin());

				}
			}

		}
		++epoch;
	}

	vector<double> result (dim);
	copy(bestGlobalPosition.begin(), bestGlobalPosition.end(), result.begin());
	return result;
} 





/*--------------------------------------------------------------------*/
/*takahama code for PSO*/

//PSO variables
/*
#define Nparticles	50
#define T_MAX		1000
#define W_0		0.9
#define W_T		0.4
#define MAX_V		2.0
#define c1		2.0
#define c2		2.0
#define Nvariables	5
#define better(y1, y2)	(y1<y2)
#define New(type, n, msg)	(type *)NewCell(sizeof(type), n, msg)

typedef struct {
double *x;
double *v;
double f;
double pbest;
double *x_star;
} ParticleRec, *Particle;
*/
//Code
/*
void Evaluate(Particle P)
{
int i;
P->f = 0.0;
for (i = 0; i<Nvariables; i++)
P->f += (P->x[i] - 1)*(P->x[i] - 1);
}
void UpdateBest(Particle P)
{
int j;
for (j = 0; j<Nvariables; j++) P->x_star[j] = P->x[j];
P->pbest = P->f;
}
int Initialize(Particle P, int n)
{
int i, j;
int G;		// the index of the best particle
G = 0;
for (i = 0; i<n; i++) {
for (j = 0; j<Nvariables; j++) {
P[i].x[j] = Rand();	// problem dependent
P[i].v[j] = 0.0;		// problem dependent
}
Evaluate(&P[i]);
UpdateBest(&P[i]);
if (better(P[i].f, P[G].f)) G = i;
}
return G;
}
void *NewCell(int size, int n, char *msg)
{
void *New;

if ((New = malloc(size*n)) == NULL) {
fprintf(stderr, "Cannot allocate memory for %d %s\n", n, msg);
exit(1);
}
return New;
}
Particle NewParticles(int n)
{
int i;
Particle P;
P = New(ParticleRec, n, "particles");
for (i = 0; i<n; i++) {
P[i].x = New(double, Nvariables, "x");
P[i].v = New(double, Nvariables, "v");
P[i].x_star = New(double, Nvariables, "x*");
}
return P;
}
void Print(Particle P)
{
int j;

for (j = 0; j<Nvariables; j++)
printf("%f ", P->x_star[j]);
printf(" = %g\n", P->pbest);
}
void PSO(){
int t, i, j;
Particle P;
int G;
double w;

P = NewParticles(Nparticles);
G = Initialize(P, Nparticles);
w = W_0;
for (t = 1; t <= T_MAX; t++) {
for (i = 0; i<Nparticles; i++) {
for (j = 0; j<Nvariables; j++) {
P[i].v[j] = w*P[i].v[j]
+ c1*Rand()*(P[i].x_star[j] - P[i].x[j])
+ c2*Rand()*(P[G].x_star[j] - P[i].x[j]);
if (P[i].v[j]<-MAX_V)
P[i].v[j] = -MAX_V;
else if (P[i].v[j]>MAX_V)
P[i].v[j] = MAX_V;
P[i].x[j] += P[i].v[j];
}
Evaluate(&P[i]);
if (better(P[i].f, P[i].pbest)) {
if (better(P[i].f, P[G].pbest)) G = i;
UpdateBest(&P[i]);
}
}
//	printf("%4d: ", t); Print(&P[G]);
w -= (W_0 - W_T) / T_MAX;
}
}
*/
/*takahama code for PSO*/
/*--------------------------------------------------------------------*/

