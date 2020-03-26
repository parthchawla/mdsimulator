#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

int main(){

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Input Parameters ////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	double dt = 0.005;
	int m = 30;
	int N = pow(m,2); // no. of particles
	double T0 = 1.0; // target temperature
	double rho = 0.5; // number density
	double A = N/rho; // area
	double L = m/pow(rho, 0.5);
	double rc = 2.5; // cut-off
	int t = 103000; // no. of timesteps
	int equi = 3000;
	int runtime = t - equi; // timesteps without equilibration

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Declare Variables and Arrays ////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	double rx[N];
	double ry[N];
	double vx[N];
	double vy[N];
	double fx[N];
	double fy[N];
	double vvx; // average x-velocity
	double vvy; // average y-velocity
	double x;
	double y;
	double r;
	double ff; // force
	double ffx;
	double ffy;
	double u; // potential energy
	double k; // kinetic energy
	double T; // temperature
	double Pc; // interaction part of pressure
	double Pk; // kinetic part of pressure
	double P; // total pressure
	int i;
	int j;
	int n;
	double sum;
	double sumx;
	double sumy;
	double diff;
	std::vector<double> energy;
	std::vector<double> pressure;
	std::vector<double> bavg1;
	std::vector<double> bavg2;
	double avg1;
	double var1;
	double avg2;
	double var2;
	double tb;
	double blockavg;
	double blockvar;
	double s;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Set up Lattice and Initialise Velocities ////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << "\nSetting up lattice." << std::endl;
	n = 0;
	for (i = 1; i <= m; i++){ // set up lattice

		for (j = 1; j <= m; j++){

			ry[n] = (L/m)*(j - 0.5);
			rx[n] = (L/m)*(i - 0.5);
			n++;

		} 
	}

	for (i = 0; i <= N-1; i++){ // give random values as initial velocities

		vx[i] = sqrt(i);
		vy[i] = cos(5*i);

	}

	std::cout << "\nInitialising Velocities." << std::endl;
	sumx = 0.0;
	sumy = 0.0;

	for (i = 0; i <= N-1; i++){ // calculate sum of all velocities

		sumx = sumx + vx[i];
		sumy = sumy + vy[i];

	}

	vvx = sumx/N; // calculate average velocity
	vvy = sumy/N;

	T = 0.0;
	for (i = 0; i <= N-1; i++){ // find intial temperature

		T = T + 0.5*(pow(vx[i], 2) + pow(vy[i], 2));

	}

	for (i = 0; i <= N-1; i++){ // initialise velocity by scaling wrt temperature

		vx[i] = vx[i] - vvx; // sets box momentum to zero
		vy[i] = vy[i] - vvy;
		vx[i] = vx[i]*sqrt(T0/T); // scale wrt target temperature and calculated temperature
		vy[i] = vy[i]*sqrt(T0/T);

	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Force Loop //////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	std::ofstream file1;
	std::ofstream file2;
	std::ofstream file3;
	std::ofstream file4;
	std::ofstream file5;
	std::ofstream file6;
	file1.open("U.txt");
	file2.open("K.txt");
	file3.open("E.txt");
	file4.open("P.txt");
	file5.open("sE.txt");
	file6.open("sP.txt");

	std::cout << "\nCalculating Velocities, Positions, Energies and Pressure."
	<< std::endl;
	std::cout << "Please wait..." << std::endl;
	sumx = 0.0;
	sumy = 0.0;
	for (int step = 0; step < t; step++){

		for (n = 0; n <= N-1; n++){ // initialise force

			fx[n] = 0;
			fy[n] = 0;

		}

		////////////////////////////////////////////////////////////////////////////////////////////
		// 	Calculate Force and Potential Energy ///////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////

		// initialise potential energy and pressure
		u = 0;
		Pc = 0;
		Pk = 0;
		P = 0;

		for (i = 0; i <= N-2; i++){

			for (j = i + 1; j <= N-1; j++){

				x = rx[i] - rx[j];
				y = ry[i] - ry[j];

				// periodic boundary conditions
				if (x >= L/2){ x = x - L; }
				if (x < -L/2){ x = x + L; }
				if (y >= L/2){ y = y - L; }
				if (y < -L/2){ y = y + L; }

				r = sqrt(pow(x, 2) + pow(y, 2));

				if (r < rc){ // check if r less than cut-off

					u = u + 4*(1/pow(r,12) - 1/pow(r,6)); // calculate potential energy
					ff = 48/pow(r, 14) - 24/pow(r, 8); // calculate force
					ffx = x*ff;
					ffy = y*ff;
					fx[i] = fx[i] + ffx;
					fy[i] = fy[i] + ffy;
					fx[j] = fx[j] - ffx;
					fy[j] = fy[j] - ffy;
					Pc = Pc + r*ff;

				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		// Apply Verlet Algorithm, calculate Kinetic Energy and Temperature/////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////
		
		k = 0; // initialise kinetic energy

		for (i = 0; i <= N-1; i++){

			vx[i] = vx[i] + dt*fx[i]; // update velocity using acceleration
			vy[i] = vy[i] + dt*fy[i];
			k = k + 0.5*(pow(vx[i], 2) + pow(vy[i], 2)); // calculate kinetic energy
			rx[i] = rx[i] + dt*vx[i]; // update position using velocity
			ry[i] = ry[i] + dt*vy[i];

		}

		u = u/N; // average potential energy
		k = k/N; // average kinetic energy
		T = k; // temperature
		Pk = rho*T;
		Pc = 0.5*Pc*rho/N;
		P = Pk + Pc;

		////////////////////////////////////////////////////////////////////////////////////////////
		// Write to Text Files /////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////

		if (step >= equi){

			sumx += u+k;
			sumy += P;
			energy.push_back(u+k);
			pressure.push_back(P);
			file1 << u << "\n";
			file2 << k << "\n";
			file3 << u+k << "\n";
			file4 << P << "\n";

		}

	}

	file1.close();
	file2.close();
	file3.close();
	file4.close();

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Analysis ////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << "\nCalculating Mean and Variance for Energy and Pressure." << std::endl;
	avg1 = sumx/runtime;
	avg2 = sumy/runtime;

	diff = 0.0;
	for (i = 0; i < energy.size(); i++) {

		diff += pow(energy[i] - avg1, 2);

	}
	var1 = diff/runtime;

	diff = 0.0;
	for (i = 0; i < pressure.size(); i++) {

		diff += pow(pressure[i] - avg2, 2);

	}
	var2 = diff/runtime;

	std::cout << "\nCalculating Block Averages and Block Variances for Energy and Pressure."
	<< std::endl;
	for (;;){

		std::cout << "Enter the block size or enter 0 to end the simulation."
		<< std::endl;
		std::cin >> tb;
		if (tb == 0){
			break;
		}

		n = 0;
		for (i = 0; i < energy.size(); i += tb) {

			sum = 0.0;
			for (j = 0; j < tb; j++) {

				sum += energy[n];
				n++;

			}
			blockavg = sum / tb;
			bavg1.push_back(blockavg);

		}
		diff = 0.0;
		for (i = 0; i < bavg1.size(); i++){

			diff += pow(bavg1[i] - avg1, 2);

		}
		blockvar = diff*tb/runtime;
		s = blockvar*tb/var1;
		file5 << s << "\n";

		n = 0;
		for (i = 0; i < pressure.size(); i += tb) {

			sum = 0.0;
			for (j = 0; j < tb; j++) {

				sum += pressure[n];
				n++;

			}
			blockavg = sum / tb;
			bavg2.push_back(blockavg);

		}
		diff = 0.0;
		for (i = 0; i < bavg2.size(); i++){

			diff += pow(bavg2[i] - avg2, 2);

		}
		blockvar = diff*tb/runtime;
		s = blockvar*tb/var2;
		file6 << s << "\n";
		bavg1.clear();
		bavg2.clear();

	}

	file5.close();
	file6.close();

}