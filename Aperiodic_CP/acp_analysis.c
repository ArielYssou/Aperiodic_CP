#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <dirent.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define UNUSED(expr) do { (void)(expr); } while (0)

#define Points 250
#define CritLambda 3.29785

void LogSpaced(double arr[], double begin, double end) {
	int index;
	double x, dx;
	if( begin == 0 )
		begin = 1.0;
	dx = (log(end) - log(begin)) / Points;
	x = log(begin);
	for(index = 0; index < Points; index++) {
		arr[index]  = exp(x);
		x += dx;
	}
	arr[Points] = end;
}

int main(int argc, char **argv) {
	char path[1000], str[1000];
	int i, k, sim, Size, sim_init, sim_tot, sim_final, inter, dead, Elems[Points+1], Died[Points+1];
	double la, lb, tmax, t, rho, rho_desv, time_desv, delta_t, survivors;
	double Tsups[Points+1], Rho_Sum[Points+1], Time_Sum[Points+1], Rho_S2[Points+1], Time_S2[Points+1];
	FILE *input_file, *rho_file, *ps_file;

	k = strtof(argv[1], NULL);
	delta_t = strtof(argv[2], NULL);
	la = strtof(argv[3], NULL);
	lb = strtof(argv[4], NULL);
	tmax = strtod(argv[5], NULL);
	Size = (int)strtod(argv[6], NULL);
	sim_init = strtod(argv[7], NULL);
	sim_final = strtod(argv[8], NULL);

	LogSpaced(Tsups, 1, tmax);

	sprintf(path, "/home/ariel/Simulacoes/data_aperiodic");
	sprintf(str, "/k=%d", k);
	strcat(path, str);
	sprintf(str, "/Delta_t=%.2f", delta_t);
	strcat(path, str);
	sprintf(str, "/Lambda_a=%.3f", la);
	strcat(path, str);
	sprintf(str, "/Lambda_b=%.8f", lb);
	strcat(path, str);
	sprintf(str, "/Tmax=%e", tmax);
	strcat(path, str);
	sprintf(str, "/Size=%d", (int)Size);
	strcat(path, str);

	for(i = 0; i < Points; i++) {
		Rho_Sum[i] = Time_Sum[i] = 0;
		Rho_S2[i] = Time_S2[i] = 0;
		Died[i] = Elems[i] = 0;
	}

	sim_tot = sim_final;
	for(sim = sim_init; sim < sim_final; sim++) {
		sprintf(str, "%s/rho_%d.dat", path, sim);

		if( access( str, F_OK ) != -1 ) {
			input_file=fopen(str, "r");
			inter = 0;
			while((fscanf(input_file, "%lf\t%lf\n", &t, &rho)) != EOF) {
				while( t > Tsups[inter]) {
					inter++;
				}
				if( rho <= 0 ) {
					Died[inter] += 1.0;
					break;
				}
				Rho_Sum[inter] += rho;
				Rho_S2[inter] += rho*rho;
				Time_Sum[inter] += t;
				Time_S2[inter] += t*t;
				Elems[inter] += 1;
			}
			fclose(input_file);
		} else {
				sim_tot--;
		}
	}

	sprintf(str, "%s/rho_av.dat", path);
	rho_file=fopen(str, "w+");

	sprintf(str, "%s/surv_prob.dat", path);
	printf("%s", str);
	ps_file=fopen(str, "w+");

	dead = 0; // Corrects the av. density acconting each dead sim
	survivors = sim_tot;
	for(i = 0; i < Points; i++) {
		if(survivors == 0)
			break;
		else {
			survivors -= Died[i];
			dead += Died[i];
		}

		if (Elems[i] != 0) {
			Rho_Sum[i] /= (double)(Elems[i] + dead);
			Time_Sum[i] /= (double)(Elems[i]);

			rho_desv = sqrt((double)(Rho_S2[i] / (Elems[i] + dead) ) - pow(Rho_Sum[i], 2.));
			time_desv = sqrt((double)(Time_S2[i] / (Elems[i]) ) - pow(Time_Sum[i], 2.));

			fprintf(rho_file, "%f,%f,%f,%f\n", Time_Sum[i], Rho_Sum[i], time_desv, rho_desv);
			fprintf(ps_file, "%f,%f\n", Time_Sum[i], (double)(survivors / (sim_tot - sim_init)));
		}
	}

	fclose(rho_file);
	fclose(ps_file);

	return 0;
}
 
