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

#define Points 1000
#define CritLambda 3.29785

int main(int argc, char **argv) {
	char path[1000], str[1000];
	int i, k, sim, Size, sim_init, sim_final, inter;
	double la, lb, tmax, t, rho, rho_desv, time_desv;
	double Rho_Sum[Points+1], Time_Sum[Points+1], Rho_S2[Points+1], Time_S2[Points+1], Survived[Points+1];
	FILE *input_file, *rho_file, *ps_file;

	if(argc >= 2)
		k = strtof(argv[1], NULL);
	else
		k = 1;

	if(argc >= 3)
		la = strtof(argv[2], NULL);
	else
		la = 2.5;

	if(argc >= 4)
		lb = strtof(argv[3], NULL);
	else
		lb = CritLambda;

	if(argc >= 5) {
		tmax = strtod(argv[4], NULL);
	} else {
		tmax = 1e4;
	}

	if(argc >= 6) {
		Size = (int)strtod(argv[5], NULL);
	} else {
		Size = 1000;
	}

	if(argc >= 7) {
		sim_init = strtod(argv[6], NULL);
	} else {
		sim_init = 0;
	}

	if(argc >= 8) {
		sim_final = strtod(argv[7], NULL);
	} else {
		sim_final = 1000;
	}	

	sprintf(path, "./data_aperiodic");
	sprintf(str, "/k=%d", k);
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
	}

	for(sim = sim_init; sim < sim_final; sim++) {
		sprintf(str, "%s/rho_%d.dat", path, sim);
		input_file=fopen(str, "r");

		inter = 0;
		while((fscanf(input_file, "%lf\t%lf\n", &t, &rho)) != EOF) {
			Rho_Sum[inter] += rho;
			Rho_S2[inter] += rho*rho;
			Time_Sum[inter] += t;
			Time_S2[inter] += t*t;
			Survived[inter] += 1;
			inter++;
		}
		fclose(input_file);
	}

	sprintf(str, "%s/rho_av.dat", path);
	printf("%s", str);
	rho_file=fopen(str, "w+");

	sprintf(str, "%s/surv_prob.dat", path);
	ps_file=fopen(str, "w+");

	i = 0;
	while(i < Points) {
		if(Time_Sum[i] == 0. && Rho_Sum[i] == 0.)
			break;
		Rho_Sum[i] /= (double)(sim_final - sim_init);
		Time_Sum[i] /= (double)(Survived[i]);

		time_desv = sqrt((double)(Time_S2[i] / (sim_final - sim_init) ) - pow(Time_Sum[i], 2.));
		rho_desv = sqrt((double)(Rho_S2[i] / (Survived[i]) ) - pow(Rho_Sum[i], 2.));

		fprintf(rho_file, "%f,%f,%f,%f\n", Time_Sum[i], Rho_Sum[i], time_desv, rho_desv);
		fprintf(ps_file, "%f,%f\n", Time_Sum[i], (double)(Survived[i] / (sim_final - sim_init)));
		i++;
	}

	fclose(rho_file);
	fclose(ps_file);

	return 0;
}
 
