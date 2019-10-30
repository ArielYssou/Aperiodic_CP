#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#pragma GCC optimize("no-move-loop-invariants", "no-unroll-loops")

/*****	Gerador de Numeros Aleatorios	*****/
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


static uint64_t s[2];

uint64_t next(void) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(void) {
	static const uint64_t JUMP[] = { 0xdf900294d8f554a5, 0x170865df4b3201fc };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			next();
		}

	s[0] = s0;
	s[1] = s1;
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump(void) {
	static const uint64_t LONG_JUMP[] = { 0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			next();
		}
	s[0] = s0;
	s[1] = s1;
}

/*seed generator */
uint64_t xs; /* The state can be seeded with any value. */

uint64_t nexts() {
	uint64_t zs = (xs += UINT64_C(0x9E3779B97F4A7C15));
	zs = (zs ^ (zs >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	zs = (zs ^ (zs >> 27)) * UINT64_C(0x94D049BB133111EB);
	return zs ^ (zs >> 31);
}

//Essa funcao leva de uint64_t para um double em [0,1]
static inline double to_double(uint64_t xd) {
       const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | xd >> 12 };
       return u.d - 1.0;
    }


#define UNUSED(expr) do { (void)(expr); } while (0)

#define CritLambda 3.29785
#define Points 500

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
}

void InitializeArrays(int Lattice[], int Infected[], int Vizd[], int Vize[], int Size, int* list_last, double rho_0) {
	int i, InitialPop, pos, site, sites[Size];
	for(i = 0; i < Size; i++) {
		Lattice[i] = 0;
		Infected[i] = 0;
		Vizd[i] = i + 1;
		Vize[i] = i - 1;
	}
	Vizd[Size-1] = 0;
	Vize[0]=Size - 1;

	for(i = 0; i < Size; i++) {
		sites[i] = i;
	}
	UNUSED(sites);
	UNUSED(site);

	if(rho_0 > 0) {
	InitialPop = (int)(Size*rho_0);
		*list_last = 0;
		while(*list_last < InitialPop) {
			pos = rand()%(Size - *list_last);
			site = sites[pos];
			Lattice[site] = 1;
			Infected[*list_last] = site;
			*list_last += 1;
			sites[pos] = sites[Size - *list_last];
		}
		*list_last -= 1;
	}
	else {
		Lattice[(int)Size/2] = 1;
		Infected[0] = (int)Size/2;
		*list_last = 0;
	}
}

int FibbTemp[120000000];

void fibonacci_word(int k, int Size, int Fibb[], int offset) {
	int i, j, pos, end;
	
	if( k == 0) {
		for(i = 0; i< Size; i++) 
			Fibb[i] = 1;
	} else {
		FibbTemp[0] = 0;
		for(i=1; i <= k; i++) {
			FibbTemp[i] = 1;
		}
		pos = 1;
		end = k + 1;
		while(end < Size  + offset) {
			if(FibbTemp[pos] == 0){
				FibbTemp[end] = 0;
				end++;
				for(i=0; i < k && end < Size + offset; i++) {
					FibbTemp[end] = 1;
					end++;
				}
			}
			else {
				FibbTemp[end] = 0;
				end++;
			}
			pos++;
		}

		j = 0;
		for(i = offset; i < Size + offset; i++) {
			Fibb[j] = FibbTemp[i];
			j++;
		}
	}
}

int main(int argc, char **argv) {
	char path[1000], str[1000], *node;
	int k, point,  pos, list_pos, list_last,
			sim, sim_init, sim_final, Size, offset, fibb_index;

	int *Lattice, *InfecList, *Vizd, *Vize, *Fibb;

	double la, lb, prob, pleft, pright, pla, pra, plb,
				 prb, t, tmax, delta_t, tinter, rho, rho_0,
				 time_spent, Times[Points];

	FILE *output_file;

	clock_t begin = clock();


	// INPUT
	k = strtod(argv[1], NULL);
	rho_0 = strtof(argv[2], NULL);
	delta_t = strtof(argv[3], NULL);

	la = strtof(argv[4], NULL);
	lb = strtof(argv[5], NULL);

	tmax = strtod(argv[6], NULL);

	Size = (int)strtod(argv[7], NULL);
	sim_init = strtod(argv[8], NULL);
	sim_final = strtod(argv[9], NULL);
	node = argv[10];

	LogSpaced(Times, 1, tmax);

	Lattice = (int*)calloc(Size, sizeof(int));
	InfecList = (int*)calloc(Size, sizeof(int));
	Fibb = (int*)calloc((int)tmax, sizeof(int));
	Vize = (int*)calloc(Size, sizeof(int));
	Vizd = (int*)calloc(Size, sizeof(int));

	srand((unsigned)time(NULL));
	offset = rand()%(int)(10 * tmax);
	//offset = 0;
	fibonacci_word(k, tmax, Fibb, offset);


	pla = (la / 2.0) / (1.0 + la);
	pra = la / (1.0 + la);
	plb = (lb / 2.0) / (1.0 + lb);
	prb = lb / (1.0 + lb);

//	if(Fibb[0] == 0) { 
//		pleft = pla;
//		pright = pra;
//	}
//	else {
//		pleft =  plb;
//		pright =  prb;
//	}
//

	pleft =  plb;
	pright =  prb;

	sprintf(path, "/home/ariel/Simulacoes/data_aperiodic");
	//sprintf(path, "./data_aperiodic");
	mkdir(path, 0777);

	sprintf(str, "/k=%d", k);
	strcat(path, str);
	mkdir(path, 0777);

	sprintf(str, "/Delta_t=%.2f", delta_t);
	strcat(path, str);
	mkdir(path, 0777);

	sprintf(str, "/Lambda_a=%.3f", la);
	strcat(path, str);
	mkdir(path, 0777);

	sprintf(str, "/Lambda_b=%.8f", lb);
	strcat(path, str);
	mkdir(path, 0777);

	sprintf(str, "/Tmax=%e", tmax);
	strcat(path, str);
	mkdir(path, 0777);

	sprintf(str, "/Size=%d", (int)Size);
	strcat(path, str);
	mkdir(path, 0777);

	xs = rand();
	s[1] = nexts();
	s[2] = nexts();

	for(sim = sim_init; sim <= sim_final; sim++) {
		sprintf(str, "%s/rho_%d.dat", path, sim);
		output_file=fopen(str, "w+");

		offset = rand()%(int)(2 * tmax);
		fibonacci_word(k, tmax, Fibb, offset);

		InitializeArrays(Lattice, InfecList, Vizd, Vize, Size, &list_last, rho_0);

		point = 0;
		fibb_index = 0;
		t = 0;
		tinter = delta_t;

		rho = (list_last + 1.0) / (double)Size;
		fprintf(output_file, "%f\t%f\n", t, rho);
		while(t < tmax && list_last >= 0) {
			pos = (next() % (list_last+1));
			list_pos = InfecList[pos];

			prob = to_double(next());
			
			if(prob < pleft) {
				if(Lattice[Vize[list_pos]] == 0) {
					Lattice[Vize[list_pos]] = 1;
					list_last++;
					InfecList[list_last] = Vize[list_pos];
				}
			}
			else if (prob < pright) {
				if(Lattice[Vizd[list_pos]] == 0) {
					Lattice[Vizd[list_pos]] = 1;
					list_last++;
					InfecList[list_last] = Vizd[list_pos];
				}
			}
			else {
				Lattice[list_pos] = 0;
				InfecList[pos] = InfecList[list_last];
				InfecList[list_last] = -1;
				list_last -= 1;
				if(list_last < 0) {
					fprintf(output_file, "%f\t%f\n", t, 0.);
					break;
				}
			}

			t += 1.0 / (list_last + 1.0);

			if(t >= tinter) {
				fibb_index++;
				tinter += delta_t;

				if(Fibb[fibb_index] == 0) { 
					pleft = pla;
					pright = pra;
				}
				else {
					pleft =  plb;
					pright =  prb;
				}
			}

			if(t >= Times[point]) {
				rho = (double)(list_last + 1) / Size;
				fprintf(output_file, "%f\t%f\n", t, rho);
				point++;
				while(t >= Times[point]) {
					point++;
				}
			}
		}
		if( list_last < 0 )
				fprintf(output_file, "%f\t0\n", t);
		fclose(output_file);
	}

	free( Lattice );
	free( InfecList );
	free( Fibb );
	free( Vizd );
	free( Vize );

	clock_t end = clock();

	time_spent = (double)(end-begin)/CLOCKS_PER_SEC;

	printf("%lf\n", time_spent);

	sprintf(str, "bash /home/ariel/aurora.sh -r %s", node);
	system(str);

	return 0;
}
