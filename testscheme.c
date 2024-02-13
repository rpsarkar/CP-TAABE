#include "bkem.h"
#include <time.h>
clock_t t,t1,t0,t2,t4,t44; 

     
int main(int argc, const char *argv[]) {

	FILE *param = fopen("a.param", "r");
	char buf[4096];
	fread(buf, 1, 4096, param);
    
    		//printf("\nSystem setup Key\n\n");

	bkem_global_params_t gps;
	setup_global_system(&gps, (const char*) buf, (argc > 1) ? atoi(argv[1]) : 2048);

		printf("Global System parameters: N = %d\n\n", gps->N);

		bkem_system_t sys;
		int w_u[Max_N][LogMax_N], L_0_u[Max_N][LogMax_N], L_1_u[Max_N][LogMax_N];
	
		setup(&sys, gps, w_u, L_0_u, L_1_u);
		
		
		get_enc_key(sys,gps);
        
          	get_decryption_key(gps, sys, sys->PK, w_u, L_0_u, L_1_u);
          	     
        

    }
    
