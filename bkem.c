#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "bkem.h"
#include <time.h>
clock_t setup_time, keygen_time, enc_time, dec_time;

// --------------------------------- Initialization --------------------------------------------
void setup_global_system(bkem_global_params_t *gps, const char *pstr, int N) 
{
	bkem_global_params_t params;
    	params = pbc_malloc(sizeof(struct bkem_global_params_s));
    	params->N = N;
    	pairing_init_set_str(params->pairing, pstr);
    	*gps = params;
}





void generateRandomArrays(int numArrays, int arrayLength, int randomArrays[Max_N][LogMax_N]) {
    srand(time(NULL));
    
    for (int i = 0; i < numArrays; i++) {
        for (int j = 0; j < arrayLength; j++) {
            randomArrays[i][j] = rand() % 2;
        }
    }
}




void setup(bkem_system_t *sys, bkem_global_params_t gps, int w_u[][LogMax_N], int L_0_u[][LogMax_N], int L_1_u[][LogMax_N]) 
{
    	setup_time = clock();
    	bkem_system_t gbs;
    	bkem_secret_key_t sk;
    	gbs = pbc_malloc(sizeof(struct bkem_system_s));
    	gbs->PK = pbc_malloc(sizeof(struct pubkey_s));

// ---------------------------------Choose random generator g_1 --------------------------------------------
    	element_init_G1(gbs->PK->g_1, gps->pairing);
    	element_random(gbs->PK->g_1);

//----------------------------------Choose another generator g_2 -------------------------------------------
    	element_init_G2(gbs->PK->g_2, gps->pairing);
    	element_random(gbs->PK->g_2);

//----------------------------------Choose another generator g_3 -------------------------------------------
    	element_init_G2(gbs->PK->g_3, gps->pairing);
    	element_random(gbs->PK->g_3);

//----------------------------------Choose another generator g_4 -------------------------------------------
    	element_init_G2(gbs->PK->g_4, gps->pairing);
    	element_random(gbs->PK->g_4);

//----------------------------------Choose another generator g_4 -------------------------------------------
    	element_init_G2(gbs->PK->g_5, gps->pairing);
    	element_random(gbs->PK->g_5);

// ---------------------------------random alpha in Zn ----------------------------------------------------
    	element_t alpha;
    	element_init_Zr(alpha, gps->pairing);
    	element_random(alpha);
	//element_printf("The component alpha is = %B\n\n", alpha);

// ---------------------------------random beta in Zn------------------------------------------------------
    	element_t beta;
   	element_init_Zr(beta, gps->pairing);
    	element_random(beta);
	//element_printf("The component beta is = %B\n\n", beta);
	element_init_Zr(gbs->PK->carry_beta, gps->pairing);
	element_set(gbs->PK->carry_beta,beta);

// --------------------------------- computes [beta^i]_1 --------------------------------------------------
	element_t temp1;
	element_init_Zr(temp1, gps->pairing);
	element_set1(temp1);
    	for(int i=1;i<=Max_N; i++)
	{
		element_init_G1(gbs->PK->g_1_beta_i_1[i], gps->pairing);
		element_mul(temp1,temp1,beta);
    		element_pow_zn(gbs->PK->g_1_beta_i_1[i],gbs->PK->g_1,temp1);
    		//element_printf("The component [beta^{%d}]_1 is = %B\n\n", i,gbs->PK->g_1_beta_i_1[i]);
    	} 

// --------------------------------- computes [beta]_2 ---------------------------------------------
	element_init_G1(gbs->PK->g_2_beta_2, gps->pairing);
    	element_pow_zn(gbs->PK->g_2_beta_2,gbs->PK->g_2,beta);
    	//element_printf("The component [beta]_2 is = %B\n\n", gbs->PK->g_2_beta_2);
 
// --------------------------------- computes e(g_1,g_1)^alpha -------------------------------------
	element_init_GT(gbs->PK->e_g_1_g_1, gps->pairing);
	element_init_GT(gbs->PK->e_g_1_g_1_alpha, gps->pairing);
    	pairing_apply(gbs->PK->e_g_1_g_1,gbs->PK->g_1,gbs->PK->g_1,gps->pairing);
	//element_printf("The component e(g_1,g_1)^alpha is = %B\n\n", gbs->PK->e_g_1_g_1_alpha);
	element_pow_zn(gbs->PK->e_g_1_g_1_alpha,gbs->PK->e_g_1_g_1,alpha);
    	//element_printf("The component e(g_1,g_1)^alpha is = %B\n\n", gbs->PK->e_g_1_g_1_alpha);  

// --------------------------------- computes W_u --------------------------------------------------


	//int w_u[Max_N][LogMax_N];
    	generateRandomArrays(Max_N, LogMax_N, w_u);
    	for (int i = 0; i < Max_N; i++) 
	{

    	/*	printf("The w_u value of %d-user: ", i);
    		for (int j = 0; j < LogMax_N; j++) 
		{
        		printf("%d", w_u[i][j]);
    		}
    		printf("\n");
	*/	
		element_init_Zr(gbs->PK->hash[i], gps->pairing);
		element_random(gbs->PK->hash[i]);
		//element_printf("The hash value of the %d-th user is = %B\n\n", i, gbs->PK->e_g_1_g_1_alpha);
    	}



/*	int w_u[Max_N][LogMax_N];
	element_t A_k[Max_N][Attri_no_u];
    	generateRandomArrays(Max_N, LogMax_N, Attri_no_u, w_u, A_k);
    	for (int i = 0; i < Max_N; i++) 
	{
    		printf("The hash value of %d-user: ", i);
    		for (int j = 0; j< LogMax_N; j++) 
		{
        		printf("%d", w_u[i][j]);
    		}
    		printf("\n");
    	}
*/


	
    	setup_time = clock() - setup_time;
    	double time_taken0 = ((double)setup_time)/CLOCKS_PER_SEC; // in seconds 
    	printf("Setup took %f seconds to execute \n\n", time_taken0);  
       

//--------------------------- MPK and MSK generation is done ----------------------------------------------------------------
  
 	
// -------------------------- To compute secret keys of users -----------------------------------------------------------------

 	keygen_time = clock();

	for(int i = 0; i < Max_N; i++)
	{
		element_init_Zr(gbs->PK->hash_plus_beta[i], gps->pairing);
		element_add(gbs->PK->hash_plus_beta[i], gbs->PK->hash[i], beta);
		//element_printf("The value of beta+hash(ID_{%d}) is= %B\n\n", i, gbs->PK->hash_plus_beta[i]);
	}

	
	//for (int i = 0; i < Max_N; i++)  		// skip the key-generation algorithm for all users
	//{	

		int i=7;				// calculate the key-generation algorithm for a single user
		

		int ct=0;
		int ct1=0; 
    		for (int k = 1; k<= Attri_no_u; k++) 
		{	
        		
			element_init_Zr(gbs->PK->A_k[k], gps->pairing);
			element_random(gbs->PK->A_k[k]);
			//element_printf("The attribute A_k[%d] set = %B\n\n", k, gbs->PK->A_k[k]); 
    		}
    		
		for (int j = 0; j< LogMax_N; j++) 
		{
			//printf("%d ", w_u[i][j]);
        		if (w_u[i][j]==0)
			{	
				L_0_u[i][ct]=j;
				//printf("L_0_u[%d][%d]=%d",i,ct,L_0_u[i][ct]);
				ct++;
			}
			else
			{	
				L_1_u[i][ct1]=j;
				//printf("L_1_u[%d][%d]=%d",i,ct1,L_1_u[i][ct1]);
				ct1++;
			}
			//printf("\n");
			
    		}
    		
    /*		printf("; "); 
		for(int j = 0; j < ct; j++) 
		{
			printf("%d ", L_0_u[i][j]);
			 	
		}
		printf("; "); 
		for(int j = 0; j < ct; j++) 
		{		 
			printf("%d ", L_1_u[i][j]);
		}
		printf("\n");
    		//printf("; ct: %d; ", ct);
    		//printf("ct1: %d\n", ct1);
    		
*/	
	
    		for (int b = 0; b<= bit; b++) 
		{
			for (int k = 0; k<= Attri_no_u; k++) 
			{	
        		
				element_init_Zr(gbs->PK->r_b[b][i][k], gps->pairing);
				element_random(gbs->PK->r_b[b][i][k]);
				//element_printf("The random values r^(%d)_[%d][%d] of the user %d is= %B\n\n", b,i, k,i, gbs->PK->r_b[b][i][k]); 
    			}
		}
		
		
		
		
		//---------------------Compute K ---------------------------------------
		element_t temp2;
		element_init_Zr(temp2, gps->pairing);
		//element_printf("The value of beta+hash(ID_{%d}) is= %B\n\n", i, gbs->PK->hash_plus_beta[i]);
		element_invert(temp2, gbs->PK->hash_plus_beta[i]);
		element_init_G1(gbs->K[i], gps->pairing);
		element_pow_zn(gbs->K[i],gbs->PK->g_2,temp2);
		//element_printf("The secret key component K[%d] of the user %d is= %B\n\n", i, i, gbs->K[i]);

		//---------------------Compute K0^(b) and K1^(b)---------------------------------------
		element_t temp3,temp4;
		for (int b = 0; b<= bit; b++) 
		{	
			element_init_G1(gbs->K0[b][i], gps->pairing);
			element_init_G1(temp3, gps->pairing);
			element_pow_zn(temp3,gbs->PK->g_1,alpha);
			element_init_G1(temp4, gps->pairing);
			element_pow_zn(temp4,gbs->PK->g_5,gbs->PK->r_b[b][i][0]);
			element_mul(gbs->K0[b][i],temp3,temp4);
			//element_printf("The secret-key component K0[%d][%d] of the user %d is= %B\n\n", b,i,i, gbs->K0[b][i]); 
			
			element_init_G1(gbs->K1[b][i], gps->pairing);
			element_pow_zn(gbs->K1[b][i],gbs->PK->g_1,gbs->PK->r_b[b][i][0]);
			//element_printf("The secret-key component K1[%d][%d] of the user %d is= %B\n\n", b,i,i, gbs->K1[b][i]);

			
			for (int k = 1; k<= Attri_no_u; k++) 
			{	
				//---------------------Compute K_{i,2}^(b) ---------------------------------------
        			element_init_G1(gbs->K2[b][i][k], gps->pairing);
				element_pow_zn(gbs->K2[b][i][k],gbs->PK->g_1,gbs->PK->r_b[b][i][k]);
				//element_printf("The secret-key component K2[%d][%d][%d] of the user %d is= %B\n\n", b,i,k,i, gbs->K2[b][i][k]);

				//---------------------Compute K_{i,3}^(b) ---------------------------------------
        			element_t temp5, temp6, temp7, temp8, temp9;
				element_init_G1(gbs->K3[b][i][k], gps->pairing);
				element_init_Zr(temp5, gps->pairing);
				element_mul(temp5,gbs->PK->A_k[k],gbs->PK->r_b[b][i][k]);
				element_init_G1(temp6, gps->pairing);
				element_pow_zn(temp6,gbs->PK->g_3,temp5);			// compute [A_i r_i^(b)]_3
				element_init_G1(temp7, gps->pairing);
				element_pow_zn(temp7,gbs->PK->g_2,gbs->PK->r_b[b][i][k]);
				element_init_Zr(temp8, gps->pairing);
				element_neg(temp8, gbs->PK->r_b[b][i][0]);
				element_init_G1(temp9, gps->pairing);
				element_pow_zn(temp9,gbs->PK->g_4,temp8);
				element_mul(gbs->K3[b][i][k],temp6,temp7);
				element_mul(gbs->K3[b][i][k],gbs->K3[b][i][k],temp9);
				//element_printf("The secret-key component K3[%d][%d][%d] of the user %d is= %B\n\n", b,i,k,i, gbs->K3[b][i][k]); 
				
				//---------------------Compute K_{i,4}^(b) ---------------------------------------
        			element_t temp10, temp11, temp12, temp13, sum, sum1;
				element_init_Zr(sum, gps->pairing);
				element_init_Zr(sum1, gps->pairing);
				element_set0(sum);
				element_set0(sum1);
				if(b==0)
				{
					for(int m=0; m<=ct; m++)
					{	
						element_init_Zr(temp10, gps->pairing);
						element_random(temp10);
						element_sub(temp10, beta, temp10);
						element_invert(temp10,temp10);
						element_add(sum,sum,temp10);	
					}
					element_init_G1(gbs->K4[b][i][k], gps->pairing);
					element_pow_zn(gbs->K4[b][i][k],gbs->PK->g_2,sum);
					//element_printf("The secret-key component K4[%d][%d][%d] of the user %d is= %B\n\n", b,i,k,i, gbs->K4[b][i][k]); 
				}
				if(b==1)
				{
					for(int m1=0; m1<=ct1; m1++)
					{
						element_init_Zr(temp12, gps->pairing);
						element_random(temp12);
						element_sub(temp12, beta, temp12);
						element_invert(temp12,temp12);
						element_add(sum1,sum1,temp12);	
					}
					element_init_G1(gbs->K4[b][i][k], gps->pairing);
					element_pow_zn(gbs->K4[b][i][k],gbs->PK->g_2,sum1);
					//element_printf("The secret-key component K4[%d][%d][%d] of the user %d is= %B\n\n", b,i,k,i, gbs->K4[b][i][k]);
				}
				
    			}
		}
			
    	//}
	
   	int size_of_MPK=(2*Max_N+1+1+5) * sizeof(element_t);		
  	element_printf("size_of_MPK = %d in bytes\n\n", size_of_MPK);
	int size_of_MSK=(2) * sizeof(element_t);		
  	element_printf("size_of_MSK = %d in bytes\n\n", size_of_MSK);
	int size_of_SK=(2*(3*Attri_no_u+2)+1) * sizeof(element_t);		
  	element_printf("size_of_SK = %d in bytes\n\n", size_of_SK);
	


    	keygen_time = clock() - keygen_time;
    	double time_taken1 = ((double)keygen_time)/CLOCKS_PER_SEC; // in seconds 
    	printf("KeyGen took %f seconds to execute \n\n", time_taken1); 


    	*sys = gbs;
    	element_clear(alpha);
    	element_clear(beta);  
    
    	//----------------------------Key Gen is done ----------------
}


//=================================Encryption =============================================


void get_enc_key( bkem_system_t gbs, bkem_global_params_t gps) 
{	
	
     	enc_time=clock();
     	//element_t s;
	element_init_Zr(gbs->s, gps->pairing); 
	element_random(gbs->s);
	for(int i=0;i< Max_N;i++)						// Vector s initialization
	{
		if (i==0)
		{
			element_init_Zr(gbs->y_vector[i], gps->pairing);
			element_set(gbs->y_vector[i],gbs->s);
		}
		else
		{
			element_init_Zr(gbs->y_vector[i], gps->pairing);
			element_random(gbs->y_vector[i]);
		}
	}
	
/*	element_t zero,one,neg_one,temp0;
	element_init_Zr(zero, gps->pairing);
	element_set0(zero);
	element_printf("The zero element = %B\n\n", zero);
	element_init_Zr(one, gps->pairing);
	element_set1(one);
	element_printf("The one element = %B\n\n", one);
	element_init_Zr(neg_one, gps->pairing);
	element_neg(neg_one,one);
	element_printf("The negative element of one= %B\n\n", neg_one);
	element_init_Zr(temp0, gps->pairing);
	element_add(temp0,neg_one,one);
	element_printf("The addition of one and negative one= %B\n\n", temp0);
*/
	for(int j=1;j<= Subs_Attri_Num;j++)     						//Access Matrix Initialization 
	{
		element_init_Zr(gbs->lambda_vector[j], gps->pairing);
		element_set0(gbs->lambda_vector[j]);
		for(int i=0;i< Max_N;i++)
		{
			element_t temp20;
			element_init_Zr(temp20, gps->pairing);
			element_init_Zr(gbs->Matrix[j][i], gps->pairing);
			element_random(gbs->Matrix[j][i]);
			element_mul(temp20, gbs->Matrix[j][i],gbs->y_vector[i]);
			element_add(gbs->lambda_vector[j],gbs->lambda_vector[j],temp20);
		}
		
	}

     	element_init_GT(gbs->M, gps->pairing);
     	element_random(gbs->M);								//Assign the broadcasting message M
	//element_printf("The message which I want to broadcast= %B\n\n", gbs->M);					
	element_init_GT(gbs->M0, gps->pairing);
     	element_random(gbs->M0);					
	element_init_GT(gbs->M00, gps->pairing);
	element_add(gbs->M00,gbs->M,gbs->M0); 						//Assign M''= M + M'
	
	
	//------------------------------------Generate C_0^(b) ---------------------------------		
		element_t temp23;
		element_init_GT(temp23, gps->pairing);
		element_pow_zn(temp23,gbs->PK->e_g_1_g_1_alpha,gbs->s);	
     		element_init_GT(gbs->C_0, gps->pairing);
     		element_mul(gbs->C_0,gbs->M0,temp23);
     		//element_printf("The ciphertext component C_0 is= %B\n\n", gbs->C_0);
     		
	
	//------------------------------------Generate C_2^(b) ---------------------------------		
		
     		element_init_G1(gbs->C_2, gps->pairing);
     		element_pow_zn(gbs->C_2,gbs->PK->g_1,gbs->s);
     		//element_printf("The ciphertext component C_2 is= %B\n\n", gbs->C_2_b[b]);
     		
	//------------------------------------Generate t^(b)_d ---------------------------------
	
			
	for(int b=0;b<=bit;b++)
	{
		
		
		
		for(int j=0;j<= Subs_Attri_Num;j++)
		{
			element_init_Zr(gbs->random_t_b[b][j], gps->pairing);
			element_random(gbs->random_t_b[b][j]);
			//element_printf("The random element t^(%d)_{%d} is= %B\n\n", b,j, gbs->random_t_b[b][j]);
		}
	
	
	//------------------------------------ Generate C^(b) ---------------------------------		
		element_t temp21;
		element_init_Zr(gbs->PK->prod_beta_plus_hash, gps->pairing);
		element_init_Zr(temp21, gps->pairing);
		element_set1(gbs->PK->prod_beta_plus_hash);
		for(int j=0;j<Subs_Num;j++)
		{
			element_mul(gbs->PK->prod_beta_plus_hash,gbs->PK->prod_beta_plus_hash,gbs->PK->hash_plus_beta[j]);
		}
		element_mul(temp21,gbs->PK->prod_beta_plus_hash,gbs->s);
		element_mul(temp21,temp21,gbs->random_t_b[b][0]);		
     		element_init_G1(gbs->C_b[b], gps->pairing);
     		element_pow_zn(gbs->C_b[b],gbs->PK->g_1,temp21);
     		//element_printf("The ciphertext component C^(%d) is= %B\n\n", b, gbs->C_b[b]);
   	 
   	
	//------------------------------------ Generate C_1^(b) ---------------------------------		
		element_t temp22;
		element_init_Zr(temp22, gps->pairing);
		element_neg(temp22,gbs->random_t_b[b][0]);
		element_mul(temp22, temp22, gbs->s);		
     		element_init_G1(gbs->C_1_b[b], gps->pairing);
     		element_pow_zn(gbs->C_1_b[b],gbs->PK->g_2,temp21);
     		//element_printf("The ciphertext component C_1^(%d) is= %B\n\n", b, gbs->C_1_b[b]);
     		
     		
     			
		
     		for(int j=1;j<= Subs_Attri_Num;j++)     						
		{
			//------------------------------------ Generate C_{i,3}^(b) ---------------------------------
			
			element_t temp24,temp25;
			element_init_G1(temp24, gps->pairing);
			element_pow_zn(temp24,gbs->PK->g_5,gbs->lambda_vector[j]);
			element_init_G1(temp25, gps->pairing);
			element_pow_zn(temp25,gbs->PK->g_4,gbs->random_t_b[b][j]);
			element_init_G1(gbs->C_3_b[b][j], gps->pairing);
			element_mul(gbs->C_3_b[b][j],temp24,temp25);
			//element_printf("The ciphertext component C_{%d,3}^(%d) is= %B\n\n", j,b, gbs->C_3_b[b][j]);
			
			//------------------------------------ Generate C_{i,4}^(b) ---------------------------------
			
			element_t temp26,temp27,temp28,temp29;
			element_init_Zr(temp26, gps->pairing);
			element_neg(temp26,gbs->PK->A_k[j]);
			element_mul(temp26, temp26,gbs->random_t_b[b][j]);
			element_init_G1(temp27, gps->pairing);
			element_pow_zn(temp27,gbs->PK->g_3,temp26);            // compute [-rho(i).t^(b)_i]_3 
			element_init_Zr(temp28, gps->pairing);
			element_neg(temp28,gbs->random_t_b[b][j]);
			element_init_G1(temp29, gps->pairing);
			element_pow_zn(temp29,gbs->PK->g_2,temp28);
			element_init_G1(gbs->C_4_b[b][j], gps->pairing);
			element_mul(gbs->C_4_b[b][j],temp27,temp29);
			//element_printf("The ciphertext component C_{%d,4}^(%d) is= %B\n\n", j,b, gbs->C_4_b[b][j]);
			
			//------------------------------------ Generate C_{i,5}^(b) ---------------------------------
			
		
			element_init_G1(gbs->C_5_b[b][j], gps->pairing);
			element_pow_zn(gbs->C_5_b[b][j],gbs->PK->g_1,gbs->random_t_b[b][j]);
			//element_printf("The ciphertext component C_{%d,5}^(%d) is= %B\n\n", j,b, gbs->C_5_b[b][j]);
			
			//------------------------------------ Generate C_{i,6}^(b) ---------------------------------
			
			element_t temp30;
			element_init_Zr(gbs->PK->prod_hash_sub_beta[j][b], gps->pairing);
			element_set1(gbs->PK->prod_hash_sub_beta[j][b]);
		
			for(int m=0;m<LogMax_N;m++)
			{
				element_init_Zr(gbs->PK->hash_of_rho_j_b[j][m][b], gps->pairing);	
				element_random(gbs->PK->hash_of_rho_j_b[j][m][b]);
				element_init_Zr(gbs->PK->hash_sub_beta[j][m][b], gps->pairing);
				element_sub(gbs->PK->hash_sub_beta[j][m][b],gbs->PK->carry_beta,gbs->PK->hash_of_rho_j_b[j][m][b]);
				element_mul(gbs->PK->prod_hash_sub_beta[j][b],gbs->PK->prod_hash_sub_beta[j][b],gbs->PK->hash_sub_beta[j][m][b]);
			}
			
			element_init_Zr(temp30, gps->pairing);
			element_mul(temp30,gbs->PK->prod_hash_sub_beta[j][b],gbs->random_t_b[b][0]);
			element_init_G1(gbs->C_6_b[b][j], gps->pairing);
			element_pow_zn(gbs->C_6_b[b][j],gbs->PK->g_1,temp30);
			//element_printf("The ciphertext component C_{%d,6}^(%d) is= %B\n\n", j,b, gbs->C_6_b[b][j]);
			
			
			//------------------------------------ Generate C_{i,7}^(b) ---------------------------------
			
			element_t temp31;
			element_init_Zr(temp31, gps->pairing);
			element_mul(temp31, gbs->random_t_b[b][0], gbs->PK->hash_sub_beta[j][tau][b]);
			element_init_G1(gbs->C_7_b[b][j], gps->pairing);
			element_pow_zn(gbs->C_7_b[b][j],gbs->PK->g_2,temp31);
			//element_printf("The ciphertext component C_{%d,7}^(%d) is= %B\n\n", j,b, gbs->C_7_b[b][j]);
			
			
			//------------------------------------ Generate C_{i,8}^(b) ---------------------------------
			
			element_t temp32,temp33,temp34,temp35,temp36,temp37;
			element_init_Zr(gbs->PK->prod_hash_sub_beta_div_tau[j][b], gps->pairing);
			element_div(gbs->PK->prod_hash_sub_beta_div_tau[j][b], gbs->PK->prod_hash_sub_beta[j][b], gbs->PK->hash_sub_beta[j][tau][b]);
			element_init_G1(temp32, gps->pairing);
			element_pow_zn(temp32,gbs->PK->g_1,gbs->PK->prod_hash_sub_beta_div_tau[j][b]);
			element_init_G1(temp33, gps->pairing);
			element_pow_zn(temp33,gbs->PK->g_2,gbs->random_t_b[b][0]);
			element_init_GT(temp36, gps->pairing);
			pairing_apply(temp36,temp32, temp33,gps->pairing);
			
			element_init_G1(temp34, gps->pairing);
			element_pow_zn(temp34,gbs->PK->g_1, gbs->s);
			element_init_GT(temp35, gps->pairing);
			pairing_apply(temp35,temp34, temp33,gps->pairing);
			element_init_GT(gbs->temp_C_i_8[b][j], gps->pairing);
			element_random(gbs->temp_C_i_8[b][j]);
			element_init_GT(gbs->C_8_b[b][j], gps->pairing);
			element_add(gbs->C_8_b[b][j],gbs->temp_C_i_8[b][j],gbs->M00);
			//element_printf("The ciphertext component C_{%d,8}^(%d) is= %B\n\n", j,b, gbs->C_8_b[b][j]);
			
			
		}	
     		
     		
     		
     		
   	} 
   	

	//int size_of_CT=(2*(3*Attri_no_u+2)+1) * sizeof(element_t);		
  	//element_printf("size_of_SK = %d in bytes\n\n", size_of_CT);


   	enc_time = clock() - enc_time; 
	double time_taken2 = ((double)enc_time)/(CLOCKS_PER_SEC); // in seconds 
  	printf("Encryption algorithm took %f seconds to execute for %d users \n\n", time_taken2, Subs_Num); 
    
       //element_clear(s);         	
}





void get_decryption_key(bkem_global_params_t gps, bkem_system_t gbs, pubkey_t PK, int w_u[][LogMax_N], int L_0_u[][LogMax_N], int L_1_u[][LogMax_N])
 {
 	dec_time=clock();
 	int i=7;
 	int b=0;
 	
 	
	element_t temp63,temp65,temp66;
 	for(int j=1;j<=Subs_Attri_Num;j++)
 	{
 		element_t temp60,temp61,temp62,temp64;
 		element_init_GT(temp60, gps->pairing);
 		pairing_apply(temp60,gbs->C_3_b[b][j],gbs->K1[b][i],gps->pairing);    		//compute e(C_{i,3}^(b), K_1^(b))
 		element_init_GT(temp61, gps->pairing);
 		pairing_apply(temp61,gbs->C_4_b[b][j],gbs->K2[b][i][j],gps->pairing);    	//compute e(C_{i,4}^(b), K_{2,i}^(b))
 		element_init_GT(temp62, gps->pairing);
 		pairing_apply(temp62,gbs->C_5_b[b][j],gbs->K3[b][i][j],gps->pairing);    	//compute e(C_{i,5}^(b), K_{3,i}^(b))
 		element_init_GT(temp63, gps->pairing);
 		element_mul(temp63,temp60,temp61);
 		element_mul(temp63,temp63,temp62);
		element_init_Zr(temp64,gps->pairing);
		element_random(temp64);
		element_pow_zn(temp63,temp63,temp64);
 		
 	}
	//element_printf("The temp63 is= %B\n\n", temp63);
	element_init_GT(temp65, gps->pairing);
	pairing_apply(temp65,gbs->C_2,gbs->K0[b][i],gps->pairing);


	element_init_GT(temp66, gps->pairing);
	element_div(temp66,temp63,temp65);
	//element_printf("The value of M'' is= %B\n\n", temp66);
	//element_printf("The value of M'' is= %B\n\n", gbs->M0);
	
	
	
 	
 	
 	//-----------------------compute F(x) and P(x) ----------------------------------------
 	
 	element_t temp40,temp_sum,temp41,temp42;
 	element_init_Zr(temp40, gps->pairing);
	element_init_Zr(temp_sum, gps->pairing);
	element_set0(temp_sum);
	
	for (int j=0;j<LogMax_N;j++)						
	{
		
		element_invert(temp40,gbs->PK->hash_sub_beta[mu][j][b]);
		element_add(temp_sum,temp_sum,temp40);
	}
	element_init_Zr(gbs->func_beta, gps->pairing);
	element_mul(gbs->func_beta, gbs->PK->prod_hash_sub_beta[mu][b], temp_sum);			// Defining F(x) at x=beta
	
	
	element_init_Zr(temp41, gps->pairing);
	element_div(temp41,gbs->PK->prod_beta_plus_hash,gbs->PK->hash_plus_beta[i]);
	element_init_Zr(gbs->PK->prod_hash, gps->pairing);
	element_set1(gbs->PK->prod_hash);
	for(int j=0;j<Subs_Num;j++)     						
	{
		element_mul(gbs->PK->prod_hash, gbs->PK->prod_hash, gbs->PK->hash[j]);
	}
	
	element_init_Zr(temp42, gps->pairing);
	element_div(temp42,gbs->PK->prod_hash,gbs->PK->hash[i]);
	
	element_init_Zr(gbs->poly_beta, gps->pairing);
	element_sub(gbs->poly_beta, temp41 ,temp42);			// Defining P(x) at x=beta
	
	
	//-----------------compute M'' ----------------------------------
	
	element_t temp43,temp44,temp45,temp46,temp47,temp48,temp49;
	element_init_Zr(temp43, gps->pairing);
	element_pow_zn(temp43,gbs->PK->g_1,gbs->func_beta);
	element_init_GT(temp44, gps->pairing);
	pairing_apply(temp44,temp43,gbs->PK->g_2,gps->pairing);    	//compute e([f(beta]_1,g_2)
	element_init_GT(temp45, gps->pairing);
	element_pow_zn(temp45,temp44,gbs->random_t_b[b][0]);		//compute e([f(beta]_1,g_2)^{t^(b)}
	
	
	element_init_Zr(temp46, gps->pairing);
	element_div(temp46,gbs->PK->prod_hash_sub_beta[mu][b],gbs->PK->hash_sub_beta[mu][tau][b]);
	element_sub(temp46, temp46, gbs->func_beta);
	element_init_G1(temp47, gps->pairing);
	element_pow_zn(temp47,gbs->PK->g_1,temp46);
	element_init_GT(temp48, gps->pairing);
	pairing_apply(temp48,temp47,gbs->PK->g_2,gps->pairing);    	//compute e([prod- f(beta]_1,g_2)
	element_init_GT(temp49, gps->pairing);
	element_pow_zn(temp49,temp48,gbs->random_t_b[b][0]);		//compute e([prod-f(beta]_1,g_2)^{t^(b)}
	
	
	element_t temp50,temp51,temp52,temp53,temp54,temp55,temp56;
	element_init_Zr(temp50, gps->pairing);
	element_mul(temp50,gbs->random_t_b[b][0],gbs->s);
	element_init_Zr(temp51, gps->pairing);
	element_sub(temp51,temp41,gbs->poly_beta);
	element_mul(temp51,temp51,temp50);
	element_init_GT(temp52, gps->pairing);
	pairing_apply(temp52,gbs->PK->g_1,gbs->PK->g_2,gps->pairing); 
	element_pow_zn(temp52,temp52,temp51);				//compute e(g_1,g_2)^{t^(b)s (-p(beta)+ prod)}
	element_init_Zr(temp53, gps->pairing);
	element_invert(temp53,temp42);
	element_init_GT(temp54, gps->pairing);
	element_pow_zn(temp54,temp52,temp53);				//compute (e(g_1,g_2)^{t^(b)s (-p(beta)+ prod)})^{prod}
	element_init_GT(temp55, gps->pairing);
	//element_printf("The value of M'' is= %B\n\n", gbs->temp_C_i_8[b][mu]);
	//element_sub(temp55,gbs->C_8_b[b][mu],gbs->temp_C_i_8[b][mu]);	// recover M''
	//element_printf("The value of M'' is= %B\n\n", temp55);
	//element_printf("The value of M'' is= %B\n\n", gbs->M00);
	
	
	
	
	



 	
/*	for(int j = 0; j < LogMax_N; j++)
	{
		printf(" %d ", w_u[i][j]);
	}
	printf("; ");
	for(int j = 0; j < LogMax_N; j++)
	{
		printf("%d ", L_0_u[i][j]);
	}
	printf("; ");
	for(int j = 0; j < LogMax_N; j++)
	{
		printf("%d ", L_1_u[i][j]);
	}
	printf("\n");
		
 	printf("ct=%d\n",ct);
 	printf("ct1=%d\n",ct1);
 */

		
 	dec_time = clock() - dec_time; 
   	double time_taken3 = ((double)dec_time)/(CLOCKS_PER_SEC); // in seconds 
  	printf("Decryption algorithm took %f seconds to execute for a subscribed user\n\n", time_taken3);
  	
  	//element_printf("The plaintext is %B\n\n",gbs->M);
 	//element_printf("The recover message is %B\n\n",Guess_M); 
 
}





/*void free_global_params(bkem_global_params_t gbs) {
    if (!gbs)
        return;

    pairing_clear(gbs->pairing);
    free(gbs);
}
*/

/*void free_pubkey(pubkey_t pk, bkem_global_params_t gbs) {
    if (!pk)
        return;

    element_clear(pk->g);

    int i;
    for (i = 0; i <= gbs->N; ++i) {
        element_clear(pk->g_i[i]);
    }

    //for (i = 0; i < gbs->A; ++i) {
       // element_clear(pk->v_i[0]);
    //}

}
*/
/*void free_bkem_system(bkem_system_t sys, bkem_global_params_t gbs) {
    if (!sys)
        return;

    free_pubkey(sys->PK, gbs);

    int i;
    /*for (i = 0; i < gbs->N; ++i) {
        element_clear(sys->d_i[i]);
    }*/
//}
