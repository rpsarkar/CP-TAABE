/****** This is the implementation for Ciphertext-policy attribute-based broadcast encryption*****/

#ifndef H_BKEM
#define H_BKEM

#include <string.h>
#include <pbc/pbc.h>


#define Max_N 2048		// Total Number of users (N)
#define LogMax_N 11		// Length of identity string or l= Log Max_n (base 2) (l)
#define Subs_Num 256		// Number of subscribed users(L)
#define Subs_Attri_Num 64	// maximum size of attribute set associated in encryption (d)
#define Attri_no_u 128		// Number of attributes associated with users(k)
#define bit 1			// bits {0,1} (b)
#define tau 3			// index tau (\tau) 
#define mu 4 			// index mu (\mu) described in Decryption

typedef struct bkem_global_params_s {
	pairing_t pairing;
	int N;
	
}* bkem_global_params_t;

typedef struct pubkey_s 
{
    	element_t g_1; 							// 1st generator: g_1
    	element_t g_2; 							// random element g_2
	element_t g_3; 							// random element g_3
	element_t g_4; 							// random element g_4
	element_t g_5; 							// random element g_5
	element_t carry_beta;						// for carry forward the value of beta
	element_t g_2_beta_2;						// the element for computing [beta]_2
	element_t g_1_beta_i_1[Max_N];					// N elements for computing [beta^i]_1 
	element_t e_g_1_g_1;						// The element for computing e(g_1,g_1)
	element_t e_g_1_g_1_alpha;					// The element for computing e(g_1,g_1)^alpha
	//element_t w_u[Max_N][LogMax_N];				// The element for computing w_u
	element_t A_k[Attri_no_u];					// Attributes for the user u
	element_t r_b[bit][Max_N][Attri_no_u];				// random values of r for calculating secret keys	
	element_t hash[Max_N];						// hash values of users 
	element_t prod_hash;						// hash values of users 
	element_t hash_plus_beta[Max_N];				// Values (beta+H(ID_u)
	element_t prod_beta_plus_hash;					// Product of values of subscribed (beta+H(ID_u)
	element_t hash_of_rho_j_b[Subs_Attri_Num][LogMax_N][bit];	// hash values of (rho(i)||j||b) 
	element_t hash_sub_beta[Subs_Attri_Num][LogMax_N][bit];		// Values (beta-H(rho(i)||j||b)
	element_t prod_hash_sub_beta[Subs_Attri_Num][bit];		// Values of product of all (beta-H(rho(i)||j||b)
	element_t prod_hash_sub_beta_div_tau[Subs_Attri_Num][bit];	// Values of product of all (beta-H(rho(i)||j||b)
	
	


     
}* pubkey_t;

typedef struct bkem_secret_key_s 
{
	/** Private key of user s */
	element_t d[Max_N][2];
}* bkem_secret_key_t;

	
typedef struct
{
    int useri,channelj;  // User user_i  subscribe to chanel channel_j
} ID;



typedef struct bkem_system_s {
	pubkey_t PK;
	
	element_t K[Max_N];
	element_t K0[bit][Max_N];
	element_t K1[bit][Max_N];
	element_t K2[bit][Max_N][Attri_no_u];
	element_t K3[bit][Max_N][Attri_no_u];
	element_t K4[bit][Max_N][Attri_no_u];
	
	
	
	element_t Matrix[Subs_Attri_Num][Max_N];
	element_t y_vector[Max_N];
	element_t lambda_vector[Subs_Attri_Num];
	element_t M;
	element_t M0;
	element_t M00;
	element_t random_t_b[bit][Subs_Attri_Num];
	
	
	element_t s;
	element_t C_b[bit];
	element_t C_0;
	element_t C_1_b[bit];
	element_t C_2;
	element_t C_3_b[bit][Subs_Attri_Num];
	element_t C_4_b[bit][Subs_Attri_Num];
	element_t C_5_b[bit][Subs_Attri_Num];
	element_t C_6_b[bit][Subs_Attri_Num];
	element_t C_7_b[bit][Subs_Attri_Num];
	element_t temp_C_i_8[bit][Subs_Attri_Num];
	element_t C_8_b[bit][Subs_Attri_Num];
	

	element_t func_beta;
	element_t poly_beta;
	
	
	
		
}* bkem_system_t;




typedef struct header_s 
{
    element_t Matrix[Subs_Attri_Num][Max_N];	// Access Matrix M from Z_p^{d*N}


}* header_t;


typedef struct kpair_s {
    element_t *HR;
}* kpair_t;

typedef struct keypair_s {
    element_t *HDR;
    element_t K;
}* keypair_t;


void free_pubkey(pubkey_t pk, bkem_global_params_t gbs);
void free_bkem_system(bkem_system_t sys, bkem_global_params_t gbs);
void free_global_params(bkem_global_params_t gbs);




void setup_global_system(bkem_global_params_t *gps, const char *params, int n);
void setup(bkem_system_t *sys, bkem_global_params_t gps, int w_u[][LogMax_N], int L_0_u[][LogMax_N], int L_1_u[][LogMax_N]);
void get_enc_key(bkem_system_t sys, bkem_global_params_t gps);
void get_decryption_key(bkem_global_params_t gbs, bkem_system_t sys, pubkey_t PK, int w_u[][LogMax_N],  int L_0_u[][LogMax_N], int L_1_u[][LogMax_N]);


#endif
