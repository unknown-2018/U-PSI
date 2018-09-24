
/*

- a test that runs both update and PSI computation in U-PSI protocol.

*/
//**********************************************************************

#include "Client.h"

//**********************************************************************
// - Function description: generates a set of random bigintegers,
// and ensures that the values are smaller than the public moduli and unequal to x-coordinates.


bigint* gen_randSet (int size, int max_bitsize, bigint* pubModuli, bigint* x_points, int xpoint_size){

	Random rd;
	mpz_t *pr_val;
	pr_val=(mpz_t*)malloc(size * sizeof(mpz_t));
	int max_bytesize = max_bitsize;
	gmp_randstate_t rand;
	bigint ran;
	rd.init_rand3(rand, ran, max_bytesize);
	bigint temp;
	mpz_init(temp);
	for(int i = 0; i < size; i++){
		mpz_urandomb(temp, rand, max_bitsize); // The last argument is in bit.
		while (mpz_cmp(temp, pubModuli[0]) > 0){ // ensures the elements are smaller than the public moduli.
			mpz_urandomb(temp, rand, max_bitsize);
			for(int j = 0; j < xpoint_size; j++){ //checks the random element is not equal to any x_points.
				if(mpz_cmp(temp, x_points[j]) == 0){
					mpz_urandomb(temp, rand, max_bitsize);
				}
			}
		}
		mpz_init_set(pr_val[i],temp);
	}
	return pr_val;
}

int main(){

	int pub_mod_bitsize = 40; // field bit-size.
	int max_setsize = 1024;
	int table_length = 30;
	// (max_setsize,table_length) can be set as follows: (1024, 30), (2048, 61), (4096, 122), (8192, 245), (16384, 491), (32768, 983), (65536, 2621), (131072, 5242),(262144, 10485), (524288, 20971), (1048576, 41943)
	int bucket_max_load = 100;
	int interSec_size = 15; // an arbitrary choice and can be anything smaller than, equal to max_setsize.
	int xsize = 201;// Note that the number of x is determined by bucket_max_load
	if(xsize < (2 * bucket_max_load) + 1) {cout<<"\nxsize must be greater than 2*bucket_max_load)+1, reset it\n";return 0;}
	int number_of_experiments = 1;
	for(int l = 0;l<number_of_experiments;l++){
		Server serv(xsize, 2, pub_mod_bitsize, max_setsize , bucket_max_load, table_length);
		Server * serv_ptr (& serv);
		int elem_bit_size = 40;
	// Assigning random values to two sets a and b.
	cout<<"\n----------------------------------------------------------------\n";
	cout<<"\t Set_size:       "<<max_setsize<<endl;
	cout<<"\t Pub_mod_bitsize:  "<<pub_mod_bitsize<<endl;
	cout<<"\t Bucket_max_load:  "<<bucket_max_load<<endl;
	cout<<"\t Table_length:     "<<table_length<<endl;
	cout<<"\n----------------------------------------------------------------\n";
	int t1,t2;
	double temp_req = 0;
	double temp_grant = 0;
	double temp_res = 0;
	double temp_intersect = 0;
	double temp_out = 0;
	mpz_t *aa,*bb;
	aa = (mpz_t*)malloc(max_setsize * sizeof(mpz_t));
	bb = (mpz_t*)malloc(max_setsize * sizeof(mpz_t));
	aa = gen_randSet (max_setsize, elem_bit_size,serv.get_pubModuli(), serv.get_xpoints(t1), xsize);
	bb = gen_randSet (max_setsize, elem_bit_size,serv.get_pubModuli(), serv.get_xpoints(t2), xsize);

	cout<<"\n----------------------------------------------------------------\n";
	for(int i = 0; i < interSec_size; i++){
		mpz_set(bb[i], aa[i]);
		cout<<"\t\t\n Elements initially in both sets:"<<bb[i]<<endl;
	}
	cout<<"\n----------------------------------------------------------------\n";
	Client B(serv_ptr, bb, max_setsize);
	string b_id="B_ID";
	Client A(serv_ptr, aa, max_setsize);
	string a_id = "A_ID";
	bigint label;
	cout<<"\n----------------- Client B outsourcing -----------------"<<endl;
	double start_out = clock();
	B.outsource_db(b_id);
	double end_out = clock();
	cout<<"\n----------------- Client A outsourcing -----------------"<<endl;
	temp_out += end_out - start_out;
	A.outsource_db(a_id);
	cout<<"\n------------------------"<<endl;
	cout<<"\n-------- Add up all the values appear below to get the total time ------------"<<endl;
	bigint ins;
	mpz_init_set_str(ins, "1", 10);
	string d = B.update(ins, "insertion", label, "B_ID");
	cout<<"\n------------------------"<<endl;
	cout<<"\n Client B's Update Status:"<<d<<endl;
	cout<<"\n------------------------"<<endl;
		
	//-----------Set intersection------------
	bigint B_tk, **q;
	int* sz;
	cout<<"\n---- Gennerate the Computation Request"<<endl;
	double start_req = clock();
	CompPerm_Request*req = B.gen_compPerm_req(B_tk);
	double end_req = clock();
	temp_req += end_req - start_req;
	cout<<"\n---- Grant the Computation"<<endl;
	double start_grant = clock();
	GrantComp_Info* ptr1 = A.grant_comp(req, q, true);
	double end_grant = clock();
	temp_grant += end_grant - start_grant;
	cout<<"\n---- Server-side Result Computation"<<endl;
	double start_res=clock();
	Server_Result * res = serv.compute_result(ptr1, B_tk);
	double end_res=clock();
	temp_res += end_res - start_res;
	cout<<"\n---- Client-side Result Retirieval"<<endl;
	double start_intersect = clock();
	vector<string>  final_res = B.find_intersection(res, sz, q);
	double end_intersect = clock();
	temp_intersect += end_intersect - start_intersect;
	cout<<"\n\n\t======= Result ======="<<endl;
	for(int i = 0; i < final_res.size(); i++){
		cout<<"\n\nFinal_res "<<i + 1<<": "<<final_res[i]<<endl;
	}
	cout<<"\n\n\t==== Run time ======"<<endl;
	double out = temp_out /number_of_experiments;
	float out_time = out / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Outsourcing-- time:"<<out_time<<endl;
	double com_req = temp_req / number_of_experiments;
	float req_time = com_req / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Computation Request-- time:"<<req_time<<endl;
	double grant = temp_grant / number_of_experiments;
	float grant_time = grant / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Computation Grant-- time:"<<grant_time<<endl;
	double res_= temp_res / number_of_experiments;
	float res_time = res_ / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Server Computation-- time:"<<res_time<<endl;
	double inter = temp_intersect / number_of_experiments;
	float inter_time = inter / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Find intersection-- time:"<<inter_time<<endl;
	cout<<"\n\n\t=================="<<endl;
	//-----------End of Set intersection------------
	}
	return 0;
}
//**********************************************************************
