/*

- Server-side computation of the U-PSI protocol

*/

#include"Server.h"

//**********************************************************************
// - Description: Constructor- It generates a set of x-coordinates and public moduli. Also, it sets a hash table parameters.

Server::Server(int num_xpoints, int dbs_size, int pub_mod_bitsize, int maxSetsize, int NoEl_bucket, int tb_size){

	max_setsize = maxSetsize;
	NoElem_in_bucket = NoEl_bucket;
	pub_moduli_bitsize = pub_mod_bitsize;
	xpoint_size = num_xpoints;
	Random rd;
	xpoints = rd.gen_randSet(num_xpoints, pub_moduli_bitsize - 20); // 20 is an arbitrary choioce. The main requirement is that x-coordinates must be non-zero.
	pu_moduli = rd.gen_randSet(1, pub_moduli_bitsize);
	mpz_nextprime(pu_moduli[0], pu_moduli[0]);
	db_size = dbs_size;
	count = 0;
	db = new Client_Dataset[dbs_size];
	table_size = tb_size;
}
//**********************************************************************
// - Function description: returns x_coordinates.

bigint* Server::get_xpoints(int& size){

	size = xpoint_size;
	bigint *ptr = xpoints;
	return ptr;
}
//**********************************************************************
// - Function description: returns public moduli.

bigint*  Server::get_pubModuli(){

	bigint *ptr;
	ptr = (mpz_t*)malloc(1 * sizeof(mpz_t));
	ptr = pu_moduli;
	return ptr;
}
//**********************************************************************
// - Function description: returns the upper bound on the set size.

int Server::get_maxSetsize(){

	return max_setsize;
}
//**********************************************************************
// - Function description: returns bin's capacity: d, as part of the hashtable parameters.

int Server::get_NoElem_in_bucket(){

	return NoElem_in_bucket;
}
//**********************************************************************
// - Function description: stores a client's dataset, given the avaialbe index at the server-side database.
// It's called by store_poly().

void Server::set_db(int index, Client_Dataset &p){
	db[index] = p;
}
//**********************************************************************
// - Function description: stores a client's dataset at the server.

void Server::store_poly(Client_Dataset& p){
	if(count <= db_size){
		if(p.poly[0].get_poly_ID() == "B_ID")
			set_db(0, p);
		else{
			++count;
			set_db(count, p);
		}
	}
	else {
		cout<<"\n Error: No space to store anymore Dataset."<<endl;
		return;
	}
}
//**********************************************************************
// - Function description: generates two pseudorandom polynomials and a set of pseudorandom values (for each bin's of a hashtable).

void Server::regen_PRpolys(bigint key_, bigint **&w_A_, bigint **& w_B_, bigint **&a_, bigint**& tmp_bl_, bigint tmp_key_){

	bigint **w_A, **w_B, **a, **tmp_bl, temp_key, *Sw1, *Sw2, **derived_key, tmp_bf_A, tmp_bf_B, temp_key2;
	w_A_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	w_B_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	a_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	Sw1 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	Sw2 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	derived_key = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	tmp_bl_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	// regenerates w_1, w_2, and a_i
	gmp_randstate_t randC, rand_Pas_C, rand02, rand03, rand04, randD, rand_Pas_D;
	gmp_randinit_default(rand02);
	gmp_randinit_default(rand03);
	gmp_randinit_default(rand04);
	gmp_randinit_default(randC);
	gmp_randinit_default(randD);
	gmp_randinit_default(rand_Pas_C);
	gmp_randinit_default(rand_Pas_D);
	gmp_randseed(randD, tmp_key_);
	gmp_randseed(randC, key_);
	mpz_init(temp_key);
	mpz_init(temp_key2);
	for(int i = 0; i < table_size; i++){
		w_A_[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		w_B_[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		a_[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		tmp_bl_[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_urandomb(temp_key, randC, pub_moduli_bitsize);// generates a key: k_i
		gmp_randseed(rand_Pas_C, temp_key);
		mpz_urandomb(temp_key2, randD, pub_moduli_bitsize);
		gmp_randseed(rand_Pas_D, temp_key2);
		derived_key[i] = (mpz_t*)malloc(3 * sizeof(mpz_t));
		// derives 3 keys from k_i. key 1 is used to generate a_i, while key 2 and 3 used to genetare w1_i and w2_i.
		for(int j = 0; j < 3; j++){
			mpz_init(derived_key[i][j]);
			mpz_urandomb(derived_key[i][j], rand_Pas_C, pub_moduli_bitsize);
		}
		gmp_randseed(rand02, derived_key[i][0]); // generates a seed for a_i.
		gmp_randseed(rand03, derived_key[i][1]); // generates a seed for sw1.
		gmp_randseed(rand04, derived_key[i][2]); // generates a seed sw2.
		for(int j = 0; j < NoElem_in_bucket + 1; j++){ // generates d+1 pseudorandom coefficients.
			mpz_init(Sw1[j]);
			mpz_init(Sw2[j]);
			mpz_urandomb(Sw1[j], rand03, pub_moduli_bitsize);
			mpz_urandomb(Sw2[j], rand04, pub_moduli_bitsize);
		}
		Polynomial pol_1, pol_2;
		w_A_[i] = pol_1.evaluate_coeffs(Sw1, xpoints, NoElem_in_bucket + 1, xpoint_size, pu_moduli[0]);
		w_B_[i] = pol_2.evaluate_coeffs(Sw2, xpoints, NoElem_in_bucket + 1, xpoint_size, pu_moduli[0]);
		for(int j = 0; j < xpoint_size; j++){ // generates v_A, v_B and a_i in step d.7 in the protocol.
				mpz_init(a_[i][j]);
				mpz_init(tmp_bl_[i][j]);
				mpz_urandomb(a_[i][j], rand02, pub_moduli_bitsize);
				mpz_urandomb(tmp_bl_[i][j], rand_Pas_D, pub_moduli_bitsize);
			}
			mpz_clear(derived_key[i][0]);
			mpz_clear(derived_key[i][1]);
			mpz_clear(derived_key[i][2]);
			free(derived_key[i]);
		}
		gmp_randclear(rand02);
		gmp_randclear(rand03);
		gmp_randclear(rand04);
		gmp_randclear(randC);
		gmp_randclear(randD);
		gmp_randclear(rand_Pas_C);
		gmp_randclear(rand_Pas_D);
		free(Sw1);
		free(Sw2);
		free(derived_key);
		mpz_clear(temp_key);
		mpz_clear(temp_key2);
	}
//**********************************************************************
// - Function description: given a client's permission to compute the PSI, it generates the result.

 Server_Result * Server::compute_result (GrantComp_Info * grantComp_info, bigint tmp_key_){

	 bigint** a, **w_A, **w_B, **tmp_bl;
	 regen_PRpolys(grantComp_info->seed, w_A,w_B, a, tmp_bl, tmp_key_);
	 int indx_A,indx_B;
	 bigint *o_A, *o_B,* *res, buf_1, buf_2, *temp_BFs, tmp_bf_A, tmp_bf_B;
	 mpz_init(buf_1);
	 mpz_init(buf_2);
	 res = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	 o_A = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	 o_B = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	 temp_BFs = (mpz_t*)malloc(table_size * sizeof(mpz_t));
	 for(int i = 0; i < table_size; i++){
		 o_A = get_client_bin(grantComp_info->pm[i][0], grantComp_info->id[1], tmp_bf_A, indx_A);
		 o_B = get_client_bin(grantComp_info->pm[i][1], grantComp_info->id[0], tmp_bf_B, indx_B);
		 mpz_init_set(temp_BFs[indx_B], tmp_bf_B);
		 mpz_clear(tmp_bf_B);
		 mpz_clear(tmp_bf_A);
		 res[indx_B] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		 // given the map,  pseudorandom values, and clients' datasets, it computes the final result.
		 for(int j = 0; j < xpoint_size; j++){
			 mpz_init(res[indx_B][j]);
			 mpz_mul(buf_1, w_A[indx_A][j], o_A[j]); // multiplies client A's dataset by the pseudorandom polynomial: w_A.
			 mpz_add(buf_2, o_B[j], tmp_bl[indx_B][j]); // adds client B's dataset by the blinding factor it sent to the server.
			 mpz_mul(buf_2, w_B[indx_B][j], buf_2);  // multiplies the value in generated one line above by the pseudorandom polynomial: w_B.
			 mpz_add(buf_2, buf_2,buf_1);
			 mpz_add(buf_2, buf_2, a[indx_A][j]); // adds all the product computed above together and with the blinding facotr a[indx_A][j].
			 mpz_mod(res[indx_B][j], buf_2, pu_moduli[0]);
			 mpz_clear(o_A[j]);
			 mpz_clear(o_B[j]);
			 mpz_clear(a[indx_A][j]);
			 mpz_clear(w_B[indx_B][j]);
			 mpz_clear(w_A[indx_A][j]);
		 }
	 }
	 Server_Result* ptr;
	 ptr = new Server_Result;
	 ptr->BF = temp_BFs;
	 ptr->result = res;
	 free(tmp_bl);
	 free(w_A);
	 free(w_B);
	 free(a);
	 free(o_A);
	 free(o_B);
	 mpz_clear(buf_1);
	 mpz_clear(buf_2);
	 return ptr;
 }
//**********************************************************************
// - Function description: given client's ID, it finds the client's dataset index at the server-side.

 int Server::find_db_index(string id){
	 int i;
	 string s;
	 for(i = 0; i < db_size; i++){ // in the server-side database containing  different clients datasets, find the clinet's index in there.
		 Client_Dataset p = get_db(i);
		 if(p.poly[0].get_poly_ID() == id){
			 return i;
			 }
	 }
	 if(i == db_size){
		 cout<<"There is exist no poly. in server with ID:"<<id<<endl;
		 return i = 1000000;// the value only is used to indicate NULL
	 }
}
//**********************************************************************
// - Function description: given client's update query, it applies the update to the client's dataset at the server-side
// It is called by update_client_bin();

void Server::update_db(bigint* vals,bigint label, string id, bigint bbf){
	bool found;
	int index = find_db_index(id); // finds the position in which the client's dataset is stored in the array of clients' datasets.
	for(int i = 0; i < table_size; i++){
		if(mpz_cmp(db[index].labels[i], label) == 0){ // finds the client's bin tagged with the label.
			//apply the changes
			mpz_set(db[index].BF[i], bbf);
			db[index].poly[i].values = vals;
			found = true;
			break;
		}
	}
	if(!found){cout<<"\n Update didn't take place as the label does not exist!"<<endl;}
}
//**********************************************************************
// - Function description: allows the client to send an update query to the server.

void Server::update_client_bin(bigint* vals, bigint label, string id, bigint bbf){

	update_db(vals, label, id, bbf);
}
//**********************************************************************
// - Function description: given a client's label, it returns the corresponding bin tagged with the label.

bigint* Server::get_client_bin(bigint label, string id, bigint& bf, int& indx){

	int j = -1;
	bigint *res;
	res = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	int index = find_db_index(id); // finds the position in which the client's dataset is stored in the array of clients' datasets.
	Client_Dataset db = get_db(index);
	for(int i = 0; i < table_size; i++){ // finds the client's bin tagged with the label.
		if(mpz_cmp(db.labels[i],label) == 0){
			j = i;
			break;
		}
	}
	if(j == -1) {cout<<"\n The label does not exist!"<<endl;return 0;}
	else{
		res = db.poly[j].get_values();
		mpz_init_set(bf, db.BF[j]);
		indx = j;
		return res;
	}
}
//**********************************************************************
