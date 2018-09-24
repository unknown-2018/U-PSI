/*

- Client-side computation of the U-PSI protocol

*/

#include"Client.h"
//**********************************************************************
// - Description:  Constructor- It generates the following private keys: seed,
// counter_key, label_key,shuffle_key, BF_key, BF_counter_key.
// It fetches the hashtable parameters, x-points, and field description from the server.
// It sets a moduli for blinding bloom filters and sets the bloom filter parameters.


Client::Client(Server*server, bigint *elements, int el_size){

	Random rd_;
	labels_bit_size = 64;
	elem = (mpz_t*)malloc(el_size * sizeof(mpz_t));
	elem = elements;
	serv = server;
	elem_size = el_size;
	//keys are generated.
	Random rd;
	gmp_randstate_t rand, rand2, rand3, rand4, rand5, rand6;
	bigint ran;
	rd.init_rand3(rand, ran, 8); // generates a truly random value.
	mpz_init_set(seed, ran); // assing the above random value to seed.
	mpz_init(ran);
	rd.init_rand3(rand2, ran, 8);
	mpz_init_set(counter_key, ran);
	mpz_init(ran);
	rd.init_rand3(rand3, ran, 8);
	mpz_init_set(label_key, ran);
	mpz_init(ran);
	rd.init_rand3(rand4, ran, 8);
	mpz_init_set(shuffle_key, ran);
	mpz_init(ran);
	rd.init_rand3(rand5, ran, 8);
	mpz_init_set(BF_key, ran);
	mpz_init(ran);
	rd.init_rand3(rand6, ran, 8);
	mpz_init_set(BF_counter_key, ran);
	int size;
	get_xpoints(size); //fetches x-coordonates from the server.
	xpoint_size = size;
	get_pubModuli(); //fetches the public moduli (description of the field) from the server.
	get_NoElem_in_bucket(); // gets the bins' capacity from the server.
	get_tablesize(); // gets the hash table length from the server.
	get_pubModuli_bitsize(); // gets the bit-size of the public moduli.
	// sets the bloom filters parameters.
	bf_parameters.projected_element_count = NoElem_in_bucket;
	bf_parameters.false_positive_probability = 0.0000000000009095; // 2^{-40}
	bf_parameters.random_seed = 0xA5A5A5A5;
	if (!bf_parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
	}
	bf_parameters.compute_optimal_parameters();
	pr_moduli_bitsize = 5780; // bit size of a bloom filter when it holds upto 100 elements, and the false positive rante is 2^{-40}.
	//generates a prime number: pr_moduli, of the above size.
	bigint* pr_mod;
	pr_mod = (mpz_t*)malloc(1 * sizeof(mpz_t));
	pr_mod = rd_.gen_randSet(1, pr_moduli_bitsize);
	mpz_nextprime(pr_mod[0], pr_mod[0]);
	mpz_init_set(pr_moduli, pr_mod[0]);
	//sets the counter.
	counter = new int [table_size];
	for(int i = 0; i < table_size; i++){
		counter[i] = 0;
	}
}
//**********************************************************************
// - Function description: Regenerates the most recent blinding factors used to blind the client's y-coordinates.

bigint** Client::regen_bl_factors(bigint seed_, bigint ck, int* counter){

	bigint** blf, bliding_key, bin_key, bin_ck ;
	blf = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	gmp_randstate_t rand,rand2,rand3,rand4;
	gmp_randinit_default(rand);
	gmp_randseed(rand,seed_);
	gmp_randinit_default(rand2);
	gmp_randseed(rand2,ck);
	gmp_randinit_default(rand3);
	gmp_randinit_default(rand4);
	mpz_init(bin_key);
	mpz_init(bin_ck);
	mpz_init(bliding_key);
	//re-generates blinding key: k_j and counter key: ck_j for each j-th bins:
	for(int i = 0;i < table_size; i++){
		blf[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_urandomb(bin_key, rand, pub_moduli_bitsize); //k_j is generated
		mpz_urandomb(bin_ck, rand2, pub_moduli_bitsize); //ck_j is generated
		// If the bin has not been updated, use only k_j.
		if(counter[i] == 0){
			mpz_set(bliding_key, bin_key);
		}
		// If the bin has been updated, then derive a pseudorandm value: c_j and
		// then compute a blinding key as: k_j+PRF(ck_j,c_j)
		else{
		bigint tmp_key;
		gmp_randseed(rand3, bin_ck);
		//compute PRF(ck_j,c_j) in step 2.a
		for(int k = 0; k < counter[i]; k++){
			mpz_init(tmp_key);
			mpz_urandomb(tmp_key, rand3, pub_moduli_bitsize);
		}
		bigint tmp2;
		mpz_init(tmp2);
		mpz_add(tmp2,bin_key, tmp_key);
		mpz_set(bliding_key, tmp2);
		mpz_clear(tmp_key);
		mpz_clear(tmp2);
	}
	gmp_randseed(rand4, bliding_key);
	// given the above blinding key, regenerate the blinding factors of the bin.
	for(int j = 0; j < xpoint_size; j++){
		mpz_init(blf[i][j]);
		mpz_urandomb(blf[i][j], rand4, pub_moduli_bitsize);
	}}
	mpz_clear(bliding_key);
	mpz_clear(bin_ck);
	gmp_randclear(rand);
	gmp_randclear(rand2);
	gmp_randclear(rand3);
	gmp_randclear(rand4);
	return blf;
}
//**********************************************************************
// - Function description: blindes a shuffled blinding factors.

bigint** Client::blind_shuffled_bl(bigint** s_bl, int table_size, bigint seed_, bigint pubmoduli_){

	bigint** blf, bin_key;
	blf= (mpz_t**)malloc(table_size * sizeof(mpz_t));
	gmp_randstate_t rand,rand2;
	gmp_randinit_default(rand);
	gmp_randinit_default(rand2);
	gmp_randseed(rand, seed_);
	mpz_init(bin_key);
	//derives a key for each bin from seed_.
	for(int i = 0; i < table_size; i++){
		blf[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_urandomb(bin_key, rand, pub_moduli_bitsize);//k_j is generated
		gmp_randseed(rand2, bin_key);
		//using the derived key, it generates "xpoint_size" pseudorandom values and uses them to blind the elemenents of s_bl.
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(blf[i][j]);
			mpz_urandomb(blf[i][j], rand2, pub_moduli_bitsize);
			mpz_add(blf[i][j], blf[i][j], s_bl[i][j]);
			mpz_mod(blf[i][j], blf[i][j], pubmoduli_);
		}
	}
	mpz_clear(bin_key);
	gmp_randclear(rand);
	gmp_randclear(rand2);
	return blf;
}
//**********************************************************************
// - Function description: updates client's outsourced data, i.e. deletes/inserts an element from/to the data.

string Client::update(bigint elem, string updt, bigint & label, string id){

	double start_0 = clock();
	string status;
	bigint *labels, *bin,zero, minus_one, bf, *unblinded_bf, *temp_elem, *temp_res,zz, *blinded_bf, *temp_ar;
	temp_ar = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	mpz_init_set_str(zero, "0", 10);
	mpz_init_set_str(minus_one, "-1", 10);
	int number_of_roots = 0;
	int num_of_elements_found = 0;
	bin = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	temp_elem = (mpz_t*)malloc(1 * sizeof(mpz_t));
	temp_res = (mpz_t*)malloc(1 * sizeof(mpz_t));
	mpz_init_set(temp_elem[0], elem);
	int tmpx;
	int j = gen_binIndx(elem, table_size); // generates the bin index to which the element: elem, belongs to.
	double end_0 = clock();
	double temp = end_0 - start_0;
	labels = gen_labels(table_size,label_key); // generates the label for that bin.
	mpz_init_set(label,labels[j]);
	bin = serv->get_client_bin(label, id, bf, tmpx); // retrieves bin: o_j, from the server.
	unblinded_bf = unblind_BF(bf, j, BF_key, BF_counter_key, 5780, pr_moduli); // unblinds the bin's bloom filters.
	double start_2 = clock();
	// re-generates keys and blinding factors.
	bigint key;
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	gmp_randstate_t rand2;
	gmp_randinit_default(rand2);
	gmp_randseed(rand2, counter_key);
	bigint bliding_key;
	bigint counter_key;
	double end_2 = clock();
	temp += end_2 - start_2;
	for(int i = 0 ; i <= j; i++){ //key derivation recall that j is bins' index. Because we need index number j, not j-th index, that's why i<=j.
		mpz_init(key);
		mpz_urandomb(key, rand, pub_moduli_bitsize); //k_j is generated.
		mpz_init(counter_key);
		mpz_urandomb(counter_key, rand2, pub_moduli_bitsize); //ck_j is generated.
	}
	double start_3 = clock();
	bigint* bld_factors;
	bld_factors = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	//sets the most recent blinding key
	if(counter[j] == 0){
		mpz_init_set(bliding_key, key);
	}
	else{
		mpz_init(bliding_key);
		bigint tmp_key;
		gmp_randstate_t rand3;
		gmp_randinit_default(rand3);
		gmp_randseed(rand3, counter_key);
		//computes PRF(ck_j,c_j) in step 2.a of the protocol.
		for(int k = 0; k < counter[j]; k++){
			mpz_init(tmp_key);
			mpz_urandomb(tmp_key, rand3, pub_moduli_bitsize);
		}
		bigint tmp2;
		mpz_init(tmp2);
		mpz_add(tmp2, key, tmp_key);
		mpz_set(bliding_key, tmp2);
	}
	bigint *un_bl;
	un_bl = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	// removes the blinding factors.
	gmp_randstate_t rand1;
	gmp_randinit_default(rand1);
	gmp_randseed(rand1, bliding_key);
	for(int k = 0; k < xpoint_size; k++){
		mpz_init(bld_factors[k]);
		mpz_init(un_bl[k]);
		mpz_urandomb(bld_factors[k], rand1, pub_moduli_bitsize);
		mpz_sub(un_bl[k], pubmoduli, bld_factors[k]);
		mpz_add(un_bl[k], un_bl[k], bin[k]);
		mpz_mod(un_bl[k], un_bl[k], pubmoduli);
	}
	counter[j] += 1; //increments the counter[j].
	int temp_counter = 0;
	temp_res = check_vals_in_BF(temp_elem, 1, unblinded_bf[0], bf_parameters, temp_counter); //checks if the elements (to be inserted/deleted) exists in the bin's bloom filters.
	bigint* new_Bigint_BF;
	new_Bigint_BF = unblinded_bf;
	if(updt == "insertion"){ // If the update is the insertion of the element: elem.
		number_of_roots = 0;
		num_of_elements_found = 0;
		if (mpz_cmp(temp_res[0],elem) != 0){
			Polynomial pol_;
			// exteracts the set elemenets from the bin by interpolating a polynomial and finding its roots.
			bigint* coeff = pol_.interpolate(xpoint_size, xpoints, un_bl, pubmoduli); // interpolate a polynomial, given x and y coordinates.
			bigint* roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli); //finds the roots.
			bigint* valid_roots = check_vals_in_BF(roots, number_of_roots, unblinded_bf[0], bf_parameters, num_of_elements_found); // extracts the valid roots, by checking which roots are in the bin's bloom filters.
			bloom_filter filter(bf_parameters); // builds a new bloom filters for the bin.
			string s = mpz_get_str(NULL, 10, elem);
			filter.insert(s); // inserts the element, to be inserted, into the new bloom filters.
			// inserts the set elements of the bin into the new bloom flters.
			if(num_of_elements_found > 0){
				for(int i = 0; i < num_of_elements_found; i++){
					string s = mpz_get_str(NULL, 10, valid_roots[i]);
					filter.insert(s);
				}
			}
			new_Bigint_BF = convert_BF_to_bigint(filter); // converts the filter to a bigint.
			// inserts the valid roots into the bin and pads it with -1 if the number of valid roots are smaller than the bin's capacity.
			// later on (when a polynomial is constructed) the -1's are replaced with random values.
			bigint* padded_bin;
			padded_bin = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
			mpz_init_set(padded_bin[0],elem);
			for(int i = 1; i < NoElem_in_bucket; i++){
				if(i < num_of_elements_found + 1){
					mpz_init_set(padded_bin[i], valid_roots[i - 1]);}
					else{mpz_init_set(padded_bin[i], minus_one);}
			}
			Polynomial pol(padded_bin, outpoly_ID, xpoints, NoElem_in_bucket, xpoint_size, pubmoduli); // builds a Polynomial. Here, -1's are replaced with random values.
			un_bl = pol.values; // sets the y-coordinates.
			string ss = mpz_get_str(NULL, 10, elem);
			status="\nElement "+ss+" has been inserted";
		}
		else{cout<<"\nElement "<<elem<<" already exists in the set, no insertion was required"<<endl;}
	}
	if(updt == "deletion"){ // If the update is the deletion of the element: elem.
		number_of_roots = 0;
		num_of_elements_found = 0;
		if (mpz_cmp(temp_res[0], elem) == 0){
			Polynomial pol_;
			// exteracts the set elemenets from the bin by interpolating a polynomial and finding its roots.
			bigint* coeff = pol_.interpolate(xpoint_size, xpoints, un_bl, pubmoduli);
			bigint* roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli);//find the roots
			bigint* valid_roots = check_vals_in_BF(roots, number_of_roots, unblinded_bf[0], bf_parameters, num_of_elements_found);
			bloom_filter filter(bf_parameters); // builds a new bloom filters.
			if(num_of_elements_found > 0){
				for(int i = 0; i < num_of_elements_found; i++){
					// inserts all set eleements of the bin into the bloom filters, except the element to be deleted: elem.
					if(mpz_cmp(valid_roots[i], elem) != 0){
						string s = mpz_get_str(NULL, 10, valid_roots[i]);
						filter.insert(s);
					}
				}
			}
			new_Bigint_BF = convert_BF_to_bigint(filter);
			// inserts the valid roots (except elem) into the bin and pads it with -1's if needed.
			bigint* padded_bin;
			padded_bin = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
			mpz_init_set(padded_bin[0], elem);
			for(int i = 0; i < NoElem_in_bucket; i++){
				if(i < num_of_elements_found){
					if(mpz_cmp(valid_roots[i], elem) != 0){
						mpz_init_set(padded_bin[i], valid_roots[i]);}
					else{mpz_init_set(padded_bin[i], minus_one);}
				}
					else{mpz_init_set(padded_bin[i], minus_one);}
			}
			Polynomial pol(padded_bin, outpoly_ID, xpoints,NoElem_in_bucket, xpoint_size,pubmoduli); // builds a Polynomial
			un_bl = pol.values; // sets y-coordinates
			string ss = mpz_get_str(NULL, 10, elem);
			status="Element "+ss+" has been deleted";
		}
		else{cout<<"Element "<<elem<<" does not exist in the set, no deletion was required"<<endl;}
	}
	double end_3 = clock();
	temp += end_3 - start_3;
	blinded_bf = blind_BF(new_Bigint_BF[0], j, BF_key, BF_counter_key, 5780, pr_moduli); // blinds the biginteger representing the blom filters.
	// blind the y-coordinates of the bin with fresh blinding factors.
	double start_4 = clock();
	bigint tmp_key2;
	gmp_randstate_t rand_x, rand_y;
	gmp_randinit_default(rand_x);
	gmp_randinit_default(rand_y);
	gmp_randseed(rand_x, counter_key);
	//compute PRF(ck_j,c_j) in step 2.a
	for(int k = 0; k < counter[j]; k++){
		mpz_init(tmp_key2);
		mpz_urandomb(tmp_key2, rand_x, pub_moduli_bitsize);
	}
	bigint tmp2;
	mpz_init(tmp2);
	mpz_add(tmp2, key, tmp_key2);
	gmp_randseed(rand_y, tmp2);
	bigint *bl2;
	bl2 = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	for(int k = 0; k < xpoint_size; k++){
		mpz_init(bld_factors[k]);
		mpz_init(bl2[k]);
		mpz_urandomb(bld_factors[k], rand_y, pub_moduli_bitsize);
		mpz_add(bl2[k], bld_factors[k], un_bl[k]);
		mpz_mod(bl2[k], bl2[k], pubmoduli);
	}
	double end_4 = clock();
	temp += end_4 - start_4;
	float diff = temp / (double) CLOCKS_PER_SEC;
	cout<<"\nIn Update-excluding: (1)blind_BF (2) Unblind_BF and (3)key_gen loop (4)gen_label:"<<diff<<endl;
	serv->update_client_bin(bl2, label, outpoly_ID, blinded_bf[0]); // asks the server to update the client's dataset.
	return status;
}
//**********************************************************************
// - Function description: generates a set of distinct pseudorandom labels.

bigint* Client::gen_labels(int size, bigint seed){
	
	double start_1 = clock();
	double end_1;
	mpz_t *labels;
	labels=(mpz_t*)malloc(size * sizeof(mpz_t));
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand,seed);
	for(int i = 0; i < size; i++){
		mpz_init(labels[i]);
		mpz_urandomb(labels[i], rand, labels_bit_size);// The last argument is in bit
		if(i == 0){end_1 = clock();} //just to time gen label in update
	}
	// checks if all the labels are distinct.
	for(int i = 0; i < size; i++){
		for(int j = i + 1 ;j < size; j++){
				if(mpz_cmp(labels[i], labels[j]) == 0) {
					cout<<"\nPick another Label key"<<endl;
					break;
				}
		}
	}
	gmp_randclear(rand);
	double temp = end_1 - start_1;
	float f_temp = temp / (double) CLOCKS_PER_SEC;
	cout<<"\n In gen_label: time to gen a label:"<<f_temp<<endl;
	return labels;
}
//**********************************************************************
// - Function description: generates a permutation map that includes a set of pairs (l_i,l_j)
// where each pair includes pseudorandom labels each blong to a different client.

bigint** Client::gen_map(int size, bigint seed1, bigint seed2){

	bigint **labels, **labels2;
	labels= (mpz_t**)malloc(size * sizeof(mpz_t));
	gmp_randstate_t rand1, rand2;
	gmp_randinit_default(rand1);
	gmp_randinit_default(rand2);
	gmp_randseed(rand1, seed1);
	gmp_randseed(rand2, seed2);
	for(int i = 0; i < size; i++){
		labels[i] = (mpz_t*)malloc(2 * sizeof(mpz_t));
		mpz_init(labels[i][0]);
		mpz_init(labels[i][1]);
		mpz_urandomb(labels[i][0], rand1, labels_bit_size);// The last argument is in bit
		mpz_urandomb(labels[i][1], rand2, labels_bit_size);
	}
	labels2 = R_shuffle(labels, size); // the label pairs is randomly permuted.
	gmp_randclear(rand1);
	gmp_randclear(rand2);
	free(labels);
	return labels2;
}
//**********************************************************************
// - Function description: given two arrays permuted (under two different keys), it finds their matching indices.
// It is called by find_matched_bins().
 
int* Client::find_matches(int* a, int* b, int size){
	 
	 int* res;
	 res=new int[size];
	 for(int i = 0; i < size; i++){
		 for(int j = 0; j < size; j++){
			 if(a[i] == b[j]){
				 res[i] = j;
				 break;
			 }
		 }
	 }
	 return res;
 }
 //**********************************************************************
// - Function description: finds the match between two arrays permuted under two different keys.
// In particular, for an index "i" in vecotr a permuted using key k_1, it finds index "j" in anoher
// vector permuted under key k_2 such that "i" and "j" belong to the same index before they were
// permuted. It allocates values from 1 to "size" to each vector val_1 and va_2. Then,
// it permutes them and then uses find() to find the matched indices in the permuted vectors.

int* Client::find_matched_bins(bigint k_1, bigint k_2, int size){

	int* val_1, *val_2, *res;
	val_1 = new int[size];
	val_2 = new int[size];
	res = new int[size];
	int *permuted_1, *permuted_2;
	for( int j = 0; j < size; j++){
		val_1[j] = j + 1;
		val_2[j] = j + 1;
	}
	permuted_1 = PR_shuffle(val_1, size, k_1);
	permuted_2 = PR_shuffle(val_2, size, k_2);
	res = find_matches(permuted_1, permuted_2, size);
	free(permuted_1);
	free(permuted_2);
	return res;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of integers, using Fisher-Yates shuffle algorithm.

int* Client::PR_shuffle(int* elem, int size, bigint seed){
	
	int *result;
	result = new int [size];
	int buffer;
	bigint r, big_j;
	mpz_init(r);
	mpz_init(big_j);
	for(int i = 0; i < size; i++){
		result[i] = elem[i];
	}
	int indx = 0;
	// use the seed to generate a random value between [1,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand,seed);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(r, rand, pub_moduli_bitsize);// Here, pub_moduli_bitsize is an arbitrary choice. It would fine as long as it's greater than the ceiling.
    // gen random value in in range [0,j]
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buffer = result[j];
		result[j] = result[indx];
		result[indx] = buffer;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return result;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of big-integers, using Fisher-Yates shuffle algorithm.

bigint*  Client::PR_shuffle(bigint* elem, int size, bigint seed){

	bigint *res,buf, r, big_j;
	mpz_init(r);
	mpz_init(buf);
	mpz_init(big_j);
	res=(mpz_t*)malloc(size * sizeof(mpz_t));
	for(int i = 0; i < size; i++){
		mpz_init(res[i]);
		mpz_set(res[i], elem[i]);
	}
	int indx = 0;
	// use the seed to generate a random value between [1,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(r, rand, pub_moduli_bitsize);// Here, pub_moduli_bitsize is an arbitrary choice. It is fine as long as it's greater than the ceiling.
    // gen random value in in range [0,j]
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		mpz_init_set(buf, res[j]);
		mpz_init(res[j]);
		mpz_set(res[j], res[indx]);
		mpz_set(res[indx], buf);
	}
	gmp_randclear(rand);
	mpz_clear(buf);
	mpz_clear(r);
	return res;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of bins, using Fisher-Yates shuffle algorithm.

bigint**  Client::PR_shuffle_bins(bigint** bins, int size, bigint seed_){

	bigint **s_bins, *buf,big_j, r;//blinding factors
	s_bins = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	buf = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	for(int i = 0; i < size; i++){
		s_bins[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(s_bins[i][j]);
			mpz_set(s_bins[i][j], bins[i][j]);
	}}
	mpz_init(r);
	mpz_init(big_j);
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed_);
	int indx = 0;
	// use the seed to generate a random value between [0,j]
	for (int  j = size - 1; j > 0; j--){
		// gen random value in in range [0,j]
		mpz_urandomb(r, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buf = s_bins[j];
		s_bins[j] = s_bins[indx];
		s_bins[indx] = buf;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return s_bins;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of polynomials, using Fisher-Yates shuffle algorithm.

Polynomial*  Client::PR_shuffle_poly(Polynomial* pol, int size, bigint seed){

	Polynomial *ply_res, buf;
	ply_res = new Polynomial[size];
	for(int i = 0; i < size; i++){
		ply_res[i] = pol[i];
	}
	int indx = 0;
	// use the seed to generate a random value between [0,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	bigint r, big_j;
	mpz_init(r);
	mpz_init(big_j);
	for (int j = size - 1; j > 0; j--){
		// gen random value in in range [0,j]
		mpz_urandomb(r, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buf = ply_res[j];
		ply_res[j] = ply_res[indx];
		ply_res[indx] = buf;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return ply_res;
}
//**********************************************************************
// - Function description: "randomly" permutes an array of label pairs, using Fisher-Yates shuffle algorithm.

bigint**  Client::R_shuffle(bigint** elem, int size){

	bigint **res;
	res=(mpz_t**)malloc(size * sizeof(mpz_t));
	for (int i = 0; i < size; i++){
		res[i]=(mpz_t*)malloc(2 * sizeof(mpz_t));
	}
	for(int i = 0; i < size; i++){
		mpz_init(res[i][0]);
		mpz_init(res[i][1]);
		mpz_set(res[i][0], elem[i][0]);
		mpz_set(res[i][1], elem[i][1]);
	}
	int indx = 0;
	bigint buf, buf1;
	Random rd;
	gmp_randstate_t rand;
	bigint ran;
	rd.init_rand3(rand, ran, 8);
	bigint temp, big_j;
	mpz_init(temp);
	mpz_init(big_j);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(temp, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(temp, temp, big_j);
		indx=mpz_get_ui(temp);
		mpz_init_set(buf, res[j][0]);
		mpz_init_set(buf1, res[j][1]);
		mpz_init(res[j][0]);
		mpz_init(res[j][1]);
		mpz_set(res[j][0], res[indx][0]);
		mpz_set(res[j][1], res[indx][1]);
		mpz_set(res[indx][0], buf);
		mpz_set(res[indx][1], buf1);
	}
	gmp_randclear(rand);
	mpz_clear(temp);
	mpz_clear(buf);
	mpz_clear(buf1);
	return res;
}
//**********************************************************************
// - Function description: given two arrays: v_a and v_b, permuted under two distinct keys: pk_1 and pk_2, it allows one to find a mateched element pairs: (e_{a,i},e_{b.j})
// such that both elements of each pair belong to the same index before the permutation. Then, it adds up the elements in each pair.

bigint** Client::combine_permuted_bins(bigint**& v_a, bigint**& v_b, bigint**& a, int v_size, int xpoint_size_, bigint pk_1, bigint pk_2, bigint pubmoduli){

	bigint** res;
	int* ar;
	res = (mpz_t**)malloc(v_size * sizeof(mpz_t));
	ar = new int[v_size];
	ar = find_matched_bins( pk_1, pk_2, v_size);
	for(int i = 0; i < v_size; i++){
		res[ar[i]] = (mpz_t*)malloc(xpoint_size_ * sizeof(mpz_t));
		// sums the elements of v_a at position i with its matched elements in v_b (at position ar[i]).
		// then sums the result with the elements of a at position i.
		// stores the result in res at position ar[i],  which is compatible with the permutation used in v_b.
		for(int j = 0; j < xpoint_size_; j++){
			mpz_init(res[ar[i]][j]);
			mpz_add(res[ar[i]][j], v_a[i][j], v_b[ar[i]][j]);
			mpz_add(res[ar[i]][j], res[ar[i]][j], a[i][j]);
			mpz_mod(res[ar[i]][j], res[ar[i]][j], pubmoduli);
			mpz_clear(v_a[i][j]);
			mpz_clear(v_b[ar[i]][j]);
			mpz_clear(a[i][j]);
		}
	}
	free(ar);
	return res;
}
//**********************************************************************
// - Function description: converts a bloom filter to a biginteger value.

bigint* Client::convert_BF_to_bigint(bloom_filter filter){
	
	int size = filter.bit_table_.size();
	bigint* res;
	res = (mpz_t*)malloc(1 * sizeof(mpz_t));
	unsigned char ar[size];
	for(int i = 0; i < size; i++){
		ar[i] = filter.bit_table_[i];
	}
	mpz_init(res[0]);
	mpz_import(res[0], sizeof(ar), 1, sizeof(ar[0]), 0, 0, ar);// converts the array of bytes to a biginteger.
	return res;
}
//**********************************************************************
// - Function description: converts a biginteger representing a bloom filter to a bloom filter.
// In the case where BF's very first bis are zero, when it is converted to a biginteger they would be lost.
// So, the function ensures they are put back (using a pad), otherwise BF would not reconstructed correctly.

 bloom_filter Client::convert_bigint_to_BF(bigint a, bloom_parameters parameters){
   // converts bigint to a bitstring and pad it if needed.
	 bloom_filter filter(bf_parameters);
	 filter.clear();
	 int size = filter.bit_table_.size();
	 int offset = 0;
	 int inx = 0;
	 string s_val, pad;
	 s_val=mpz_get_str(NULL, 2, a);
	 //padds the string if needed
	 if(s_val.length() < size * 8){
		 int dif = (size * 8) - s_val.length();
		 for (int j = 0;j < dif; j++){pad += "0";}
	 }
	 s_val = pad + s_val;
   // stores each 8-bits of the string in the hex form in each index of a new bf vector.
	 while (offset < s_val.length() / 8 + 1){
	 	string tmp_binr = s_val.substr(offset * 8,8);
		// converts tmp_binr to hex.
		const unsigned g_unMaxBits = 8;
		bitset <g_unMaxBits> bs(tmp_binr);
		unsigned n = bs.to_ulong();
		stringstream ss;
		ss << hex << n;
		string tmp_hex = "0x" + boost::to_upper_copy(ss.str());
		// stores hex in the filter.
	 	std::istringstream str(tmp_hex);
	 	int num;
	 	str >> std::hex >> num;
	 	filter.bit_table_[inx] = num;
	 	inx++;
	 	offset++;
	 }
	 return filter;
}
//**********************************************************************
// - Function description: given an element, it determines its bin's index in the hash table.

int Client::gen_binIndx(bigint elem, int table_size){
  
	bigint b, zz;
	mpz_init(zz);
	mpz_set_ui(zz, table_size);
	string s_val;
	CryptoPP::SHA512 hash2;
 	s_val=mpz_get_str(NULL, 10, elem);
	unsigned int nDataLen = s_val.length();
	byte digest[CryptoPP::SHA512::DIGESTSIZE];
	hash2.CalculateDigest(digest, (byte*)s_val.c_str(), nDataLen); // hashes the element.
	s_val.clear();
	mpz_init(b);
	mpz_import(b, sizeof(digest), 1, sizeof(digest[0]), 0, 0, digest);
	mpz_mod(b, b, zz);
	int j = mpz_get_ui(b); // converts the hash value to an integer.
	return j;
}
//**********************************************************************
// - Function description: given a hashtable containing set elements, it assigns a bloom filter to each bin of the table.

bigint* Client::assing_BFs2HT(Hashtable HT, int NoElem_in_bucket, int table_size, bloom_parameters parameters){
	
	bloom_filter filter(bf_parameters);
	bigint minus_one;
	mpz_init_set_str(minus_one, "-1", 10);
	bigint* temp_ar, *bigint_BF, *temp_bigint;
	temp_ar = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
	bigint_BF = (mpz_t*)malloc(table_size * sizeof(mpz_t));
	temp_bigint = (mpz_t*)malloc(1 * sizeof(mpz_t));
	// retrives the set elements of each bin, and inserts them to bloom filters.
	for(int i = 0; i < table_size; i++){
		temp_ar = HT.get_bucket(i);
		for(int j = 0; j < NoElem_in_bucket; j++){
			if(mpz_cmp(temp_ar[j], minus_one) > 1){
				string s = mpz_get_str(NULL, 10, temp_ar[j]); // this is done to make query and insertion compatible.
				filter.insert(s);
				s.clear();
			}
		}
		temp_bigint = convert_BF_to_bigint(filter); // converts the bloom filter to a biginteger.

		mpz_init_set(bigint_BF[i], temp_bigint[0]); // stores the biginteger in an array.
		filter.clear();
	}
	return bigint_BF;
}
//**********************************************************************
// - Function description: given an array of values and a biginteger representing a Bloom filter, it returns those
// elements of the array that exist in the bloom filter. In particular, given a polynomial's roots and Bloom filter belonging to the same bin,
// it returns the roots that have been encoded in the bloom filter.

bigint* Client::check_vals_in_BF(bigint* vals, int val_size, bigint bf, bloom_parameters parameters, int& counter){

	bigint *res, minus_one, zero;
	mpz_init_set_str(minus_one, "-1", 10);
	mpz_init_set_str(zero, "0", 10);
	if(mpz_cmp(bf,zero) == 0){
		res = (mpz_t*)malloc(1 * sizeof(mpz_t));
		mpz_init_set(res[0], minus_one);
		counter = 0;
		return res;
	}
	int indx = 0;
	res = (mpz_t*)malloc(val_size * sizeof(mpz_t));
	bloom_filter filter(bf_parameters);
	filter = convert_bigint_to_BF(bf, bf_parameters); //convers biginteger (representation of a bloom filters) to a bloom filters.
	for (int j = 0; j < val_size; j++){
		string s = mpz_get_str(NULL, 10, vals[j]);
		if(filter.contains(s)){
			mpz_init_set(res[indx], vals[j]);
			indx++;
		}
	}
	filter.clear();
	counter = indx;
	return res;
	}
//**********************************************************************
// - Function description: fetches bin's capacity: d, as apart of the hashtable parameter, from the server.

void Client::get_NoElem_in_bucket(){
	
	NoElem_in_bucket = serv->get_NoElem_in_bucket();
}
//**********************************************************************
// - Function description: fetches an array of x-coordinates from the server.

void Client::get_xpoints(int& size){
	
	xpoints = serv->get_xpoints(size);
	xpoint_size = size;
}
//**********************************************************************
// - Function description: retrives the public moduli bit-size from the server.

void Client::get_pubModuli_bitsize(){
	
	pub_moduli_bitsize = serv->get_pubModuli_bitsize();
}
//**********************************************************************
// - Function description: fetches the public moduli from the server.

void Client::get_pubModuli(){
	
	bigint *ptr = (mpz_t*)malloc(1 * sizeof(mpz_t));
	ptr = serv->get_pubModuli();
	mpz_init_set(pubmoduli, ptr[0]);
}
//**********************************************************************
// - Function description: fetches the hashtable length: h, as apart of the hashtable parameter, from the server.

void Client::get_tablesize(){
	
	table_size = serv->get_table_size();
}
//**********************************************************************
// - Function description: prepare the set elements and sends a blinded dataset to the server.

void Client::outsource_db(string& poly_ID){
	
	Client_Dataset db;
	bigint minus_one, *blinded_BF, *permuted_BBF, *bigint_BF, tmp_key, key, ck, bliding_key, tmp2;
	mpz_init_set_str(minus_one, "-1", 10);

	Hashtable HT(NoElem_in_bucket, elem, elem_size, table_size); // contructs a hash table and inserts the element into it.
	bigint_BF = assing_BFs2HT(HT, NoElem_in_bucket, table_size, bf_parameters); // assigns a bloom filters to each bin. It returns an array of bigint representing bloom filters.
	blinded_BF = blind_BFs(bigint_BF, table_size , BF_key, 5780, pr_moduli); // blinds the bigintegers representing the bloom filters.
	db.BF = PR_shuffle(blinded_BF, table_size, shuffle_key); // permutes the array of bigintegers and stores the result in Client_Dataset that will be sent to the server.
	//sets parameters to represent each bin by a polynomial
	Polynomial *poly;
	poly = new Polynomial [table_size];
	outpoly_ID = poly_ID;
	gmp_randstate_t rand,rand2, rand3;
	gmp_randinit_default(rand);
	gmp_randinit_default(rand2);
	gmp_randinit_default(rand3);
	gmp_randseed(rand, seed);
	gmp_randseed(rand2, counter_key);
	mpz_init(key);
	mpz_init(ck);
	mpz_init(bliding_key);// bliding_key is the most recent key used to generate blinding factors.
	// for every index in the hash table, it contructs a polynomial (in poly is decided whether dummy values shuold be used).
	for(int i = 0; i < table_size; i++){
		Polynomial pol(HT.get_bucket(i), poly_ID, xpoints, NoElem_in_bucket, xpoint_size, pubmoduli);
		poly[i] = pol;
		// assigns a seed to every index of HT.	Each seed is used to blind corresponding poly.
		mpz_urandomb(key, rand, pub_moduli_bitsize); // k_j is generated
		mpz_urandomb(ck, rand2, pub_moduli_bitsize); // ck_j is generated
		mpz_set(bliding_key, key);
		poly[i].blind_poly(bliding_key, pubmoduli, pub_moduli_bitsize); // blind every poly.
	}
	bigint* labels;
	labels = (mpz_t*)malloc(table_size * sizeof(mpz_t));
 	labels = gen_labels(table_size, label_key);
	db.labels = PR_shuffle(labels, table_size, shuffle_key);
	db.poly = PR_shuffle_poly(poly, table_size, shuffle_key);
	serv->store_poly(db);
}
//**********************************************************************
// - Function description: generates a request for PSI computation. This request is created by the result recipient client and
// sent to the other client for its permission.

CompPerm_Request * Client::gen_compPerm_req(bigint & tmp_key_){

	bigint **bl, **s_bl, **r, tmp_key;
	gmp_randstate_t rand;
	Random rd;
	rd.init_rand3(rand, tmp_key, 8);
	CompPerm_Request* ptr;
	ptr = new CompPerm_Request;
	bl = regen_bl_factors(seed, counter_key, counter); // regenerates the blinding factors.
	s_bl = PR_shuffle_bins(bl, table_size, shuffle_key); // pseudorandomly permutes the blinding factors.
	r = blind_shuffled_bl(s_bl, table_size, tmp_key, pubmoduli); // blinds the shuffled blinding factors.
	mpz_init_set(tmp_key_, tmp_key);
	// sets the values to be sent to the other client.
	mpz_init_set(ptr->shuffle_key_, shuffle_key);
	mpz_init_set(ptr->label_key_, label_key);
	ptr->r = r;
	ptr->id = outpoly_ID;
	return ptr;
}
//**********************************************************************
// - Function description: given a request for PSI computation, it grants (or does not grant) the computation.

GrantComp_Info * Client::grant_comp(CompPerm_Request* com_req,bigint **&qq, bool accept){
	// if the client does not grant the computation, then it returns NULL.
	GrantComp_Info * ptr;
	ptr = new GrantComp_Info;
	if(!accept){
		ptr = NULL;
		return ptr;
	}
	bigint **derived_key, **a, *Sw1, *Sw2, **q, **pm, **bl, **s_bl, **v_A, **v_B, temp_key, tk;
	v_A = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	v_B = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	a = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	pm = gen_map (table_size, label_key, com_req->label_key_); // generates a randomly permuted mapping vector.
	bl = regen_bl_factors(seed, counter_key, counter); // regenerates the blinding factors.
	s_bl = PR_shuffle_bins(bl, table_size, shuffle_key); // permuted the blinding factors.
	gmp_randstate_t randC, rand_Pas_C;
	gmp_randinit_default(randC);
	gmp_randinit_default(rand_Pas_C);
	// gen a_i, w1_i and w2_i
	q = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	derived_key = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	a = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	Sw1 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	Sw2 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	gmp_randstate_t rand02, rand03, rand04;
	gmp_randinit_default(rand02);
	gmp_randinit_default(rand03);
	gmp_randinit_default(rand04);
	gmp_randstate_t rand;
	Random rd;
	rd.init_rand3(rand, tk, 8);// generate a fresh master seed.
	gmp_randseed(randC, tk);
	for(int i = 0; i < table_size; i++){
		v_A[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		v_B[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		a[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		q[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_init(temp_key);
		mpz_urandomb(temp_key, randC, pub_moduli_bitsize); // generates a key: k_i
		gmp_randseed(rand_Pas_C, temp_key);
		derived_key[i] = (mpz_t*)malloc(3 * sizeof(mpz_t));
		// derives 3 keys from k_i. key 1 is used to generate a_i, while key 2 and 3 used to genetare w1_i and w2_i (i.e. pseudorandom polynomials).
		for(int j = 0; j < 3; j++){
			mpz_init(derived_key[i][j]);
			mpz_urandomb(derived_key[i][j], rand_Pas_C, pub_moduli_bitsize);
		}
		gmp_randseed(rand02, derived_key[i][0]); // generates a seed for a_i.
		gmp_randseed(rand03, derived_key[i][1]); // generates a seed Sw1.
		gmp_randseed(rand04, derived_key[i][2]); // generates a seed Sw2.
		for(int j = 0; j < NoElem_in_bucket + 1; j++){ // generates two sets: Sw1 and Sw2, each of which contains d+1 random coefficients.
			mpz_init(Sw1[j]);
			mpz_init(Sw2[j]);
			mpz_urandomb(Sw1[j], rand03, pub_moduli_bitsize);
			mpz_urandomb(Sw2[j], rand04, pub_moduli_bitsize);
		}
		Polynomial pol_1, pol_2;
		bigint* temp_w1 = pol_1.evaluate_coeffs(Sw1, xpoints, NoElem_in_bucket + 1, xpoint_size, pubmoduli); // evaluates the random coefficients in Sw1 at x-coordinates.
		bigint* temp_w2 = pol_2.evaluate_coeffs(Sw2, xpoints, NoElem_in_bucket + 1, xpoint_size, pubmoduli); // evaluates the random coefficients in Sw2 at x-coordinates.
		// generates v_A, v_B and a_i in step d.7 in the protocol.
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(v_A[i][j]);
			mpz_init(v_B[i][j]);
			mpz_init(a[i][j]);
			mpz_mul(v_A[i][j], s_bl[i][j], temp_w1[j]);
			mpz_mod(v_A[i][j], v_A[i][j], pubmoduli);
			mpz_mul(v_B[i][j], com_req->r[i][j], temp_w2[j]);
			mpz_mod(v_B[i][j], v_B[i][j], pubmoduli);
			mpz_urandomb(a[i][j], rand02, pub_moduli_bitsize);
			}
	}
	q = combine_permuted_bins(v_A, v_B, a, table_size, xpoint_size, shuffle_key, com_req->shuffle_key_, pubmoduli); // computes q in step d.8 in the protocol.
	// sets the values to be sent to the server.
	qq = q;
	mpz_init_set(ptr->seed, tk);
	ptr->pm = pm;
	ptr->id = new string[2];
	ptr->id[0] = com_req->id;// Client B's ID
	ptr->id[1] = outpoly_ID;// Client A's ID
	return ptr;
}
//**********************************************************************
// - Function description: given an array of polynomial's coefficients, it finds and returns the polynomials roots.

bigint* Client::findroots(bigint *coeff, int coeff_size, int& number_of_roots, bigint pubmoduli){

	int counter_roots = 0;
	bigint *res;
	res=(mpz_t*)malloc(coeff_size * sizeof(mpz_t));
	char * tmp_mod = mpz_get_str(NULL, 10, pubmoduli);
	ZZ p = to_ZZ(tmp_mod);
	ZZ_p::init(p);
	ZZ_pX P;
	ZZ one(1);
	for(int j = 0; j < coeff_size; j++){
		char * tmp = mpz_get_str(NULL, 10, coeff[j]);
		ZZ_p dd = to_ZZ_p(conv<ZZ> (tmp));
		SetCoeff(P, j, dd);
	}
	ZZ_p a = LeadCoeff(P);
	ZZ aa = rep(a);
	if(aa > one){
		MakeMonic(P);
	}
	Vec< Pair < ZZ_pX, long > > factors;
	CanZass(factors, P);
	vec_ZZ_p root;
	for(int j = 0; j < factors.length(); j++){
		if(factors[j].a.rep.length() == 2){
			root = FindRoots(factors[j].a);
			for(int k = 0; k < root.length(); k++){
				stringstream ss;
				ss << root[k];
				string tmpm = ss.str();
				char cv[tmpm.length()];
				strcpy(cv,tmpm.c_str());
				mpz_init_set_str(res[counter_roots], cv, 10);
				counter_roots++;
			}
		}
	}
	number_of_roots = counter_roots;
	return res;
}
//**********************************************************************
// - Function description: given an array of blinded bigintegers, representing bloom filters, it unblinds them.

bigint* Client::unblind_BFs(bigint* BF, int table_size_, bigint BF_key, bigint BF_counter_key, int bit_size, bigint pr_moduli_){

	bigint *unblinded_BF, *bld_factor, *shuffled_bld;
	bld_factor = (mpz_t*)malloc(table_size_ * sizeof(mpz_t));
	gmp_randstate_t rand_y, rand_x;
	gmp_randinit_default(rand_y);
	gmp_randinit_default(rand_x);
	gmp_randseed(rand_y, BF_key);
	gmp_randseed(rand_x, BF_counter_key);
	bigint bf_key, bf_ck;
	for(int i = 0;i < table_size_; i++){
		mpz_init(bf_key);
		mpz_urandomb(bf_key, rand_y, bit_size);
		mpz_init(bf_ck);
		mpz_urandomb(bf_ck, rand_x, bit_size);
		if(counter[i] == 0){
			mpz_init_set(bld_factor[i], bf_key);
		}
		else{
			bigint tmp_key;
			gmp_randstate_t rand3;
			gmp_randinit_default(rand3);
			gmp_randseed(rand3, bf_ck);
			for(int k = 0;k < counter[i]; k++){
				mpz_init(tmp_key);
				mpz_urandomb(tmp_key, rand3, bit_size);
			}
			mpz_init(bld_factor[i]);
			mpz_add(bld_factor[i], tmp_key, bf_key);
			mpz_mod(bld_factor[i], bld_factor[i], pr_moduli_);
		}
	}
	shuffled_bld = PR_shuffle(bld_factor, table_size_, shuffle_key); //shuffles the blinding factors, as the bloom filters are in the shuffled form.
	unblinded_BF = (mpz_t*)malloc(table_size_ * sizeof(mpz_t));
	for(int j = 0; j < table_size_; j++){
		mpz_init(unblinded_BF[j]);
		mpz_sub(unblinded_BF[j], pr_moduli_, shuffled_bld[j]);
		mpz_add(unblinded_BF[j], unblinded_BF[j], BF[j]);
		mpz_mod(unblinded_BF[j], unblinded_BF[j], pr_moduli_);
	}
	return unblinded_BF;
}
//**********************************************************************
// - Function description: given the server's response, it find the intersection. In particular,
// it  unblinds the result, finds the polynomials roots and retrives the valid ones.

vector <string> Client::find_intersection(Server_Result* res, int*& size, bigint** q){

	vector <string> all_valid_roots;
	bigint *un_bl, *unbl_BFs;
	bigint* valid_roots, *roots, *coeff;
	un_bl = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	coeff = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	roots = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	unbl_BFs = unblind_BFs(res->BF, table_size, BF_key, BF_counter_key, 5780, pr_moduli);
	string tempstr;
	for(int i = 0; i < table_size; i++){
		// removes the blinding factors
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(un_bl[j]);
			mpz_sub(un_bl[j], pubmoduli, q[i][j]);
			mpz_add(un_bl[j], un_bl[j], res->result[i][j]);
			mpz_mod(un_bl[j], un_bl[j], pubmoduli);
		}
		int number_of_roots = 0;
		int num_of_elements_found = 0;
		coeff = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		roots = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		Polynomial pol;
		coeff = pol.interpolate(xpoint_size, xpoints, un_bl,pubmoduli); // interpolates a polynomial given x and y coordinates.
		roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli); // finds the roots of the interpolated polynomial.
		free(coeff);
		valid_roots = check_vals_in_BF(roots, number_of_roots, unbl_BFs[i], bf_parameters, num_of_elements_found); // extracts the valid roots.
		for(int k = 0; k < num_of_elements_found; k++){ // stores the valid roots.
			string tempstr;
			tempstr = mpz_get_str(NULL, 10, valid_roots[k]);
			all_valid_roots.push_back(tempstr);
			tempstr.clear();
		}
	}
	return all_valid_roots;
}
//**********************************************************************
// - Function description: given an array of bigintegers representing Bloom filters, it blinds and returns the blinded ones.
// Used in outsourcing phase, so no counter (to regenerate the most recent blinding factors) is needed.

bigint* Client::blind_BFs(bigint* bf, int bf_size, bigint BF_key, int bit_size, bigint pr_moduli){
	//genrates table size pseudorandom values using the BF_key.
	bigint* blinded_BF;
	blinded_BF = (mpz_t*)malloc(bf_size * sizeof(mpz_t));
	gmp_randstate_t rand_x;
	gmp_randinit_default(rand_x);
	gmp_randseed(rand_x, BF_key);
	bigint blinding_fac;
	for(int i = 0; i < bf_size; i++){
		mpz_init(blinding_fac);
		mpz_urandomb(blinding_fac, rand_x, bit_size);
		mpz_init(blinded_BF[i]);
		mpz_add(blinded_BF[i], blinding_fac, bf[i]); // blinds each biginteger representing a bloom filters.
		mpz_mod(blinded_BF[i], blinded_BF[i], pr_moduli);
	}
	mpz_clear(blinding_fac);
	gmp_randclear(rand_x);
	return blinded_BF;
}
//**********************************************************************
// - Function description: unblinds a biginteger representing a Bloom filters.

bigint* Client::unblind_BF(bigint BF, int  indx , bigint BFkey, bigint BF_counterkey, int bit_size, bigint pr_moduli){
	// regenertes the corresponding blinding factor.
	double start_1 = clock();
	bigint *unblinded_BF, tmp_key;
	unblinded_BF = (mpz_t*)malloc(1 * sizeof(mpz_t));
	gmp_randstate_t rand_y, rand_x, rand3;
	gmp_randinit_default(rand3);
	gmp_randinit_default(rand_y);
	gmp_randinit_default(rand_x);
	gmp_randseed(rand_y, BFkey);
	gmp_randseed(rand_x, BF_counterkey);
	bigint bf_key, bf_ck;
	double end_1;
	for(int i = 0; i <= indx; i++){
		mpz_init(bf_key);
		mpz_urandomb(bf_key, rand_y, bit_size);
		mpz_init(bf_ck);
		mpz_urandomb(bf_ck, rand_x, bit_size);
		if(i == 0){end_1 = clock();}
	}
	double temp = end_1 - start_1;
	double start_2 = clock();
	gmp_randstate_t rand1;
	gmp_randinit_default(rand1);
	bigint bld_factor;
	mpz_init(bld_factor);
	if(counter[indx] == 0){
		mpz_init_set(bld_factor, bf_key);
	}
	else{
		gmp_randseed(rand3, bf_ck);
		for(int k = 0; k < counter[indx]; k++){
			mpz_init(tmp_key);
			mpz_urandomb(tmp_key, rand3, bit_size);
		}
		mpz_add(bld_factor, tmp_key, bf_key);
		mpz_mod(bld_factor, bld_factor, pr_moduli);
	}
		mpz_init(unblinded_BF[0]);
		mpz_sub(unblinded_BF[0], pr_moduli, bld_factor);
		mpz_add(unblinded_BF[0], unblinded_BF[0], BF); // unblinds the biginteger representing a bloom filters.
		mpz_mod(unblinded_BF[0], unblinded_BF[0], pr_moduli);
		mpz_clear(bld_factor);
		double end_2 = clock();
		temp += end_2 - start_2;
		float f_temp = temp/(double) CLOCKS_PER_SEC;
		cout<<"\nIn update-- unblind_BF time: "<<f_temp<<endl;
		return unblinded_BF;
}
//**********************************************************************
// - Function description: blinds a biginteger representing a Bloom filter.

bigint* Client::blind_BF(bigint BF, int  indx , bigint BFkey, bigint BF_counterkey, int bit_size, bigint pr_moduli){
	// generates a fresh blinding factor.
	double start_1 = clock();
	double temp,end_1;
	bigint *blinded_BF, bf_key, bf_ck;
	blinded_BF = (mpz_t*)malloc(1 * sizeof(mpz_t));
	gmp_randstate_t rand_y,rand_x;
	gmp_randinit_default(rand_y);
	gmp_randinit_default(rand_x);
	gmp_randseed(rand_y,BFkey);
	gmp_randseed(rand_x,BF_counterkey);
	for(int i = 0; i <= indx; i++){
		mpz_init(bf_key);
		mpz_urandomb(bf_key,rand_y, bit_size);
		mpz_init(bf_ck);
		mpz_urandomb(bf_ck, rand_x, bit_size);
		if(i == 0){end_1 = clock();}
	}
	temp = end_1 - start_1;
	double start_2 = clock();
	gmp_randstate_t rand1;
	gmp_randinit_default(rand1);
	bigint bld_factor;
	if(counter[indx] == 0){
		mpz_init_set(bld_factor, bf_key);
	}
	else{
		bigint tmp_key;
		mpz_init(bld_factor);
		gmp_randstate_t rand3;
		gmp_randinit_default(rand3);
		gmp_randseed(rand3, bf_ck);
		for(int k = 0; k < counter[indx]; k++){
			mpz_init(tmp_key);
			mpz_urandomb(tmp_key, rand3, bit_size);
		}
		mpz_add(bld_factor, tmp_key, bf_key);
		mpz_mod(bld_factor, bld_factor, pr_moduli);
	}
		mpz_init(blinded_BF[0]);
		mpz_add(blinded_BF[0], bld_factor, BF); // blinds the biginteger representing a bloom filters.
		mpz_mod(blinded_BF[0], blinded_BF[0], pr_moduli);
		double end_2 = clock();
		temp += end_2 - start_2;
		float f_temp = temp / (double) CLOCKS_PER_SEC;
		cout<<"\n In update - blind_BF time:"<<f_temp<<endl;
		return blinded_BF;
}
//**********************************************************************
