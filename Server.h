/*

- Server-side computation of the U-PSI protocol

*/

/*
 Variables description:

        xpoints: x-coordinates generated by the server.
	pu_moduli: allows to define the field as F_{pubmoduli}.
	db_size: current length of server-side database.
	xpoint_size: number of x-coordinates.
	pub_moduli_bitsize: bit-size of pu_moduli.
	table_size: the hash table length.
	max_setsize: sets size upper bound.
	NoElem_in_bucket: each bin's capacity.

*/
//**********************************************************************

#include"Polynomial.h"

//**********************************************************************
// - Description: a request message for PSI computation. It is sent from the result recipient client to the authorizer client.

struct CompPerm_Request{
	string id;
	bigint**r;//blinded blinding factors
	bigint label_key_;
	bigint shuffle_key_;
};
//**********************************************************************
// - Description: a message for granting PSI computation. It is sent from the authorizer client to the server.

struct GrantComp_Info{
	string* id;// the result reciepent id is in id[0].
	bigint seed;// temporary key: tk
	bigint** pm; // permutation map
};
//**********************************************************************
// - Description: a message containing the result computed by the server and it is sent to the result recipient.

struct Server_Result{
	bigint** result;// an array of permuted bins of a hash table.
	bigint* BF;// an array of blinded Bloom filters.
};
//**********************************************************************
// - Description: A client's dataset maintained by the server.

struct Client_Dataset{
	Polynomial* poly;// an array of polynomials
	bigint* labels;// an array of labels
	bigint* BF;//array of blinded bloom filters
};
//**********************************************************************

class Server{

public:
	Server(int num_xpoints, int dbs_size, int pu_mod_bitsize, int maxSetsize, int NoEl_bucket, int tb_size);
	Server_Result* compute_result (GrantComp_Info * grantComp_info, bigint tmp_key_);
	bigint* get_xpoints(int& size);
	bigint* get_pubModuli();
	bigint* get_client_bin(bigint label, string ID, bigint& bf, int& indx);
	int get_maxSetsize();
	int get_NoElem_in_bucket();
	int get_table_size(){return table_size;}
	int get_pubModuli_bitsize(){return pub_moduli_bitsize;}
	void store_poly(Client_Dataset&);
	void update_client_bin(bigint* vals, bigint label, string id, bigint bbf);

private:
	void set_db(int index, Client_Dataset &p);
	void update_db(bigint* vals, bigint label, string id, bigint bbf);
	void regen_PRpolys(bigint key_, bigint **&w_A_, bigint **&w_B_, bigint **&a_, bigint**& tmp_bl_, bigint tmp_key_);
	Client_Dataset get_db(int index){return db[index];}
	int find_db_index(string id);
	//variables
	bigint * xpoints, * pu_moduli ;
	Client_Dataset *db;
	int db_size, xpoint_size, count, table_size, pub_moduli_bitsize, max_setsize, NoElem_in_bucket;
};
//**********************************************************************