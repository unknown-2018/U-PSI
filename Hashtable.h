/*

- Hash table class used in the U-PSI protocol.

*/


#include"Rand.h"

//**********************************************************************

class Hashtable{

public:
	Hashtable ( int elem_in_bucket, bigint* elemen, int elem_size,int table_size);
	bigint* get_bucket(int index);
private:
	bigint **T; // hash table
};
