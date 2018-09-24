# U-PSI
U-PSI: An Efficient Updatable Delegated Private Set intersection

* An efficient protocol that allows clients independently prepare and outosurce their private data
to a cloud server. Then, they can efficiently update or run private set intersection on the outsourced data.

# Dependencies
 * GMP: https://gmplib.org/
 * Cryptopp: https://www.cryptopp.com
 * NTL: https://www.shoup.net/ntl
 * Bloom filter: http://www.partow.net/programming/bloomfilter/index.html

# Runnig a Test
* First, clone the above libraries, and the U-PSI file. Then, install the libraries and unzip the U-PSI file. Next:

    <>cd U-PSI
    
* run: g++  -c  Rand.cpp -c Hashtable.cpp -c Polynomial.cpp -c Server.cpp -c Client.cpp
