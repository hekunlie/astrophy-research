1. ''invalid pointer'' or segment fault： the index accesses the regoin out of bound (non-allocated location in memory).
     
     2018-12-22: Reading a hdf5 contains an array which is larger than that has been declared in the program 
     will give rise to the ''invalid pointer'' bug when it comes to delete the array on the end of the program.