MT 3D Inversion Workshop: Dublin Secret Model 1 (DSM1)

The 3D inversion test data set consists of 100 sites which are equally distributed over an area of about 50 km x 50 km. There are 10 profiles with 10 sites each. The spacing between the sites and also between the profiles is 5 km.

The meaning of the site names is as follows: letters for the profile - A (most northern profile) to J (most southern profile) and numbers for the sites along the profile - 01 (most western site) to 10 (most eastern site). The file ‘sitelocation.txt’ is a list of all the site names and site coordinates.

Each site has its own data file (A01.asc to J10.asc) which contains two header lines (1st: site name, easting and northing; 2nd: header of the data columns) and a data block. The data block consists of 9 columns: 1st column contains the periods in s and is followed by 8 columns - one for each impedance element.

The period range of the data set is 0.56 s to 10 000 s and the expected resistivity range of the structure is 1 Ωm - 1000 Ωm.

The impedance values (in SI units - Ω) are based on the eiωt convention and the assumption that x is pointing North and y is pointing East.

It is your choice how many and which sites and periods you include in your inversion (as it would be for real data). 
