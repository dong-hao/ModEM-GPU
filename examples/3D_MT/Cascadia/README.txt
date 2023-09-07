This example is provided for the convenience of the users, and makes it possible to explore ModEM inversion capabilities on the real 
data. The model and data set presented in this example are published in the following reference:

Patro, P. K., and G. D. Egbert (2008), Regional conductivity structure of Cascadia: Preliminary results from 3D inversion of USArray transportable array magnetotelluric data, Geophys. Res. Lett., 35, L20311, doi:10.1029/2008GL035326. 

Please respect the authors' copyright.

The published inversion was performed using WSINV3D. For your convenience, I have also provided a script to convert Weerachai 
Siripunvaraporn's data format to the one used in ModEM (readWSdata.pl). However, in using this script please be aware that the site names 
and geographic locations inserted in the data file will not be the true values, since this information is absent from WS data files.

I have provided four data files that can be used directly with the Modular System:

cascad_errfl5_published.dat  - impedances for 8 periods identical to the published data with 5% error floor
cascad_errfl5.dat            - all data (impedance + Hz) recreated independently for 10 periods (5% errors)
cascad_errfl3.dat            - all data (impedance + Hz) recreated independently for 10 periods (3% errors)
cascad_noerrfl.dat           - all data (impedance + Hz) recreated independently for 10 periods

Except for the data type headers, each line in the data files is an independent value. If you delete any of the lines, the program should 
still work. Please let me know if that is not the case.
Note that the files recreated independently also contain the names and geographic coordinates of the sites. These values are not directly 
used by the inversion at present, and can be replaced by stubs (such as in cascad_errfl5_published.dat).

The two prior models that are provided are in WS format, except they've been converted to natural log (see model headers).
The dimensions are 48x46x34 and 80x78x34, respectively. The models include the ocean, and the model covariance files are also provided.

A. Kelbert, Jan 2011
