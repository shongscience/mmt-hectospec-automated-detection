# mmt-hectospec-automated-detection
Automated detection routines written in IDL for MMT/Hectospec (Hong et al. 2014, PASP, 126, 1048)


1. Description : 

  This program depends on "hsred" and some hidden IDL routines, as many IDL routines have similar issues. This is the reason why python is getting more common in astronomy. Hence, this repository might be for back-up purpose, rather sharing purpose, though some people still can try this. 

2. How this works : 

  To reduce MMT/Hectospec data, there are two pipelines, one from Harvard/cfa and the other from HSRED https://www.mmto.org/node/536. I have used HSRED and added my own IDL routines to this HSRED package for implementing my detection algrothim.   
  
  2.1 Install : 
  
  Copy all IDL files in ~/spec2d/ to ~/hsred/idl/spec2d/, where HSRED is installed.
  
  2.2 How to use : 
  
  We assume that we have a fully reduced data from HSRED, "ex.fits", in "~/work/" directory. This ex.fits is supposed to have 300 spectra.
 
    A. Compile my "shong" files
      IDL> @compileAllShong.pro
    B. aa
      IDL> 
  

3. 
