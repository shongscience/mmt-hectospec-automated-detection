# mmt-hectospec-automated-detection
Automated detection routines written in IDL for MMT/Hectospec (Hong et al. 2014, PASP, 126, 1048)


1. Description : 

  This program depends on "hsred" and some hidden IDL routines, as many IDL routines have similar issues. This is the reason why python is getting more common in astronomy. Hence, this repository might be for back-up purpose, rather for share, though some people still can try this. 

2. How this works : 

  To reduce MMT/Hectospec data, there are two pipelines, one from Harvard/cfa and the other from HSRED https://www.mmto.org/node/536. I have used HSRED and added my own IDL routines to this HSRED package for implementing my detection algrothim.   
  
  2.1 Install : 
  
  Copy all IDL files in ~/spec2d/ to ~/hsred/idl/spec2d/, where HSRED is installed.
  
  2.2 How to use : 
  
  We assume that we have a fully reduced data from HSRED, "ex.fits", in "~/work/" directory. This ex.fits is supposed to have 300 spectra.
 
    A. Compile my "shong" files : 
    
      IDL> @compileAllShong.pro

    B. Lauch the interactive investigating mode : (This is a hacked verson of hs_page_file.pro) 
      
      IDL> shong_page_file, 'ex.fits' 
      
    C. type "help" to see all possible menus for shong_page_file mode :
    
      "x" will set the searching range; and many other useful commands. 
      "dumpdetecnosmooth" will dump all potential detections on the "./detection" folder. 
      "stop" to exit this "shong_page_file" mode. 
      
    D. move on to the ./detection/ folder and run mmtExtractDump.pro : This will generate all pictures shown in Hong et al. 2014
    
      IDL>.comp mmtExtractDump.pro
      IDL> mmtExtractDump,'pm1lae.list’,dldp=1.1,histymax=10  ;; dldp is lambda per pixel
      
      ** There are hard-wired montecarlo simulation data 
      of reliability and completeness of detections, applicable to MMT/Hectospec only. 

    E. The attached "eps" files are automatically generated. 

