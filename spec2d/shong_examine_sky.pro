;+
; NAME:
;   shong_examine_sky.pro
;
; PURPOSE:
;   Examine the extracted file. Modified by S. Hong. for my customized purpose. 
;
; CALLING SEQUENCE:
;  hs_page_file, filename, nstart=, plugmap=, ivar=, nsmooth=, xrange=, 
;
; INPUTS:
;  filename 6   - filename for for the Hectospec data
;
; OPTIONAL KEYWORDS:
;   nstart - starting aperture
;   nsmooth - number of pixels to smooth over
;   xrange - plotting xrange
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA

;------------------------------------------------------------------------------

PRO shong_examine_sky, filename, nstart=nstart, plugmap=plugmap, ivar=ivar, $
                  zmap=zmap, autoline=autoline, smin = smin, synth=synth, $
                  nsmooth=nsmooth, psym=psym,title=title, xrange=xrange, $
                  clipfact=clipfact, dosnr=dosnr, pseudoflux=pseudoflux, $
                  synthmag=synthmag, thick=thick
  
  if not keyword_set(xrange) then xrange=[3800,8500]
  if not keyword_set(psym) then psym=0
  if not keyword_set(smin) then smin = 3650
  if NOT keyword_set(nstart) then nstart = 0
  if NOT keyword_set(clipfact) then clipfact =0.6
  
  ;;Read the file
  lam = mrdfits(filename, 0)
  object = mrdfits(filename, 1)
  ivar = mrdfits(filename, 2)
  mask = mrdfits(filename, 4)*0.0
  original_mask = mrdfits(filename, 4)
  sky_one_struct = mrdfits(filename, 6)
  sky_two_struct = mrdfits(filename, 8)
  plugmap = mrdfits(filename, 5)

  sky_one = bspline_valu(lam,sky_one_struct)
  sky_two = bspline_valu(lam,sky_two_struct)
     
     
  if keyword_set(pseudoflux) then begin
     fcalfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
     tcalfile = getenv('HSRED_DIR') + '/etc/average_tell.fits'
     
     fcal = mrdfits(fcalfile, 1)
     tcal = mrdfits(tcalfile, 1)
     flux_factor = bspline_valu(alog10(lam), fcal)
     tell_factor = bspline_valu((lam), tcal)
     divideflat, object, flux_factor, invvar=ivar
     divideflat, object, tell_factor, invvar=ivar
  endif
  
  
  test = ''
  splog, '>>>>>>>>>>TYPE HELP for Commands<<<<<<<<<'
  
  insmooth=0
  shong_object_junk = 0.0D
  shong_ivar_junk = 0.0D
  shong_gaussian_kernel_sigma = 1.0D
  shong_gaussian_ker = psf_gaussian([1.0,floor(2.5*shong_gaussian_kernel_sigma),shong_gaussian_kernel_sigma],$
		npixel=floor(5.0*shong_gaussian_kernel_sigma))
  shong_gaussian_ker = shong_gaussian_ker/total(shong_gaussian_ker) ; normalize 
  shong_gaussian_ker2 = shong_gaussian_ker * shong_gaussian_ker 

  shong_isYRangeManual = 0
  yminn = 0.0D
  ymaxx = 10.0D

  i = nstart
  goMasked = 0; show unmasked data as default
  goUnsmooth = 1; show un-smoothed profile as default
  goGaussianSmooth = 0; Showing top hat smooth as default: 1 means gaussian smooth or top hat smooth
  stop = 'go'
  
  if keyword_set(dosnr) then hs_snr, lam, object, plugmap, snr=snr
  if not keyword_set(dosnr) then snr = fltarr(n_elements(object(0,*)))
  
  if NOT keyword_set(nsmooth) then nsmooth=0
  while stop ne 'stop' do begin
	  
     ;;Mask pixels were the ivar is lt 0.6 the smoothed value
     temp = ivar[*,i] / smooth(ivar[*,i], 21)
     newmask = temp lt clipfact
     
     junk = min(abs(lam[*,i]-5588),pix)
     mask[(pix-10):(pix+10),i]=mask[(pix-10):(pix+10),i]+$
        newmask[(pix-10):(pix+10)]
     
     junk = min(abs(lam[*,i]-5577), pix)
     
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     junk = min(abs(lam[*,i]-6300), pix)
     
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     
     junk = min(abs(lam[*,i]-6363), pix)

     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     
     pixels = where(mask[*,i] eq 0)  ; un-masked indice
     masked_pixels = where(mask[*,i] ne 0)  ; masked indice

	;;;; For more functions ;; added by Hong.
	;;;; MASK ************
	shong_object_junk = total(object[pixels,i])/n_elements(pixels)
	shong_ivar_junk = total(ivar[pixels,i])/n_elements(pixels)
	unmasked_object = object[*,i]
	unmasked_var = 1.0D/ivar[*,i]
	masked_object = object[*,i]
	masked_var = 1.0D/ivar[*,i]
	masked_object[masked_pixels] = shong_object_junk ; assign a mean for the masked object pixels
    masked_var[masked_pixels] = 1.0D/shong_ivar_junk ; assing a mean for ivar masked
;	print, "unmasked objects: ", total(object[pixels,i])," number of unmasked pixels: ", n_elements(pixels), " avg: ", shong_object_junk
;	print, "unmasked ivar: ", total(ivar[pixels,i])," number of unmasked pixels: ", n_elements(pixels), " avg: ", shong_ivar_junk
;	print, "total objects: ", total(object[*,i])," number of unmasked pixels: ", n_elements(object[*,i]), " avg: "
;	print, "total ivar: ", total(ivar[*,i])," number of unmasked pixels: ", n_elements(ivar[*,i]), " avg: "

	;;;; SMOOTH ************
    shong_gaussian_ker = psf_gaussian([1.0,floor(2.5*shong_gaussian_kernel_sigma),shong_gaussian_kernel_sigma],$
    		npixel=floor(5.0*shong_gaussian_kernel_sigma))
    shong_gaussian_ker = shong_gaussian_ker/total(shong_gaussian_ker) ; normalize 
    shong_gaussian_ker2 = shong_gaussian_ker * shong_gaussian_ker 

	sky_one_tsmooth = smooth(sky_one,nsmooth)
	sky_two_tsmooth = smooth(sky_two,nsmooth)
	sky_one_gsmooth = convol(sky_one,shong_gaussian_ker,/center)
	sky_two_gsmooth = convol(sky_two,shong_gaussian_ker,/center)
	masked_object_tsmooth = smooth(masked_object,nsmooth)
	masked_object_gsmooth = convol(masked_object,shong_gaussian_ker,/center)
	unmasked_object_tsmooth = smooth(unmasked_object,nsmooth)
	unmasked_object_gsmooth = convol(unmasked_object,shong_gaussian_ker,/center)

	masked_var_tsmooth = smooth(masked_var,nsmooth)/double(nsmooth) ; additional 1/N factor for normalization
	masked_var_gsmooth = convol(masked_var,shong_gaussian_ker2,/center)
	unmasked_var_tsmooth = smooth(unmasked_var,nsmooth)/double(nsmooth) ; the same with the above
	unmasked_var_gsmooth = convol(unmasked_var,shong_gaussian_ker2,/center)



     l1 =  min(abs(lam[pixels,i]-xrange(0)),l)
     m1 =  min(abs(lam[pixels,i]-xrange(1)),m)
     
     l = l(0)
     m = m(0)
     
     
     if keyword_set(title) then title1=title(i)
     if not keyword_set(title) then title1=''



     junk1 =  min(abs(lam[pixels,i]-xrange(0)),xmini)
     junk2 =  min(abs(lam[pixels,i]-xrange(1)),xmaxi)
 

	 ymax = max(sky_one[xmini:xmaxi]) * 1.1
	 ymin = min(sky_one[xmini:xmaxi])
	 diff = ymax-ymin

     djs_plot, lam[*,i], sky_one,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
		
     djs_oplot, lam[*,i], sky_two, color='blue'
     ;djs_oplot, lam[*,i], 1.0D*sqrt(masked_var), color='red'
   
    	; djs_oplot, lam[pixels,i], 50*mask, color='blue'






		;; showing skylines
           skylines= [5577.,6300.64]
           skylabels = ['OI(5577)','OI(6301)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='magenta', linestyle=2
              djs_xyouts, (1+0.0002)*skylines(j), ymax-diff*(0.2+0.02*(-1)^j), $
                 skylabels(j), charsize=2, color='green'
           endfor



     if keyword_set(plugmap) then begin
        
        djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
           ymin + (ymax-ymin)*0.90  ,  plugmap[i].objtype, $
           color='yellow', thick=thick, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.05,$
           ymin + (ymax-ymin)*0.85, $
           'snr = ' + strn(snr(i)), $
           color='yellow', thick=thick, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
           ymin + (ymax-ymin)*0.95, $
           'icode = ' + strn(plugmap[i].icode), $
           color='yellow', thick=thick, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
           ymin + (ymax-ymin)*0.9, $
           'ap mag = ' + strn(plugmap[i].rapmag), $
           color='yellow', thick=thick, charsize=2, /isolatin1
        
        ;;Here, I am going to do something a bit strange and go ahead and
        ;;integrate the spectrum at the redshift observed.  This only works 
        ;;if the /synthmag is set
        
        if keyword_set(synthmag) then begin
           
           k_load_filters, 'sdss_r0.par', nlam, lam1, trans
           integrand = object[*,i]*1e-17
           intertrans = interpol(trans, lam1, lam[*,i])
           k = where(intertrans lt 0)
           intertrans(k) = 0
           rmag = int_tabulated(lam[*,i], integrand*intertrans)
           djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.45,$
              ymin + (ymax-ymin)*0.95, $
              'synth mag = ' + strn(-2.5*alog10(rmag)), $
              color='yellow', thick=thick, charsize=2, /isolatin1
        endif
        
        
     endif


     
     loop = 0 
     while loop eq 0 do begin
        read, prompt='Command: ', test
        
        if strupcase(test) eq 'N' then begin
           i = i+ 1
           loop = 1
        endif
        
        if strupcase(test) eq 'P' then begin
           i = i - 1
           loop = 1
        endif
        
        if strupcase(test) eq 'C' then begin
           i = i 
           loop = 1
        endif

        if strupcase(test) eq 'STOP' then begin
           loop = 1
           stop = 'stop'
        endif
 


        if strupcase(test) eq 'X' then begin
           splog, 'Enter New xrange'
           read, prompt='xmin : ', xmin
           read, prompt='xmax : ', xmax
		   xrange(0) = xmin
		   xrange(1) = xmax
		   i = i
           loop = 1
        endif

        if strupcase(test) eq 'Y' then begin
           splog, 'Enter New yrange'
           read, prompt='ymin : ', yminn
           read, prompt='ymax : ', ymaxx
		   i = i
           loop = 1
        endif
	    if strupcase(test) eq 'YON' then begin
		   shong_isYRangeManual = 1
		   i = i
           loop = 1
        endif
	    if strupcase(test) eq 'YOFF' then begin
		   shong_isYRangeManual = 0
		   i = i
           loop = 1
        endif

       if strupcase(test) eq 'SHOWLOG' then begin
           i = i 
           loop = 1
			 ;;; print the related quantities
			 print, '>>Current xrange = ', xrange, '; ymin, ymax = ', ymin, ymax
			 print, '>>Current tsmooth = ', nsmooth, ' ;gsmooth = ', shong_gaussian_kernel_sigma,$
				 ' ;unsmooth = ', goUnsmooth, ' ;Gaussian Smooth = ', goGaussianSmooth
       endif
 

       if strupcase(test) eq 'HELP' then begin
           i = i 
           loop = 1
		  splog, 'Commands'
		  splog, 'N - next'
		  splog, 'P - previous'
		  splog, 'Stop - stop'
		  splog, 'A - plot  absorption lines'
		  splog, 'E - plot  emission lines'
		  splog, 'C - clear lines'
		  splog, 'MASK - show masked profile'
		  splog, 'UNMASK - show unmasked profile'
		  splog, 'GSMOOTH - show gaussian smoothed '
		  splog, 'TSMOOTH - show top-hat smoothed '
		  splog, 'UNSMOOTH - show un-smoothed '
		  splog, '*X - change xrange'
		  splog, '*Y - set manual yrange'
		  splog, '*YON - On manual yrange'
		  splog, '*YOFF - Off manual yrange'
		  splog, '*S - change nsmooth'
		  splog, '*SUP - smooth up +'
		  splog, '*SDW - smooth down -'
		  splog, '*G - change sigma of Gaussian smooth'
		  splog, '*GUP - gaussian smooth up +'
		  splog, '*GDW - gaussian smooth down -'
		  splog, '***SHOWSKYMAX - Show the maximum position of sky OI6301 using unmasked-var-gsmooth profile'
		  splog, '***SHOWLOG - Show internal variables'
       endif


   
         if strupcase(test) eq 'S' then begin
           splog, 'Enter New Top Hat Smooth Width'
           read, prompt='nsmooth : ', insmooth
		   nsmooth = insmooth
		   i = i
           loop = 1
        endif

          if strupcase(test) eq 'SUP' then begin
           splog, 'Top Hat Smooth Width Up'
		   nsmooth=nsmooth+2
		   i = i
           loop = 1 
        endif
        
           if strupcase(test) eq 'SDW' then begin
           splog, 'Top Hat Smooth Width Down'
		   nsmooth=nsmooth-2
		   if (nsmooth lt 0) then begin
				nsmooth=0			
		   endif
		   i = i
           loop = 1 
        endif

          if strupcase(test) eq 'G' then begin
           splog, 'Enter Gaussian Smooth Kernel Width'
           read, prompt='nsmooth : ', isigma
		   shong_gaussian_kernel_sigma = double(isigma)
		   i = i
           loop = 1
        endif

          if strupcase(test) eq 'GUP' then begin
           splog, 'Gaussian Smooth Width Up'
		   shong_gaussian_kernel_sigma = shong_gaussian_kernel_sigma + 1.0D
		   i = i
           loop = 1 
        endif
        
           if strupcase(test) eq 'GDW' then begin
           splog, 'Gaussian Smooth Width Down'
		   shong_gaussian_kernel_sigma = shong_gaussian_kernel_sigma - 1.0D
		   if (shong_gaussian_kernel_sigma lt 1.0) then begin
				shong_gaussian_kernel_sigma = 1.0			
		   endif 
		   i = i
           loop = 1 
        endif
       
         if strupcase(test) eq 'MASK' then begin
           i = i 
           loop = 1
		   goMasked = 1
        endif

          if strupcase(test) eq 'UNMASK' then begin
           i = i 
           loop = 1
		   goMasked = 0
        endif
 
          if strupcase(test) eq 'UNSMOOTH' then begin
           i = i 
           loop = 1
		   goUnsmooth = 1 
        endif
        if strupcase(test) eq 'GSMOOTH' then begin
           i = i 
           loop = 1
		   goUnsmooth = 0 
		   goGaussianSmooth = 1
        endif
        if strupcase(test) eq 'TSMOOTH' then begin
           i = i 
           loop = 1
		   goUnsmooth = 0 
		   goGaussianSmooth = 0
        endif
 

        if strupcase(test) eq 'SHOWSKYMAX' then begin
		ii=0
		for ii=0,299 do begin
			   	unmasked_var = 1.0D/ivar[*,ii]				
				unmasked_var_gsmooth = convol(unmasked_var,shong_gaussian_ker2,/center)

			   lam_min = min(abs(lam[*,ii]-6296),pix_lammin) ; 6290 
			   lam_max = min(abs(lam[*,ii]-6305),pix_lammax) ; 6310 
			   var_max = max(unmasked_var_gsmooth[pix_lammin:pix_lammax],pix_max) 
			   pix_max = pix_max + pix_lammin

				;djs_oplot,[lam[pix_max,i]],[unmasked_var_gsmooth[pix_max]],psym=2,thick=3,symsize=3,color='red'

			  ; oplot,[lam[pix_max,ii]],[sqrt(unmasked_var_gsmooth[pix_max])],psym=2,thick=3,symsize=3,color=100

			   print, plugmap[ii].objtype,lam[pix_max,ii],sqrt(unmasked_var_gsmooth[pix_max]),shong_gaussian_kernel_sigma
				
			   junky = mpfitpeak(lam[pix_lammin:pix_lammax],unmasked_var_gsmooth[pix_lammin:pix_lammax],outpara)
			   print, plugmap[ii].objtype,outpara(1),outpara(0)+outpara(3),outpara(2),shong_gaussian_kernel_sigma
			   

	;		   print,"lammin_i, lammax_i,lammin, lammax, idx, lam, stdev: ",$
	;				pix_lammin,pix_lammax,lam[pix_lammin,i], lam[pix_lammax,i], pix_max,lam[pix_max,i],unmasked_var_gsmooth[pix_max]
		 endfor
		
           i = i
           loop = 1
 		endif

     if strupcase(test) eq 'SHOWDLAMBDA' then begin
		ii=0
		sum = 0.0D
		leng=n_elements(lam[*,i])
		leng = leng - 2
		dlambda = 0.0D
		print, "Delta Lambda:"
		for ii=0, leng do begin
			dlambda = abs(lam[ii,i]-lam[ii+1,i])	
			sum= sum + dlambda
			print,dlambda ; 6290 				
		endfor
		print,"Deleta Lambda Average: ", (sum/double(leng+2))
        i = i
        loop = 1
 	endif





       
        if strupcase(test) eq 'A' then begin
           splog, 'Enter the Redshift you would like plotted'
           read, prompt='Redshift : ', z
           
           lines= [3797.90,3834.00,3889.00,3933.70,3968.50,4101.70,$
                   4304.40,4340.48,4861.34,5175.36,5268.98,5892.50,$
                   4172.00,4226.74]
           labels = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
                     'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
                     'Na',  'CaI',  'CaI']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='magenta', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor
        endif
        
        if strupcase(test) eq 'E' then begin
           splog, 'Enter the Redshift you would like plotted'
           read, prompt='Redshift : ', z
           
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='yellow', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor 
        endif
        
     endwhile
     
  endwhile
  
END

