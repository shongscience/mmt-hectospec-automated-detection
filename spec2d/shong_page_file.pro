;+
; NAME:
;   shong_page_file
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
;   January 2012 - Written by S Hong - NOAO
;------------------------------------------------------------------------------

PRO shong_page_file, filename, nstart=nstart, plugmap=plugmap, ivar=ivar, $
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
  sky_one = mrdfits(filename, 6)
  sky_two = mrdfits(filename, 8)
  plugmap = mrdfits(filename, 5)
  skyoneflux = bspline_valu(lam, sky_one)
     
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
  shong_isXRangeManual = 0
  shong_isShowSkyEM = 1
  shong_isShowSky = 0
  shong_isRunDetect = 0
  shong_isDumpDetect = 0
  shong_isDumpDetectNoSmooth = 0
  shong_isDumpNoDetectNoSmooth = 0
  yminn = -50.0D
  ymaxx = 100.0D
  xminn = 4000.0D
  xmaxx = 5000.0D

  i = nstart
  goMasked = 0; show unmasked data as default
  goUnsmooth = 1; show un-smoothed profile as default
  goGaussianSmooth = 0; Showing top hat smooth as default: 1 means gaussian smooth or top hat smooth
  stop = 'go'
  
  if keyword_set(dosnr) then hs_snr, lam, object, plugmap, snr=snr
  if not keyword_set(dosnr) then snr = fltarr(n_elements(object(0,*)))
  
  if NOT keyword_set(nsmooth) then nsmooth=0


  
  window, 0, xpos=2100, ypos=1000 ;make the plots on window 0 not on window 10 
  !p.multi = 0 
  while stop ne 'stop' do begin
	
	wset,0
  	!p.multi = 0 

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
     
     ymax = max(smooth(object[pixels(l):pixels(m),i],nsmooth)) *1.1
     ymin = min(smooth(object[pixels(l):pixels(m),i],nsmooth))
     diff = ymax-ymin    

     
     if keyword_set(title) then title1=title(i)
     if not keyword_set(title) then title1=''



     junk1 =  min(abs(lam[pixels,i]-xrange(0)),xmini)
     junk2 =  min(abs(lam[pixels,i]-xrange(1)),xmaxi)


 
	 ymax_mo = max(masked_object[xmini:xmaxi]) * 1.4
	 ymin_mo = min(masked_object[xmini:xmaxi]) * 1.2
	 ymin_mv = min(sqrt(masked_var[xmini:xmaxi]))
	 ymax_mv = max(sqrt(masked_var[xmini:xmaxi])) * 5.4

	 ymax_uo = max(unmasked_object[xmini:xmaxi]) * 1.4
	 ymin_uo = min(unmasked_object[xmini:xmaxi]) * 1.2
	 ymin_uv = min(sqrt(unmasked_var[xmini:xmaxi]))
	 ymax_uv = max(sqrt(unmasked_var[xmini:xmaxi])) * 5.4

	 ymax_mog = max(masked_object_gsmooth[xmini:xmaxi]) * 1.4
	 ymin_mvg = min(sqrt(masked_var_gsmooth[xmini:xmaxi]))

	 ymax_uog = max(unmasked_object_gsmooth[xmini:xmaxi]) * 1.4
	 ymin_uvg = min(sqrt(unmasked_var_gsmooth[xmini:xmaxi]))

	 ymax_mot = max(masked_object_tsmooth[xmini:xmaxi]) * 1.4
	 ymin_mvt = min(sqrt(masked_var_tsmooth[xmini:xmaxi]))
	
	 ymax_uot = max(unmasked_object_tsmooth[xmini:xmaxi]) * 1.4
	 ymin_uvt = min(sqrt(unmasked_var_tsmooth[xmini:xmaxi]))

     ymax = max([ymax_mo,ymax_uo,ymax_mog,ymax_uog,ymax_mot,ymax_uot,ymax_uv,ymax_mv])
     ymin = min([ymin_mv,ymin_uv,ymin_mvg,ymin_uvg,ymin_mvt,ymin_uvt,ymin_uo,ymin_mo])
	 diff = ymax-ymin

	 if (shong_isYRangeManual) then begin
		ymax = ymaxx	
		ymin = yminn
     	diff = ymax-ymin    
	 endif
	 if (shong_isXRangeManual) then begin
		xrange(0) = xminn	
		xrange(1) = xmaxx
	 endif







	 if (goUnsmooth) then begin ;*********AA

	 if (goMasked) then begin ;**AA


     djs_plot, lam[*,i], masked_object,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
     djs_oplot, lam[*,i], 3.0D*sqrt(masked_var), color='red'
 
   		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif

  

	 endif else begin ;**AA
 

     djs_plot, lam[*,i], unmasked_object,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
		
     djs_oplot, lam[*,i], 3.0D*sqrt(unmasked_var), color='red'

   		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif





	 endelse ;**AA

	 endif else begin ; ********* AA

	 if (goGaussianSmooth) then begin ; ******BB

	 if (goMasked) then begin ;**AA

     djs_plot, lam[*,i], masked_object_gsmooth,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick	
     djs_oplot, lam[*,i], 3.0D*sqrt(masked_var_gsmooth), color='red'
 
   		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif

  


	 endif else begin ;**AA



     djs_plot, lam[*,i], unmasked_object_gsmooth,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
		
     djs_oplot, lam[*,i], 3.0D*sqrt(unmasked_var_gsmooth), color='red'

   		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif





	 endelse ;**AA

	 endif else begin ; ******BB

	 if (goMasked) then begin ;**AA


     djs_plot, lam[*,i], masked_object_tsmooth,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
		
     djs_oplot, lam[*,i], 3.0D*sqrt(masked_var_tsmooth), color='red'

   		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif





	 endif else begin ;**AA

     djs_plot, lam[*,i], unmasked_object_tsmooth,$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=thick
		
     djs_oplot, lam[*,i], 3.0D*sqrt(unmasked_var_tsmooth), color='red'

		if (shong_isShowSky) then begin
			;; showing skylines
     		djs_oplot, lam[*,i], 0.1*skyoneflux[*,i], color='magenta', linestyle=0
		endif



	 endelse ;**AA

	 endelse ;******BB


     endelse ; ********* AA



	if (shong_isShowSkyEM) then begin
		;; showing skylines
	       skylines= [4046.56,4165,4168,4358.34,4420,4423,$
				4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5889.95,5895.92,6300.64,6363.78]
	       skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
				'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
				'OI(5577)','NaI(5890)',$
				'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='magenta', linestyle=2
              djs_xyouts, (1+0.0002)*skylines(j), ymax-diff*(0.3+0.04*(j mod 5)), $
                 skylabels(j), charsize=2, color='green'
           endfor
	endif


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



	 if (shong_isRunDetect) then begin
		xminIdx = 0
		xmaxIdx = 0
		junk = min(abs(lam[*,i]-xrange(0)),xminIdx)
		junk = min(abs(lam[*,i]-xrange(1)),xmaxIdx)
		;print, "In showdetect: xmin_idx, xmax_idx = ", xminIdx, xmaxIdx
		shong_runDetection,plugmap[i],1.0,4.0,$
			lam[xminIdx:xmaxIdx,i], object[xminIdx:xmaxIdx,i], (sqrt(1.0D/ivar[xminIdx:xmaxIdx,i]))
	 endif

	 if (shong_isDumpDetect) then begin
		xminIdx = 0
		xmaxIdx = 0
		junk = min(abs(lam[*,i]-xrange(0)),xminIdx)
		junk = min(abs(lam[*,i]-xrange(1)),xmaxIdx)
		;print, "In showdetect: xmin_idx, xmax_idx = ", xminIdx, xmaxIdx
		
		ibackup = i
		for i=0L,299 do begin
			shong_dumpDetection,plugmap[i],1.5,4.0,$
				lam[xminIdx:xmaxIdx,i], object[xminIdx:xmaxIdx,i], (sqrt(1.0D/ivar[xminIdx:xmaxIdx,i])), goplot=1
			print,(i+1),"/300 Done..."
		endfor	
		i = ibackup

		shong_isDumpDetect = 0 ;; do not dump in the next time	
	 endif

	 if (shong_isDumpDetectNoSmooth) then begin
		xminIdx = 0
		xmaxIdx = 0
		junk = min(abs(lam[*,i]-xrange(0)),xminIdx)
		junk = min(abs(lam[*,i]-xrange(1)),xmaxIdx)
		;print, "In showdetect: xmin_idx, xmax_idx = ", xminIdx, xmaxIdx
		
		wholeXminIdx = xminIdx
		wholeXmaxIdx = xmaxIdx
		wholelam = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		wholeprofile = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		wholeivar = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		tempname=''
		ibackup = i


		for i=0L,299 do begin
			shong_dumpDetectionNoSmoothGoodFit,plugmap[i],1.5,4.0,$
				lam[xminIdx:xmaxIdx,i], object[xminIdx:xmaxIdx,i], (sqrt(1.0D/ivar[xminIdx:xmaxIdx,i])), goplot=0, interactive=0

			wholelam = lam[*,i]
			wholeprofile = object[*,i]
			wholeivar = ivar[*,i]
			tempname=strcompress(plugmap[i].objtype,/remove_all)+'specWholeProfile'+'.sav'
			tempname='./detection/'+tempname
			if direxist('detection') eq 0L then spawn, 'mkdir detection'
			save,wholeXminIdx,wholeXmaxIdx,wholelam,wholeprofile,wholeivar,filename=tempname
			print,(i+1),"/300 Done..."
		endfor	
		i = ibackup

		shong_isDumpDetectNoSmooth = 0 ;; do not dump in the next time	
	 endif


	 if (shong_isDumpNoDetectNoSmooth) then begin
		xminIdx = 0
		xmaxIdx = 0
		junk = min(abs(lam[*,i]-xrange(0)),xminIdx)
		junk = min(abs(lam[*,i]-xrange(1)),xmaxIdx)
		;print, "In showdetect: xmin_idx, xmax_idx = ", xminIdx, xmaxIdx
		
		wholeXminIdx = xminIdx
		wholeXmaxIdx = xmaxIdx
		wholelam = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		wholeprofile = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		wholeivar = dindgen(n_elements(lam[*,1])) ;; initialize wholelam
		tempname=''
		ibackup = i


		for i=0L,299 do begin
			shong_dumpNoDetectionNoSmoothGoodFit,plugmap[i],1.5,4.0,$
				lam[xminIdx:xmaxIdx,i], object[xminIdx:xmaxIdx,i], (sqrt(1.0D/ivar[xminIdx:xmaxIdx,i])), goplot=0, interactive=0

			wholelam = lam[*,i]
			wholeprofile = object[*,i]
			wholeivar = ivar[*,i]
			tempname=strcompress(plugmap[i].objtype,/remove_all)+'specWholeProfile'+'.sav'
			tempname='./nodetection/'+tempname
			if direxist('nodetection') eq 0L then spawn, 'mkdir nodetection'
			save,wholeXminIdx,wholeXmaxIdx,wholelam,wholeprofile,wholeivar,filename=tempname
			print,(i+1),"/300 Done..."
		endfor	
		i = ibackup

		shong_isDumpNoDetectNoSmooth = 0 ;; do not dump in the next time	
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
           read, prompt='xmin : ', xminn
           read, prompt='xmax : ', xmaxx
		   xrange(0) = xminn
		   xrange(1) = xmaxx
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

	    if strupcase(test) eq 'XON' then begin
		   shong_isXRangeManual = 1
		   i = i
           loop = 1
        endif
	    if strupcase(test) eq 'XOFF' then begin
		   shong_isXRangeManual = 0
		   xrange(0) = min(lam[*,i])
		   xrange(1) = max(lam[*,i])
		   i = i
           loop = 1
        endif




		if strupcase(test) eq 'SHOWSKYEM' then begin
		   shong_isShowSkyEM = 1
		   i = i
           loop = 1
        endif
		if strupcase(test) eq 'HIDESKYEM' then begin
		   shong_isShowSkyEM = 0
		   i = i
           loop = 1
        endif

		if strupcase(test) eq 'SHOWSKY' then begin
		   shong_isShowSky = 1
		   i = i
           loop = 1
        endif
		if strupcase(test) eq 'HIDESKY' then begin
		   shong_isShowSky = 0
		   i = i
           loop = 1
        endif


		if strupcase(test) eq 'RUNDETECT' then begin
		   shong_isRunDetect = 1
		   i = i
           loop = 1
    	endif

		if strupcase(test) eq 'ENDDETECT' then begin
		   shong_isRunDetect = 0
		   wdelete, 10
		   !p.multi = 0
		   i = i
           loop = 1
    	endif

		if strupcase(test) eq 'DUMPDETECT' then begin
		   shong_isDumpDetect = 1
		   i = i
           loop = 1
    	endif

		if strupcase(test) eq 'DUMPDETECTNOSMOOTH' then begin
		   shong_isDumpDetectNoSmooth = 1
		   i = i
           loop = 1
    	endif
		if strupcase(test) eq 'DUMPNODETECTNOSMOOTH' then begin
		   shong_isDumpNoDetectNoSmooth = 1
		   i = i
           loop = 1
    	endif



		if strupcase(test) eq 'FIND' then begin
			i = i
			itmp=0L
			iFound=0L
			searchname=''
			read,searchname, prompt='Input the number of LAE :'
			searchname='LAEA_'+strcompress(searchname,/remove_all)
			print,">> Searching ",searchname," ..."

			currentname=''
			for itmp=0L,299 do begin 
				currentname=strcompress(plugmap[itmp].objtype,/remove_all)
				print,"      searchname, currentname: ",searchname," ",currentname
				if searchname eq currentname then begin 
					iFound = itmp
					itmp = 5000
					print,"      Found: ",searchname," ",currentname
					i = iFound
				endif
				;itmp++
			endfor 
           loop = 1
    	endif

		if strupcase(test) eq 'PRINTLIST' then begin
			i = i
			itmp=0L
			iFound=0L
			outfilename=''
			outfilename=filename+'.list'
			openw, outunit, outfilename,/get_lun
			print,">> Printing the object name list in ",outfilename," ..."

			currentname=''
			for itmp=0L,299 do begin 
				currentname=strcompress(plugmap[itmp].objtype,/remove_all)
				print,"Index : ",itmp,"  Name : ",currentname," ",plugmap[itmp].ra," ",plugmap[itmp].dec," ",plugmap[itmp].fiberid
				printf,outunit,itmp," ",currentname," ",plugmap[itmp].objtype," ",plugmap[itmp].ra," ",plugmap[itmp].dec,$
					" ",plugmap[itmp].fiberid
				;itmp++
			endfor 
           	loop = 1
			
			free_lun, outunit 
    	endif






       if strupcase(test) eq 'SHOWLOG' then begin
           i = i 
           loop = 1
			 ;;; print the related quantities
			 print, '>>Current xrange = ', xrange, '; ymin, ymax = ', ymin, ymax, '; xmini, xmaxi',xmini, xmaxi
			 print, '>>Current tsmooth leng = ', nsmooth, ' ;gsmooth leng = ', shong_gaussian_kernel_sigma
			 print,	'>>Current flags>> ','Unsmooth = ', goUnsmooth, ' ; G Smooth = ', goGaussianSmooth, $
					' ;Masked = ', goMasked
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
		  splog, 'LAE - plot  emission lines '
		  splog, 'C - clear lines'
		  splog, 'MASK - show masked profile'
		  splog, 'UNMASK - show unmasked profile'
		  splog, 'GSMOOTH - show gaussian smoothed '
		  splog, 'TSMOOTH - show top-hat smoothed '
		  splog, 'UNSMOOTH - show un-smoothed '
		  splog, 'SHOWSKYEM - show sky emission lines '
		  splog, 'HIDESKYEM - hide sky emission lines '
		  splog, '*X - change xrange'
		  splog, '*Y - set manual yrange'
		  splog, '*XON - On manual xrange'
		  splog, '*XOFF - Off manual xrange'
		  splog, '*YON - On manual yrange'
		  splog, '*YOFF - Off manual yrange'
		  splog, '*T - change nsmooth'
		  splog, '*TUP - smooth up +'
		  splog, '*TDW - smooth down -'
		  splog, '*G - change sigma of Gaussian smooth'
		  splog, '*GUP - gaussian smooth up +'
		  splog, '*GDW - gaussian smooth down -'
		  splog, '***CHECKSKY - Show the positions of sky OI6301 using the object profile'
		  splog, '***SHOWLOG - Show internal variables'
		  splog, '***PRINTLIST - Print Object Name List'
		  splog, '***SHOWDLAMBDA - Show wavelength intervals'
		  splog, '***RUNDETECT - Run interactive detection '
		  splog, '***ENDDETECT - End interactive detection '
		  splog, '***DUMPDETECT - Run batch detection mode : dumping all the detections to files '
		  splog, '***DUMPDETECTNOSMOOTH - Run batch detection mode : dumping all the detections to files (no gaussian smoothing) '
		  splog, '***DUMPNODETECTNOSMOOTH - Run batch detection mode : dumping all the detections to files (no gaussian smoothing) '
		  splog, '***FIND - Find the right LAE using ID number '
		  splog, '***HELP - Show this'
       endif


   
         if strupcase(test) eq 'T' then begin
           splog, 'Enter New Top Hat Smooth Width'
           read, prompt='nsmooth : ', insmooth
		   nsmooth = insmooth
		   i = i
           loop = 1
        endif

          if strupcase(test) eq 'TUP' then begin
           splog, 'Top Hat Smooth Width Up'
		   nsmooth=nsmooth+2
		   i = i
           loop = 1 
        endif
        
           if strupcase(test) eq 'TDW' then begin
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
 

        if strupcase(test) eq 'CHECKSKY' then begin
        splog, '>>>>>>> SKY CHECK MODE:                                                 <<<<<<<<<'
        splog, '>>>>>>> 1. find the suitable "gsmooth length" before running this       <<<<<<<<<'
        splog, '>>>>>>>    ,because this routine use the gsmooth kernel of gsmooth mode <<<<<<<<<'
        splog, '>>>>>>> 2. type stop to exit       <<<<<<<<<'
        splog, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

		central_wave = 0.
		width_wave = 0.
        read, prompt='>>>> Input Central Wavelength (ex. 6563.2) : ', central_wave
        read, prompt='>>>> Input Delta Lambda (ex. 10): ', width_wave
		wavemin = central_wave - width_wave
		wavemax = central_wave + width_wave
		isLoop='yes'
		ii=i
		while strupcase(isLoop) ne 'STOP' do begin
			   	temp_obj = object[*,ii]				
				temp_obj_gsmooth = convol(temp_obj,shong_gaussian_ker2,/center)
				pix_max=0

   			   lam_min = min(abs(lam[*,ii]-wavemin),pix_lammin) ; 6290 
			   lam_max = min(abs(lam[*,ii]-wavemax),pix_lammax) ; 6310 
			   var_max = max(temp_obj_gsmooth[pix_lammin:pix_lammax],pix_max) 
			   pix_max = pix_max + pix_lammin
				temp_ymax = max(temp_obj_gsmooth[pix_lammin:pix_lammax])* 1.1
				temp_ymin = min([0.0, min(temp_obj_gsmooth[pix_lammin:pix_lammax])])
				temp_diff = temp_ymax - temp_ymin

				;djs_oplot,[lam[pix_max,i]],[unmasked_var_gsmooth[pix_max]],psym=2,thick=3,symsize=3,color='red'

			  ; oplot,[lam[pix_max,ii]],[sqrt(unmasked_var_gsmooth[pix_max])],psym=2,thick=3,symsize=3,color=100

;			   print,">>Maximum Position>> ", plugmap[ii].objtype,$
;					pix_lammin, pix_lammax,pix_max,$
;					lam[pix_lammin,ii],lam[pix_lammax,ii], lam[pix_max,ii], $
;					temp_obj_gsmooth[pix_max],shong_gaussian_kernel_sigma
				
			   print,">>Maximum Position>> ", plugmap[ii].objtype,$
					lam[pix_max,ii], temp_obj_gsmooth[pix_max],shong_gaussian_kernel_sigma


			   junky = mpfitpeak(lam[pix_lammin:pix_lammax,ii],temp_obj_gsmooth[pix_lammin:pix_lammax],outpara)
;			   print,">>Fit the Peak Position>> ", plugmap[ii].objtype,outpara(1),outpara(0)+outpara(3),outpara(2),shong_gaussian_kernel_sigma
			   print,">>Fit the Peak Position>> ", plugmap[ii].objtype,outpara(1),outpara(0)+outpara(3),outpara(2),shong_gaussian_kernel_sigma
			   print,">>Fit Gaussian: max, centroid, sigma, offset>> ", plugmap[ii].objtype,outpara
				
				fitgaussian = gaussian(lam[pix_lammin:pix_lammax,ii],outpara)
				
			;;;; PLOT
		     djs_plot, lam[*,ii], temp_obj_gsmooth,$
		        /xstyle , yrange=[temp_ymin, temp_ymax],  ps=psym, title=title1, $
		        xrange=xrange,/nstyle, /ystyle, thick=thick
		     djs_oplot, lam[pix_lammin:pix_lammax,ii], temp_obj_gsmooth[pix_lammin:pix_lammax], color='blue'
		     djs_oplot, lam[pix_lammin:pix_lammax,ii], fitgaussian, color='red'
	       ;; showing skylines
	           skylines= [4046.56,4358.34,5577.34,5889.95,5895.92,6300.64,6363.78]
	           skylabels = ['HgI(4046)','HgI(4358)','OI(5577)','NaI(5890)','NaI(5896)','[OI](6301)','[OI](6363)']
	           for j = 0, n_elements(skylines) -1 do begin
   		           djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
   	    	          color='magenta', linestyle=2
    	          djs_xyouts, (1+0.0002)*skylines(j), temp_ymax-temp_diff*(0.2+0.02*(-1)^j), $
            	     skylabels(j), charsize=2, color='green'
        	   endfor
		        djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
		           temp_ymin + (temp_ymax-temp_ymin)*0.90  ,  plugmap[ii].objtype, $
		           color='yellow', thick=thick, charsize=2, /isolatin1

				oplot,[outpara(1)],[outpara(0)+outpara(3)],psym=2,thick=3,symsize=3,color=200			
				
				ii=ii+1
				read, prompt='stop to exit; p for the previous; enter for the next : ', isLoop
				if strupcase(isLoop) eq 'P' then begin
					ii=ii-2
					if ii lt 0 then begin
						ii=0
					endif
				endif
				ii = ii mod 300
		endwhile

		read, prompt='Print the positions?(y/n): ', isLoop
			if strupcase(isLoop) eq 'Y' then begin
				ii=0
				for	ii=0,299 do begin		
			   	temp_obj = object[*,ii]				
				temp_obj_gsmooth = convol(temp_obj,shong_gaussian_ker2,/center)
				pix_max=0

   			   lam_min = min(abs(lam[*,ii]-wavemin),pix_lammin) ; 6290 
			   lam_max = min(abs(lam[*,ii]-wavemax),pix_lammax) ; 6310 
			   var_max = max(temp_obj_gsmooth[pix_lammin:pix_lammax],pix_max) 
			   pix_max = pix_max + pix_lammin
				temp_ymax = max(temp_obj_gsmooth[pix_lammin:pix_lammax])* 1.1
				temp_ymin = min([0.0, min(temp_obj_gsmooth[pix_lammin:pix_lammax])])
				temp_diff = temp_ymax - temp_ymin

				
;			   print,">>Maximum Position>> ", plugmap[ii].objtype,$
;					lam[pix_max,ii], temp_obj_gsmooth[pix_max],shong_gaussian_kernel_sigma

			   junky = mpfitpeak(lam[pix_lammin:pix_lammax,ii],temp_obj_gsmooth[pix_lammin:pix_lammax],outpara)
			   print,plugmap[ii].objtype, $
					outpara(1),outpara(0)+outpara(3),$
					outpara(2),shong_gaussian_kernel_sigma, $
					lam[pix_max,ii], temp_obj_gsmooth[pix_max]
			   endfor
			endif
        splog, '>>>>>>>GOODBYE: SKY CHECK MODE  <<<<<<<<<'

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

        if strupcase(test) eq 'LAE' then begin
           splog, 'Enter the Redshift you would like plotted'
           read, prompt='Observed wavelength of Lyman Alpha : ', z  ;; get an input as a wavelength
			z = (z - 1215.67)/1215.67 ;; change it into redshift scale
           
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

