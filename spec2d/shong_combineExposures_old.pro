;+
; NAME:
;   shong_combineExposures.pro
;
; PURPOSE:
;   combine two separate day's exposures into a new single one. 
;
; CALLING SEQUENCE:
;  hs_page_file, refdaySpec.fits, anotherdaySpec.fits ;; fits files are outputs from Richard Cool's IDL routines 
;
; INPUTS:
;  two Spec.fits files
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

Pro shong_combineExposures, reffile, subjfile, goplot=goplot,goeach=goeach, xrange=xrange

	if not keyword_set(goplot) then goplot=0 ;; default is not visualizing the combining processes
	if not keyword_set(goeach) then goeach=0 ;; default is not visualizing the combining processes
	if not keyword_set(xrange) then xrange=[3800,8500]

  ;;Read the reference file
  lam1 = mrdfits(reffile, 0, hdr)
  object1 = mrdfits(reffile, 1)
  ivar1 = mrdfits(reffile, 2)
  original_mask1 = mrdfits(reffile, 4)
  plugmap1 = mrdfits(reffile, 5)
  
  mask = mrdfits(reffile, 4)*0.0
  sky_one = mrdfits(reffile, 6)
  sky_two = mrdfits(reffile, 8)
  sky_one_scl = mrdfits(reffile, 7)
  sky_two_scl = mrdfits(reffile, 9)
  skyoneflux = bspline_valu(lam1, sky_one)

	npix = n_elements(lam1[*,100]) ;; read the number of pixels of the 100th fiber
	print,">>NPIX at 100th = ", npix
	npix = n_elements(lam1[*,200]) ;; read the number of pixels of the 100th fiber
	print,">>NPIX at 200th = ", npix



  ;;Read the subject file
  lam2 = mrdfits(subjfile, 0)
  object2 = mrdfits(subjfile, 1)
  ivar2 = mrdfits(subjfile, 2)
  original_mask2 = mrdfits(subjfile, 4)
  plugmap2 = mrdfits(reffile, 5)
 

	;;Make the output fields 
	newlam = dblarr(npix,300)
	newobject = dblarr(npix,300)
	newivar = dblarr(npix,300)
	newormask = dblarr(npix,300)
	newandmask = dblarr(npix,300)


	fiberidx = 13L
	;; if goplot = 1 
	if goplot eq 1 then begin
	  window, 0, xsize=1000, ysize=1100 ;make the plots on window 0 not on window 10 
  	!p.multi = [0,1,4]

;		djs_plot,lam1[*,fiberidx],abs(lam1[*,fiberidx] - lam2[*,fiberidx]),color='red'
;		djs_oplot,lam1[*,fiberidx],abs(object1[*,fiberidx] - object2[*,fiberidx]),color='blue'
;		djs_oplot,lam1[*,fiberidx],abs(ivar1[*,fiberidx] - ivar2[*,fiberidx]),color='yellow'

;		djs_plot,lam1[*,fiberidx],object1[*,fiberidx],color='blue',xrange=[4000,5000]
;		djs_oplot,lam2[*,fiberidx],object2[*,fiberidx],color='yellow'

;		djs_plot,lam1[*,fiberidx],ivar1[*,fiberidx],color='blue'
;		djs_oplot,lam2[*,fiberidx],ivar2[*,fiberidx],color='yellow'

;		djs_plot,lam1[*,fiberidx],original_mask1[*,fiberidx],color='blue'
;		djs_oplot,lam2[*,fiberidx],original_mask2[*,fiberidx],color='yellow'
	endif





	;;; read each fiber : 0 - 299

	inlam = dblarr(2,npix)
	influx = dblarr(2,npix)
	inivar = dblarr(2,npix)
	inmask = dblarr(2,npix)
	tmpflux = dblarr(npix)
	tmpivar = dblarr(npix)
	tmpormask = dblarr(npix)
	tmpandmask = dblarr(npix)
	dummystring=''
	
	tmplamspacing = dblarr(npix)
	spline_binsize = 0.0d
	tmpidxmin = 0L
	tmpidxmax = 0L
	for fiberidx=0,299 do begin
	
		inlam[0,*] = lam1[*,fiberidx]
		inlam[1,*] = lam2[*,fiberidx]
		influx[0,*]= object1[*,fiberidx]
		influx[1,*]= object2[*,fiberidx]
		inivar[0,*] = ivar1[*,fiberidx]
		inivar[1,*] = ivar2[*,fiberidx]
		inmask[0,*] = original_mask1[*,fiberidx]
		inmask[1,*] = original_mask2[*,fiberidx]
	
		tmplam = lam1[*,fiberidx]

		;; determine binsize for spline_interpolation
		tmplamspacing[*] = 0.0d ;;; measureing binsize
		for iL=0,299 do begin
			tmplamspacing[iL] = alog10(inlam[1,iL]) - alog10(inlam[0,iL])
		endfor 
		print,">> spacing of alog10(lamda) = ",moment(tmplamspacing)
		print,"inloglam[1] - inloglam[0] =",(alog10(inlam[101]) - alog10(inlam[100]))
		spline_binsize=(alog10(inlam[1]) - alog10(inlam[0]))/4.0d
		;;;;;;;;;

		combine1fiber,alog10(inlam),influx,inivar,finalmask= inmask,newloglam=alog10(tmplam),$ 
			 newflux=tmpflux, newivar=tmpivar, ormask=tmpormask,andmask=tmpandmask,binsz=spline_binsize,nord=6

		newlam[*,fiberidx] = tmplam
		newobject[*,fiberidx] = tmpflux
		newivar[*,fiberidx] = tmpivar
		newandmask[*,fiberidx] = tmpandmask
		newormask[*,fiberidx] = tmpormask

		print,">> fiberidx, npix1-2, plugmap1-2, num(newlam) = ",fiberidx,",",n_elements(tmplam)," ",$
			n_elements(tmpflux),"  ,",plugmap1[fiberidx].objtype," ",plugmap1[fiberidx].objtype," ,",n_elements(newlam[*,fiberidx])


		if goplot eq 1 then begin
			junk = min(abs(newlam[*,fiberidx]-xrange(0)),tmpidxmin)
			junk = min(abs(newlam[*,fiberidx]-xrange(1)),tmpidxmax)

			djs_plot,lam1[*,fiberidx],abs(lam1[*,fiberidx] - lam2[*,fiberidx]),color='red',xrange=xrange
			djs_oplot,lam1[*,fiberidx],abs(lam1[*,fiberidx] - newlam[*,fiberidx]),color='yellow',xrange=xrange	

			djs_plot,lam2[*,fiberidx],object2[*,fiberidx],color='blue',xrange=xrange,$
				yrange=[min(tmpflux[tmpidxmin:tmpidxmax]),max(tmpflux[tmpidxmin:tmpidxmax])],thick=1.5
			djs_oplot,lam1[*,fiberidx],object1[*,fiberidx],color='yellow'
			;djs_oplot,tmplam,tmpflux,color='red',thick=2
			djs_oplot,newlam[*,fiberidx],newobject[*,fiberidx],color='red',thick=2.0

			djs_plot,lam1[*,fiberidx],ivar1[*,fiberidx],color='blue',xrange=xrange,yrange=[min(ivar1[*,fiberidx]),max(tmpivar)]
			djs_oplot,lam2[*,fiberidx],ivar2[*,fiberidx],color='yellow'
			;djs_oplot,tmplam,tmpivar,color='red',thick=2
			djs_oplot,newlam[*,fiberidx],newivar[*,fiberidx],color='red',thick=1.0

			djs_plot,lam1[*,fiberidx],original_mask1[*,fiberidx],color='blue',xrange=xrange,yrange=[min(tmpormask),max(tmpormask)]
			djs_oplot,lam2[*,fiberidx],original_mask2[*,fiberidx],color='yellow'
			;djs_oplot,tmplam,tmpmask,color='red',thick=2
			djs_oplot,newlam[*,fiberidx],newormask[*,fiberidx],color='red',thick=2.0
			djs_oplot,newlam[*,fiberidx],newandmask[*,fiberidx],color='magenta',thick=0.5

			if goeach eq 1 then begin ;; if want to stop at each object
				read, dummystring, prompt='Press anykey for the next :'
			endif
		endif

	endfor


;  lam1 = mrdfits(reffile, 0, hdr)
;  object1 = mrdfits(reffile, 1)
;  ivar1 = mrdfits(reffile, 2)
;  original_mask1 = mrdfits(reffile, 4)
;  plugmap1 = mrdfits(reffile, 5)  
;  mask = mrdfits(reffile, 4)*0.0
;  sky_one = mrdfits(reffile, 6)
;  sky_two = mrdfits(reffile, 8)
;  sky_one_scl = mrdfits(reffile, 7)
;  sky_two_scl = mrdfits(reffile, 9)
;  skyoneflux = bspline_valu(lam1, sky_one)



;		outname='shong_combined.fits'

;    mwrfits, float(lam1), outname, /create, objhdr
;    mwrfits, float(newobject), outname
;    mwrfits, float(newivar), outname
;    mwrfits, float(newandmask), outname
;    mwrfits, float(newormask), outname
;    mwrfits, plugmap1, outname
;    mwrfits, sky_one, outname
;    mwrfits, sky_one_scl, outname
;    mwrfits, sky_two, outname
;    mwrfits, sky_two_scl, outname

end

