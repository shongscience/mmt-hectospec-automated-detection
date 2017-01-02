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
;	goplot = 1 means show all steps in graphs 
;	goeach = 1 means stop at each fiber when visualizing the rescaling process 
;	resolution = 1.- or 2.0 ... means increasing the refinments for binning. higher value gives higher refinement; default=1.0 
;	spord = spline fit order; default = 3rd polynormial
;
;	mask: I do not use mask information. All mask related quantities are not reliable. It is necessary to edit it if you need it. 
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   January 2012 - Written by S Hong - NOAO
;------------------------------------------------------------------------------

Pro shong_combineExposures, reffile, subjfile, goplot=goplot,goeach=goeach, xrange=xrange, resolution=resolution,spord=spord

	if not keyword_set(goplot) then goplot=0 ;; default is not visualizing the combining processes
	if not keyword_set(goeach) then goeach=0 ;; default is not visualizing the combining processes
	if not keyword_set(xrange) then xrange=[3800,8500]
	if not keyword_set(resolution) then resolution=1.0d
	if not keyword_set(spord) then spord=3


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
	  window, 0, xsize=900, ysize=1200 ;make the plots on window 0 not on window 10 
  	!p.multi = [0,1,6]
		!p.charsize=2.0

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

;	inlam = dblarr(2,npix)
;	influx = dblarr(2,npix)
;	inivar = dblarr(2,npix)
;	inmask = dblarr(2,npix)
	inlam = dblarr(npix)
	influx = dblarr(npix)
	inivar = dblarr(npix)
	inmask = dblarr(npix)
	tmpflux = dblarr(npix)
	tmpivar = dblarr(npix)
	tmpormask = dblarr(npix)
	tmpandmask = dblarr(npix)
	dummystring=''
	otmpflux = dblarr(npix)
	otmpivar = dblarr(npix)
	otmpmask = dblarr(npix)
	
	tmplamspacing = dblarr(npix)
	spline_binsize = 0.0d
	tmpidxmin = 0L
	tmpidxmax = 0L
	for fiberidx=0,299 do begin
	
		print,">> size(ivar1[*,fiberidx]), size(otmpivar[*]) = ",size(ivar1[*,fiberidx]),",", size(otmpivar)
		print,">> size(object1[*,fiberidx]), size(otmpflux[*]) = ",size(object1[*,fiberidx]),",", size(otmpflux)
		print,">> size(mask1[*,fiberidx]), size(otmpmask[*]) = ",size(original_mask1[*,fiberidx]),",", size(otmpmask)
		inlam[*] = lam2[*,fiberidx]
		influx[*]= object2[*,fiberidx]
		inivar[*] = ivar2[*,fiberidx]
		inmask[*] = original_mask2[*,fiberidx]

		otmpflux[*] = double(object1[*,fiberidx])
		otmpivar[*] = double(ivar1[*,fiberidx])
		otmpmask[*] = double(original_mask1[*,fiberidx])
	
		tmplam = lam1[*,fiberidx]

		;; determine binsize for spline_interpolation
		;tmplamspacing[*] = 0.0d ;;; measureing binsize
		;for iL=0,299 do begin
		;	tmplamspacing[iL] = alog10(inlam[1,iL]) - alog10(inlam[0,iL])
		;endfor 
		;print,">> spacing of alog10(lamda) = ",moment(tmplamspacing)
		;print,"inloglam[1] - inloglam[0] =",(alog10(inlam[101]) - alog10(inlam[100]))
		spline_binsize=(alog10(inlam[1]) - alog10(inlam[0]))/resolution
		;;;;;;;;;

		combine1fiber,alog10(inlam),influx,inivar,finalmask= inmask,newloglam=alog10(tmplam),$ 
			 newflux=tmpflux, newivar=tmpivar, ormask=tmpormask,andmask=tmpandmask,$
			 binsz=spline_binsize, nord=spord,brkptbin=(1.2*spline_binsize*resolution),maxsep=(2.0*spline_binsize*resolution)
;			 newflux=tmpflux, newivar=tmpivar, ormask=tmpormask,andmask=tmpandmask,binsz=spline_binsize,nord=6

		;; combine and write the rescaled results
		newlam[*,fiberidx] = tmplam
		newobject[*,fiberidx] = tmpflux + otmpflux
		newivar[*,fiberidx] = tmpivar * otmpivar / (tmpivar + otmpivar)
		;newandmask[*,fiberidx] = min(tmpandmask[*],otmpmask[*])
		newandmask[*,fiberidx] = tmpandmask
		newormask[*,fiberidx] = tmpormask + otmpmask
		

		
		

		print,">> fiberidx, npix1-2, plugmap1-2, num(newlam) = ",fiberidx,",",n_elements(tmplam)," ",$
			n_elements(tmpflux),"  ,",plugmap1[fiberidx].objtype," ",plugmap1[fiberidx].objtype," ,",n_elements(newlam[*,fiberidx])


		if goplot eq 1 then begin
			junk = min(abs(newlam[*,fiberidx]-xrange(0)),tmpidxmin)
			junk = min(abs(newlam[*,fiberidx]-xrange(1)),tmpidxmax)

			djs_plot,lam1[*,fiberidx],abs(lam1[*,fiberidx] - lam2[*,fiberidx]),$
				title='Wavelength Difference: Red line (difference from rescaled lambda) should be zero',color='blue',xrange=xrange, thick=4.0
			djs_oplot,lam1[*,fiberidx],abs(lam1[*,fiberidx] - newlam[*,fiberidx]),color='red',xrange=xrange	
			djs_oplot,lam1[*,fiberidx],abs(lam2[*,fiberidx] - newlam[*,fiberidx]),color='yellow',xrange=xrange	

			djs_plot,lam2[*,fiberidx],object2[*,fiberidx],color='yellow',xrange=xrange,$
				title='Yellow = original spectrum; Red  = rescaled spectrum',$
				yrange=[min(influx[tmpidxmin:tmpidxmax]),max(influx[tmpidxmin:tmpidxmax])],thick=1.5
			;djs_oplot,lam1[*,fiberidx],object1[*,fiberidx],color='blue'
			djs_oplot,tmplam,tmpflux,color='red'
			;djs_oplot,newlam[*,fiberidx],newobject[*,fiberidx],color='blue'

			djs_plot,lam2[*,fiberidx],ivar2[*,fiberidx],color='yellow',$
				title='Yellow = original ivar; Red  = rescaled ivar',$
				xrange=xrange,yrange=[min(ivar2[*,fiberidx]),max(tmpivar)]
			;djs_oplot,lam1[*,fiberidx],ivar1[*,fiberidx],color='yellow'
			djs_oplot,tmplam,tmpivar,color='red'
			;djs_oplot,newlam[*,fiberidx],newivar[*,fiberidx],color='red'

			djs_plot,lam2[*,fiberidx],original_mask2[*,fiberidx],color='yellow',$
				title='Masks: Yellow = original, Red = ormask, Magenta = andmask; currently not working well',$
				xrange=xrange,yrange=[min(tmpormask),max(tmpormask)]
			;djs_oplot,lam1[*,fiberidx],original_mask1[*,fiberidx],color='yellow'
			djs_oplot,tmplam,tmpormask,color='red'
			djs_oplot,tmplam,tmpandmask,color='magenta'
			;djs_oplot,newlam[*,fiberidx],newandmask[*,fiberidx],color='magenta',thick=2.0
			;djs_oplot,newlam[*,fiberidx],newormask[*,fiberidx],color='red',thick=2.0

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			

			djs_plot,lam1[*,fiberidx],object1[*,fiberidx],color='blue',xrange=xrange,$
				title='Blue = original reference profile; Red  = rescaled subject profile; Magenta = combined profile',$
				yrange=[min(newobject[tmpidxmin:tmpidxmax,fiberidx]),max(newobject[tmpidxmin:tmpidxmax,fiberidx])]
			djs_oplot,tmplam,tmpflux, color='red'
			djs_oplot,newlam[*,fiberidx],newobject[*,fiberidx],color='magenta'

			djs_plot,lam1[*,fiberidx],ivar1[*,fiberidx],color='blue',$
				title='Blue = original reference ivar; Red  = rescaled subject ivar; Magenta = combined ivar',$
				xrange=xrange,yrange=[0.0,max(inivar)]
			;djs_oplot,lam1[*,fiberidx],ivar1[*,fiberidx],color='yellow'
			djs_oplot,tmplam,tmpivar, color='red'
			djs_oplot,newlam[*,fiberidx],newivar[*,fiberidx],color='magenta'


			if goeach eq 1 then begin ;; if want to stop at each object
				read, dummystring, prompt='Press anykey for the next (e for exit) :'
				if strupcase(dummystring) eq 'E' then begin 
					fiberidx=400
				endif
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



		outname='shong_combined.fits'

    mwrfits, float(lam1), outname, /create, objhdr
    mwrfits, float(newobject), outname
    mwrfits, float(newivar), outname
    mwrfits, float(newandmask), outname
    mwrfits, float(newormask), outname
    mwrfits, plugmap1, outname
    mwrfits, sky_one, outname
    mwrfits, sky_one_scl, outname
    mwrfits, sky_two, outname
    mwrfits, sky_two_scl, outname

end

