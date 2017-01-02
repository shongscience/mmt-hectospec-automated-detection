;+
; NAME:
;   shong_combine
;
; PURPOSE:
;   combine two separately reduced Hect.fits files. Modified by S. Hong. for my customized purpose. 
;
; CALLING SEQUENCE:
;  hs_combine, filename1, filename2, nstart=, plugmap=, ivar=, nsmooth=, xrange=, 
;
; INPUTS:
;  filename    - filename for for the Hectospec data
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
;   June 2012 - Written by S Hong - NOAO
;------------------------------------------------------------------------------

PRO shong_combine, filename1, filename2, nstart=nstart, plugmap=plugmap, ivar=ivar, $
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
  lam1 = mrdfits(filename1, 0)
  object1 = mrdfits(filename1, 1)
  ivar1 = mrdfits(filename1, 2)
  mask = mrdfits(filename1, 4)*0.0
  original_mask = mrdfits(filename1, 4)
  sky_one = mrdfits(filename1, 6)
  sky_two = mrdfits(filename1, 8)
  plugmap = mrdfits(filename1, 5)
  skyoneflux = bspline_valu(lam1, sky_one)

  lam2 = mrdfits(filename2, 0)
  object2 = mrdfits(filename2, 1)
  ivar2 = mrdfits(filename2, 2)

	xrange=[4000,5000]
  
  window, 0, xs=900, ys=1000 ;make the plots on window 0 not on window 10 
  !p.multi = [0,1,3]
	junk=''

	for i=0, 299 do begin
		print, i 
	     djs_plot, lam1[*,i], object1[*,i], xrange=xrange 
	     djs_oplot, lam1[*,i], 3.0D*sqrt(1.0d/ivar1[*,i]), color='red'
	     djs_plot, lam2[*,i], object2[*,i], xrange=xrange 
	     djs_oplot, lam2[*,i], 3.0D*sqrt(1.0d/ivar2[*,i]), color='red'
	     djs_plot, lam1[*,i],(lam2[*,i]-lam1[*,i]) 
		read,prompt='Enter for the next fiber.',junk
		;print,size(lam1)
		;print,size(lam2)
	endfor

END

