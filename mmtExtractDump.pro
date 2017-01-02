;-----------------------------------
;+
; NAME:
;   mmtExtractDumpNoSmooth.pro
;
; PURPOSE:
;   After running "dumpdetect", we will get /detection/ directory storing all of detected dumps
;	This script will make a table having D001, D002, D003... 
;
; CALLING SEQUENCE:
;  mmtExtractDump,'lae.list'
;	"dump.summary" will be created.
; INPUTS:
;  
;
; OPTIONAL KEYWORDS:
;
;
; OUTPUTS:
;	dump.summary : print 1st, 2nd, and 3rd detections
;	dump.detection : print final detection, 1st, 2nd detections 
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS: 
;		Depends on raString.pro and decString.pro (included)
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2012 - Written by S Hong - NOAO
;------------------------------------------------------------------------------

pro mmtExtractDump, dumpcatlist, histymax=histymax, skipdump=skipdump, dldp=dldp

	
	if NOT keyword_set(histymax) then histymax=10.0d
	if NOT keyword_set(skipdump) then skipdump=0
	if NOT keyword_set(dldp) then dldp=1.0                 ;; dldp = dLambda / dpixel = pixel sampling rate ... in this case 1.1A/pix
	if keyword_set(skipdump) then skipdump=1
	

	;; read the list file
	readcol,dumpcatlist,f='a',filelist
	numlist = n_elements(filelist)
	;	print,filelist 
	;	print,n_elements(filelist)



	;;------------
	readcol,'/Users/shong/idlsrc/mypro/ib445_avg.bp',f='f,f',nflam,nfthru
	normthru = max(nfthru)
	normthru = histymax/normthru
	nfthru[*] = nfthru[*] * normthru

	;;-------------------
	;; read D001, D002, D003. 
	;;	if there is no secondary detection, we set all values "-1"
	numlist--
	numEachDetection=0L
	info = ''
	info2nd =''
	info3rd =''
	infodtcode = '' ;; D100, D110, D111 codes
	openw,dumplun, 'dump.summary',/get_lun
	for iL=0,numlist do begin 
		print,"***Reading...",filelist[iL]
		readcol,filelist[iL],f='a,x,f,f,i,i,x,x,x,f,f,f',objname,ra,dec,idxstart,idxwidth,centroid,fmax,fwhm
			print,string(format='(" **ObjName = ",a,"; detections = ",i02)',objname[0],n_elements(objname))
		numEachDetection = n_elements(objname)

		info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ")',$
  				objname[0],ra[0],dec[0],centroid[0],fmax[0],fwhm[0],idxstart[0],idxwidth[0])

		infodtcode='Dxxx'
		if numEachDetection ge 3 then begin 		
			info2nd = string(format='(3f12.3," ",i5," ",i3," ")',centroid[1],fmax[1],fwhm[1],idxstart[1],idxwidth[1])
			info3rd = string(format='(3f12.3," ",i5," ",i3," ")',centroid[2],fmax[2],fwhm[2],idxstart[2],idxwidth[2])
			infodtcode='D111'
		endif
		if numEachDetection eq 2 then begin 		
			info2nd = string(format='(3f12.3," ",i5," ",i3," ")',centroid[1],fmax[1],fwhm[1],idxstart[1],idxwidth[1])
			info3rd = string(format='(3f12.3," ",i5," ",i3," ")',-1,-1,-1,-1,-1)
			infodtcode='D110'
		endif
		if numEachDetection eq 1 then begin 		
			info2nd = string(format='(3f12.3," ",i5," ",i3," ")',-1,-1,-1,-1,-1)
			info3rd = string(format='(3f12.3," ",i5," ",i3," ")',-1,-1,-1,-1,-1)
			infodtcode='D100'
		endif
	
		print,info+info2nd+info3rd+infodtcode
		printf,dumplun,info+info2nd+info3rd+infodtcode
	endfor
	print,"***Wrinting dump.summary ..." 
	free_lun,dumplun ;; create a dumping file "dump.summary"
	;;==========================

	;;==========================
	;; Specific Section (only valid for lyman alpha detection in 4000-5000 for mmt data


	;--- make dump.summary
	readcol,'dump.summary',f='a,f,f,f,f,f,i,i,f,f,f,i,i,f,f,f,i,i,a',objname,ra,dec,fcent,ffmax,ffwhm,$
			fidx,fwd,scent,sfmax,sfwhm,sidx,swd,tcent,tfmax,tfwhm,tidx,twd,dtcode
	idxNoneDetect = where(scent lt 0.0)
	scent[idxNoneDetect] = 2905.0; for the case of no secondary detection, we set the centroid lambda 2905.0
	idxNoneDetect = where(tcent lt 0.0)
	tcent[idxNoneDetect] = 2905.0; for the case of no secondary detection, we set the centroid lambda 2905.0

	;-------- dtcode
	idxtff = where(dtcode eq 'D100')
	numd100 = n_elements(idxtff) 
	idxtff = where(dtcode eq 'D110')
	numd110 = n_elements(idxtff) 
	idxtff = where(dtcode eq 'D111')
	numd111 = n_elements(idxtff) 




	numDump = n_elements(objname)
	numDump = numDump-1

;	print,">>>>dtcode<<<<"
;	print,dtcode[0]
;	print,dtcode[1]

;	window,0,retain=2,xsize=1100,ysize=1000








	charangstrom = '('+STRING(197B)+')'

	;;==========================================================
	; plot cut of all dumps... cut function is y = x - 8 
	!p.charsize=1.5
   	!p.multi=[0,4,3]
    ps_open,'dumpcutsummaryinfo',/ps,/encap,/color
    device,/times,xsiz=11.0,ysiz=6.5,/inch
	lambdatitle = TeXtoIDL("\lambda ")+charangstrom


	;-------- plot the detection probability
	restore,'/Users/shong/work/newdata/monte/nosmoothfiterror/DFcontour.sav' ;; hard-wired
		;;;save,xdensfunc,ydensfunc,densfunc,falseDetectionOneSig,trueFlux,fwhm,dProb,filename='DFcontour.sav' 
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;-------- plot the detection probability
	restore,'/Users/shong/work/newdata/monte/nosmoothresults/monte_results_refine_fluxscale_nosmooth.sav' ;; hard-wired
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''




	;;;---------
	;; measure dp,fp ;; dp: detection prob, fp: false prob
	tmpfluxidx = 0L
	tmpfwhmidx = 0L
	numDump = n_elements(ffmax)
	fdp = dblarr(numDump)
	sdp = dblarr(numDump)
	tdp = dblarr(numDump)
	ffp = dblarr(numDump)
	sfp = dblarr(numDump)
	tfp = dblarr(numDump)
	numDump--
	for iL=0L,numDump do begin
		junk = min(abs(ffwhm[iL]/dldp - fwhm),tmpfwhmidx)	
		junk = min(abs(ffmax[iL]*ffwhm[iL]*1.06447/dldp - trueFlux),tmpfluxidx)	
		fdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		ffp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,ffwhm[iL],(ffwhm[iL]*ffmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor
	for iL=0L,numDump do begin
		junk = min(abs(sfwhm[iL]/dldp - fwhm),tmpfwhmidx)	
		junk = min(abs(sfmax[iL]*sfwhm[iL]*1.06447/dldp - trueFlux),tmpfluxidx)	
		sdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		sfp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,sfwhm[iL],(sfwhm[iL]*sfmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor
	for iL=0L,numDump do begin
		junk = min(abs(tfwhm[iL]/dldp - fwhm),tmpfwhmidx)	
		junk = min(abs(tfmax[iL]*tfwhm[iL]*1.06447/dldp - trueFlux),tmpfluxidx)	
		tdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		tfp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,tfwhm[iL],(tfwhm[iL]*tfmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor

;--- new criteria  1
	fcutidx = where(fdp gt 0.5)
	scutidx = where(sdp gt 0.5)
	tcutidx = where(tdp gt 0.5)
	fcutoidx = where(fdp le 0.5)
	scutoidx = where(sdp le 0.5)
	tcutoidx = where(tdp le 0.5)
;--- new criteria  2
;	fcutidx = where((ffp lt 0.05) and (fdp gt 0.01) )
;	scutidx = where((sfp lt 0.05) and (sdp gt 0.01) )
;	tcutidx = where((tfp lt 0.05) and (tdp gt 0.01) )
;	fcutoidx = where((ffp ge 0.05) and (fdp gt 0.01) )
;	scutoidx = where((sfp ge 0.05) and (sdp gt 0.01) )
;	tcutoidx = where((tfp ge 0.05) and (tdp gt 0.01) )
;--- old criteria
;	scutidx = where((sfmax*sfwhm*1.06447 - 8.0d) gt sfwhm)
;	tcutidx = where((tfmax*tfwhm*1.06447 - 8.0d) gt tfwhm)
;	fcutoidx = where((ffmax*ffwhm*1.06447 - 8.0d) le ffwhm)
;	scutoidx = where((sfmax*sfwhm*1.06447 - 8.0d) le sfwhm)
;	tcutoidx = where((tfmax*tfwhm*1.06447 - 8.0d) le tfwhm)
;	cutfunc = dblarr(n_elements(trueFlux))
;	cutfunc = trueFlux - 8.0d
	
	writecol,'finaldump.summary',objname,ra,dec,fcent,ffmax,ffwhm,$
			fidx,fwd,scent,sfmax,sfwhm,sidx,swd,tcent,tfmax,tfwhm,tidx,twd,$
			fdp,ffp,sdp,sfp,tdp,tfp,dtcode,$
			fmt='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",3f12.3," ",i5," ",i3," ",3f12.3," ",i5," ",i3," ",5(f6.4," "),f6.4," ",a-6)'

	print,"size filelist = ", size(filelist)
	print,"size fcutidx = ", size(fcutidx)
	print,"size fcutoidx = ", size(fcutoidx)
	print,"size scutidx = ", size(scutidx)
	print,"size scutoidx = ", size(scutoidx)
	print,"size tcutidx = ", size(tcutidx)
	print,"size tcutoidx = ", size(tcutoidx)


	; 111111111111111111111111111111111111111
	;-------- plot histogram of dump.summary 
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent,binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0


	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Noise (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[fcutoidx],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Detection (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[fcutidx],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0






	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	labelLevels = [0.5,0.95,0.99]
	labelNames = ['0.5','0.95','0.99']
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)



	stitle = string(format='("D: (T,W) = (",f4.1,"," ,f4.1,")")',(dThreshold+iThres*0.1),(dWidth+iWds))
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	djs_oplot,ffmax*ffwhm*1.06447/dldp,ffwhm/dldp,psym=7,symsize=0.3,color='red'
;	djs_oplot,trueFlux,cutfunc,color='yellow'



	;2222222222222222222222222222222222222222
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent,binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Noise (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[scutoidx],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Detection (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[scutidx],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0





	;-------- plot the detection probability

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)


;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	djs_oplot,sfmax*sfwhm*1.06447/dldp,sfwhm/dldp,psym=7,symsize=0.3,color='grey'



	;33333333333333333333333333333333333333333333
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent,binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Noise (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent[tcutoidx],binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Detection (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent[tcutidx],binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0




	;-------- plot the detection probability
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)

;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow

	djs_oplot,tfmax*tfwhm*1.06447/dldp,tfwhm/dldp,psym=7,symsize=0.3,color='black'

	ps_close 
	!p.multi=0
	;write_png,'./dumpsummaryINFO.png',tvrd(true=1)




	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	!p.charsize=1.8
   	!p.multi=[0,3,3]
    ps_open,'dumpsummaryinfo',/ps,/encap,/color
    device,/times,xsiz=8.0,ysiz=6.5,/inch
	lambdatitle = TeXtoIDL("\lambda ")+charangstrom

	; 111111111111111111111111111111111111111
	;-------- plot histogram of dump.summary 
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Histogram (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent,binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	;-------- plot the detection probability
	restore,'/Users/shong/work/newdata/monte/nosmoothresults/monte_results_refine_fluxscale_nosmooth.sav' ;; hard-wired
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''


	djs_plot,fdp,fcent,xtitle='Completeness',psym=7,symsize=0.5,color='red',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0

	djs_plot,(1.0-ffp),fcent,xtitle='Reliability',psym=7,symsize=0.5,color='red',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0




	;2222222222222222222222222222222222222222
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Histogram (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent,binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,sdp,scent,xtitle='Completeness',psym=7,symsize=0.5,color='grey',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0

	djs_plot,(1.0-sfp),scent,xtitle='Reliability',psym=7,symsize=0.5,color='grey',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0








	;33333333333333333333333333333333333333333333
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='Histogram (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent,binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,tdp,tcent,xtitle='Completeness',psym=7,symsize=0.5,color='black',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0

	djs_plot,(1.0-tfp),tcent,xtitle='Reliability',psym=7,symsize=0.5,color='black',yr=[4000,5000]
	djs_oplot,[0,100],[4047,4047],color='red',linestyle=1,thick=2.0
	djs_oplot,[0,100],[4358,4358],color='red',linestyle=1,thick=2.0


	itdp = where(tdp ge 0.5)
	print,">>>>>>>>>>>>>>>>>>>>>>>> third detection > detection prob >= 0.5 <<<<<<<<<<<<<<<<<<<<"
	print,objname[itdp]
	print,tdp[itdp]
	print,tcent[itdp]
	print,">>>>>>>>>>>>>>>>>>>>>>>> END: third detection > detection prob >= 0.5 <<<<<<<<<<<<<<<<<<<<"



	ps_close 
	!p.multi=0
	;write_png,'./dumpsummaryINFO.png',tvrd(true=1)













	;;==========================================================
	; plot cut of all dumps... cut function is y = x - 8 
	!p.charsize=1.5
   	!p.multi=[0,3,2]
    ps_open,'dumphistoinfo',/ps,/encap,/color
    device,/times,xsiz=7.5,ysiz=4.0,/inch
	lambdatitle = TeXtoIDL("\lambda ")+charangstrom


	;-------- plot the detection probability
	restore,'/Users/shong/work/newdata/monte/nosmoothfiterror/DFcontour.sav' ;; hard-wired
		;;;save,xdensfunc,ydensfunc,densfunc,falseDetectionOneSig,trueFlux,fwhm,dProb,filename='DFcontour.sav' 
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;;;---------
	;; measure dp,fp ;; dp: detection prob, fp: false prob
	tmpfluxidx = 0L
	tmpfwhmidx = 0L
	numDump = n_elements(ffmax)
	fdp = dblarr(numDump)
	sdp = dblarr(numDump)
	tdp = dblarr(numDump)
	ffp = dblarr(numDump)
	sfp = dblarr(numDump)
	tfp = dblarr(numDump)
	numDump--
	for iL=0L,numDump do begin
		junk = min(abs(ffwhm[iL]/dldp - fwhm),tmpfwhmidx)	
		junk = min(abs(ffmax[iL]*ffwhm[iL]*1.06447/dldp - trueFlux),tmpfluxidx)	
		fdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		ffp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,ffwhm[iL],(ffwhm[iL]*ffmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor
	for iL=0L,numDump do begin
		junk = min(abs(sfwhm[iL]/dldp - fwhm),tmpfwhmidx)	
		junk = min(abs(sfmax[iL]*sfwhm[iL]*1.06447/dldp - trueFlux),tmpfluxidx)	
		sdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		sfp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,sfwhm[iL],(sfwhm[iL]*sfmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor
	for iL=0L,numDump do begin
		junk = min(abs(tfwhm[iL] - fwhm),tmpfwhmidx)	
		junk = min(abs(tfmax[iL]*tfwhm[iL]*1.06447 - trueFlux),tmpfluxidx)	
		tdp[iL] = dProb[0,0,tmpfluxidx,tmpfwhmidx]
		tfp[iL]= falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
		print,tfwhm[iL],(tfwhm[iL]*tfmax[iL]*1.06447),fwhm[tmpfwhmidx],trueFlux[tmpfluxidx],$
			dProb[0,0,tmpfluxidx,tmpfwhmidx],falseDetectionOneSig[tmpfluxidx,tmpfwhmidx]
	endfor

;--- new criteria  1
	fcutidx = where(fdp gt 0.5)
	scutidx = where(sdp gt 0.5)
	tcutidx = where(tdp gt 0.5)
	fcutoidx = where(fdp le 0.5)
	scutoidx = where(sdp le 0.5)
	tcutoidx = where(tdp le 0.5)
;--- new criteria  2
;	fcutidx = where((ffp lt 0.05) and (fdp gt 0.01) )
;	scutidx = where((sfp lt 0.05) and (sdp gt 0.01) )
;	tcutidx = where((tfp lt 0.05) and (tdp gt 0.01) )
;	fcutoidx = where((ffp ge 0.05) and (fdp gt 0.01) )
;	scutoidx = where((sfp ge 0.05) and (sdp gt 0.01) )
;	tcutoidx = where((tfp ge 0.05) and (tdp gt 0.01) )
;--- old criteria
;	scutidx = where((sfmax*sfwhm*1.06447 - 8.0d) gt sfwhm)
;	tcutidx = where((tfmax*tfwhm*1.06447 - 8.0d) gt tfwhm)
;	fcutoidx = where((ffmax*ffwhm*1.06447 - 8.0d) le ffwhm)
;	scutoidx = where((sfmax*sfwhm*1.06447 - 8.0d) le sfwhm)
;	tcutoidx = where((tfmax*tfwhm*1.06447 - 8.0d) le tfwhm)
;	cutfunc = dblarr(n_elements(trueFlux))
;	cutfunc = trueFlux - 8.0d
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	labelLevels = [0.5,0.95,0.99]
	labelNames = ['0.5','0.95','0.99']
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	iThres = 0 
	iWds = 0 
	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)

	stitle = string(format='("D: (T,W) = (",f4.1,"," ,f4.1,")")',(dThreshold+iThres*0.1),(dWidth+iWds))
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	djs_oplot,ffmax*ffwhm*1.06447/dldp,ffwhm/dldp,psym=7,symsize=0.3,color='red'
;	djs_oplot,trueFlux,cutfunc,color='yellow'






		;-------- plot the detection probability

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)


;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	djs_oplot,sfmax*sfwhm*1.06447/dldp,sfwhm/dldp,psym=7,symsize=0.3,color='grey'


	;-------- plot the detection probability
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)

;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow

	djs_oplot,tfmax*tfwhm*1.06447/dldp,tfwhm/dldp,psym=7,symsize=0.3,color='black'


	; 111111111111111111111111111111111111111
	;-------- plot histogram of dump.summary 
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent,binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	;2222222222222222222222222222222222222222
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent,binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0


	;33333333333333333333333333333333333333333333
	;-------- plot histogram of dump.summary
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='All (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent,binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0



	ps_close 
	!p.multi=0
	;write_png,'./dumpsummaryINFO.png',tvrd(true=1)












	;;==========================================================

	;;;;;;;;;;;;;;
	; make a plot for 1st detection with three category 
	; 1. only 1 detection , 2. has 2nd, 3. has third

	charangstrom = '('+STRING(197B)+')'
	!p.charsize=1.4
   	!p.multi=[0,4,3]
    ps_open,'firstdetectsummaryinfo',/ps,/encap,/color
    device,/times,xsiz=8.0,ysiz=6.0,/inch
	lambdatitle = TeXtoIDL("\lambda ")+charangstrom

	idxttt = where((scent gt 3000.0d) and (tcent gt 3000.0d))
	idxttf = where((scent gt 3000.0d) and (tcent lt 3000.0d))
	idxtff = where(scent lt 3000.0d)

	;-- old criteria
;	cfidxtff = where((ffmax[idxtff]*ffwhm[idxtff]*1.06447 - 8.0d) le ffwhm[idxtff])
;	cfidxttf = where((ffmax[idxttf]*ffwhm[idxttf]*1.06447 - 8.0d) le ffwhm[idxttf])
;	cfidxttt = where((ffmax[idxttt]*ffwhm[idxttt]*1.06447 - 8.0d) le ffwhm[idxttt])
;	csidxttf = where((sfmax[idxttf]*sfwhm[idxttf]*1.06447 - 8.0d) le sfwhm[idxttf])
;	csidxttt = where((sfmax[idxttt]*sfwhm[idxttt]*1.06447 - 8.0d) le sfwhm[idxttt])
;	ctidxttt = where((tfmax[idxttt]*tfwhm[idxttt]*1.06447 - 8.0d) le tfwhm[idxttt])

	cfidxtff = where(fdp[idxtff] le 0.5)
	cfidxttf = where(fdp[idxttf] le 0.5)
	cfidxttt = where(fdp[idxttt] le 0.5)
	csidxttf = where(sdp[idxttf] le 0.5)
	csidxttt = where(sdp[idxttt] le 0.5)
	ctidxttt = where(sdp[idxttt] le 0.5)

	cfidxttt = idxttt[cfidxttt]
	csidxttt = idxttt[csidxttt]
	ctidxttt = idxttt[ctidxttt]
	cfidxttf = idxttf[cfidxttf]
	csidxttf = idxttf[csidxttf]
	cfidxtff = idxtff[cfidxtff]


	print,"ttt f, ttt cf = ",n_elements(idxttt)," ",n_elements(cfidxttt)
	print,"ttt s, ttt cs = ",n_elements(idxttt)," ",n_elements(csidxttt)
	print,"ttt t, ttt ct = ",n_elements(idxttt)," ",n_elements(ctidxttt)
	print,"ttf f, ttf cf = ",n_elements(idxttf)," ",n_elements(cfidxttf)
	print,"ttf s, ttf sf = ",n_elements(idxttf)," ",n_elements(csidxttf)
	print,"tff f, tff cf = ",n_elements(idxtff)," ",n_elements(cfidxtff)
	; 111111111111111111111111111111111111111
	;-------- plot histogram 1,1,1
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxttt],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[cfidxttt],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0



	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[idxttt],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[csidxttt],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0



	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent[idxttt],binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 1 (3rd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent[ctidxttt],binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0




	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 0 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxttf],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 0 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[cfidxttf],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0






	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 0 (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[idxttf],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0


	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 1 0 (2nd detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[csidxttf],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0





	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 0 0 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxtff],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='1 0 0 (1st detection)',linestyle=1,thick=2.0
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[cfidxtff],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

;	print,tcent[idxtff]
;	print,scent[idxtff]
;	print,tcent[idxttf]
;	print,scent[idxttf]

    ps_close 
    !p.multi=0







	;;==========================================================

	;;;;;;;;;;;;;;
	; make a plot for 1st detection with three category 
	; 1. only 1 detection , 2. has 2nd, 3. has third

	charangstrom = '('+STRING(197B)+')'
	!p.charsize=1.4
   	!p.multi=[0,3,3]
    ps_open,'firsttrimsummaryinfo',/ps,/encap,/color
    device,/times,xsiz=7.0,ysiz=5.6,/inch
	lambdatitle = TeXtoIDL("\lambda ")+charangstrom

	idxttt = where((scent gt 3000.0d) and (tcent gt 3000.0d))
	idxttf = where((scent gt 3000.0d) and (tcent lt 3000.0d))
	idxtff = where(scent lt 3000.0d)

	;-- old criteria
;	cfidxtff = where((ffmax[idxtff]*ffwhm[idxtff]*1.06447 - 8.0d) le ffwhm[idxtff])
;	cfidxttf = where((ffmax[idxttf]*ffwhm[idxttf]*1.06447 - 8.0d) le ffwhm[idxttf])
;	cfidxttt = where((ffmax[idxttt]*ffwhm[idxttt]*1.06447 - 8.0d) le ffwhm[idxttt])
;	csidxttf = where((sfmax[idxttf]*sfwhm[idxttf]*1.06447 - 8.0d) le sfwhm[idxttf])
;	csidxttt = where((sfmax[idxttt]*sfwhm[idxttt]*1.06447 - 8.0d) le sfwhm[idxttt])
;	ctidxttt = where((tfmax[idxttt]*tfwhm[idxttt]*1.06447 - 8.0d) le tfwhm[idxttt])

	cfidxtff = where(fdp[idxtff] le 0.5)
	cfidxttf = where(fdp[idxttf] le 0.5)
	cfidxttt = where(fdp[idxttt] le 0.5)
	csidxttf = where(sdp[idxttf] le 0.5)
	csidxttt = where(sdp[idxttt] le 0.5)
	ctidxttt = where(sdp[idxttt] le 0.5)

	cfidxttt = idxttt[cfidxttt]
	csidxttt = idxttt[csidxttt]
	ctidxttt = idxttt[ctidxttt]
	cfidxttf = idxttf[cfidxttf]
	csidxttf = idxttf[csidxttf]
	cfidxtff = idxtff[cfidxtff]


	print,"ttt f, ttt cf = ",n_elements(idxttt)," ",n_elements(cfidxttt)
	print,"ttt s, ttt cs = ",n_elements(idxttt)," ",n_elements(csidxttt)
	print,"ttt t, ttt ct = ",n_elements(idxttt)," ",n_elements(ctidxttt)
	print,"ttf f, ttf cf = ",n_elements(idxttf)," ",n_elements(cfidxttf)
	print,"ttf s, ttf sf = ",n_elements(idxttf)," ",n_elements(csidxttf)
	print,"tff f, tff cf = ",n_elements(idxtff)," ",n_elements(cfidxtff)
	; 111111111111111111111111111111111111111
	;-------- plot histogram 1,1,1
	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D111 (1st detection)',linestyle=1,thick=2.0,position=[0.05,0.05,0.3,0.3]
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4165,4165],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4168,4168],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4420,4420],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4665,4665],[0,100],color='red',linestyle=1,thick=2.0
;	djs_oplot,[4669,4669],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxttt],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0



	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D111 (2nd detection)',linestyle=1,thick=2.0,position=[0.383,0.05,0.633,0.3] 
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[idxttt],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0


	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D111 (3rd detection)',linestyle=1,thick=2.0,position=[0.713,0.05,0.963,0.3]
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,tcent[idxttt],binsize=5.0,/fill, polycolor='black',/oplot,datacolorname='black'
	djs_oplot,nflam,nfthru,color='green',thick=2.0



	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D110 (1st detection)',linestyle=1,thick=2.0,position=[0.05,0.383,0.3,0.633]
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxttf],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0


	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D110 (2nd detection)',linestyle=1,thick=2.0,position=[0.383,0.383,0.633,0.633]
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,scent[idxttf],binsize=5.0,/fill, polycolor='grey',/oplot,datacolorname='grey'
	djs_oplot,nflam,nfthru,color='green',thick=2.0




	djs_plot,[4358,4358],[0,100],color='red',yrange=[0,histymax],xrange=[4000,5000],$
        xtitle=lambdatitle, ytitle='Counts',title='D100 (1st detection)',linestyle=1,thick=2.0,position=[0.05,0.716,0.3,0.966]
	djs_oplot,[4047,4047],[0,100],color='red',linestyle=1,thick=2.0
	cgHistoplot,fcent[idxtff],binsize=5.0,/fill, polycolor='red',/oplot,datacolorname='red'
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	djs_xyouts,11500,12500,string(format='("# of D100 = ",i5)',numd100),charsize=1.1,/device
	djs_xyouts,11500,12000,string(format='("# of D110 = ",i5)',numd110),charsize=1.1,/device
	djs_xyouts,11500,11500,string(format='("# of D111 = ",i5)',numd111),charsize=1.1,/device

;	print,tcent[idxtff]
;	print,scent[idxtff]
;	print,tcent[idxttf]
;	print,scent[idxttf]



    ps_close 
    !p.multi=0
























;	interrupt dummy
;	dummystring =''
;	read,dummystring , prompt='e for exit : '
;	if strupcase(dummystring) eq 'E' then stop 







	;; postdetection criteria: determine the final detection among D001, D002, D003
	dSL = 5.0 ;; skyline +- deltaSkyLine (dSL) : defines the redshift screening window
	info='empty'
	tailinfo=''
	openw,dumplun, 'dump.detection',/get_lun

	tmpdp = 0.0d
	tmpfp = 0.0d

	for iL=0L, numDump do begin  ;;;;; post detection criteria
		info='empty'



		boolallreject = 0

		;----- select the first detection, not within sky lines
		if fcent[iL] le (4047.0 - dSL) then begin 
			info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
  				objname[iL],ra[iL],dec[iL],fcent[iL],ffmax[iL],ffwhm[iL],fidx[iL],fwd[iL],fdp[iL],ffp[iL])
			tmpdp=fdp[iL]
			tmpfp=ffp[iL]
		endif
		if (fcent[iL] ge (4047.0 + dSL)) && (fcent[iL] le (4358.0 - dSL)) then begin 
			info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
  				objname[iL],ra[iL],dec[iL],fcent[iL],ffmax[iL],ffwhm[iL],fidx[iL],fwd[iL],fdp[iL],ffp[iL])
			tmpdp=fdp[iL]
			tmpfp=ffp[iL]
		endif
		if fcent[iL] ge (4358.0 + dSL) then begin 
			info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
  				objname[iL],ra[iL],dec[iL],fcent[iL],ffmax[iL],ffwhm[iL],fidx[iL],fwd[iL],fdp[iL],ffp[iL])
			tmpdp=fdp[iL]
			tmpfp=ffp[iL]
		endif

		;----- select the second detection, when within sky lines
		if (fcent[iL] lt (4047.0 + dSL)) && (fcent[iL] gt (4047.0 - dSL)) $    ;;in the 4047 sky line
			&&( $
				(scent[iL] le (4047.0 - dSL)) || $
				( (scent[iL] ge (4047.0 + dSL)) && (scent[iL] le (4358.0 - dSL)) ) || $
				(scent[iL] ge (4358.0 + dSL)) $
			) $
			then begin 
			info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
  				objname[iL],ra[iL],dec[iL],scent[iL],sfmax[iL],sfwhm[iL],sidx[iL],swd[iL],sdp[iL],sfp[iL])
			tmpdp=sdp[iL]
			tmpfp=sfp[iL]
		endif	
	
		if (fcent[iL] lt (4358.0 + dSL)) && (fcent[iL] gt (4358.0 - dSL)) $  ;; in the 4358 sky line 
			&&(  $
				(scent[iL] le (4047.0 - dSL)) || $
				( (scent[iL] ge (4047.0 + dSL)) && (scent[iL] le (4358.0 - dSL)) ) || $
				(scent[iL] ge (4358.0 + dSL)) $
			)$
			then begin 
			info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
  				objname[iL],ra[iL],dec[iL],scent[iL],sfmax[iL],sfwhm[iL],sidx[iL],swd[iL],sdp[iL],sfp[iL])
			tmpdp=sdp[iL]
			tmpfp=sfp[iL]
		endif	




		tailinfo = string(format='(" ",3f12.3," ",i5," ",i3," ",3f12.3," ",i5," ",i3," ",5(f6.4," "),f6.4)',$
  				fcent[iL],ffmax[iL],ffwhm[iL],fidx[iL],fwd[iL],$
				scent[iL],sfmax[iL],sfwhm[iL],sidx[iL],swd[iL],$
				fdp[iL],ffp[iL],sdp[iL],sfp[iL],tdp[iL],tfp[iL])





		;----- select the third detection, when 1st and 2nd are rejected 

		if ((abs(fcent[iL] - 4047.0) lt dSL) && (abs(scent[iL] - 4358.0) lt dSL)) then boolallreject = 1
		if ((abs(scent[iL] - 4047.0) lt dSL) && (abs(fcent[iL] - 4358.0) lt dSL)) then boolallreject = 1

		if boolallreject eq 1 then begin
	       	info = string(format='(a-20," ",2f12.7," ",3f12.3," ",i5," ",i3," ",f6.4," ",f6.4)',$
                objname[iL],ra[iL],dec[iL],tcent[iL],tfmax[iL],tfwhm[iL],tidx[iL],twd[iL],tdp[iL],tfp[iL])
	       	tmpdp=tdp[iL]
   		   	tmpfp=tfp[iL]
			print,">> First and second are rejected : ",info+tailinfo
		endif

		





		if info ne 'empty' then begin  ;; this is the post-detection criteria
			if tmpdp gt 0.5 then begin  ;; detection criteria
				if boolallreject eq 1 then print,">>>>>>>>>>>>>>>>>>>>>>>> Take the third detection > detection prob >= 0.5 and sky residuals  <<<<<<<<<<<<<<<<<<<<"
				print,info+tailinfo
				if boolallreject eq 1 then print,">>>>>>>>>>>>>>>>>>>>>>>> END : Take the third detection > detection prob >= 0.5 and sky residuals  <<<<<<<<<<<<<<<<<<<<"
				printf,dumplun,info+tailinfo
			endif else begin 
				print,"--No Detection: ",objname[iL]," ",info
			endelse
		endif else begin
			print,"--No Detection: ",objname[iL]," ",info
		endelse
	endfor
	free_lun, dumplun









;*********
if skipdump eq 0 then begin   ;; skipdump not to produce these eps files 

	;------- plot profiles of detections to eps

;window,1,retain=2,xsize=900,ysize=1000
!p.charsize=0.7
!p.multi=[0,4,8]

readcol,'dump.detection',f='a,x,x,f,f,f,i,i,f,f',name,gcentroid,gmax,gsigma,idxStart,ddidxWidth,ddprob,dfprob
;readcol,'LAE1p5First.cat',f='a,x,i,i,x,x,x,f,f,f',name,idxStart,idxWidth,gcentroid,gmax,gsigma

filename = name+'Profile.sav'
maxLine = n_elements(name)
maxLine--
idxLower=0L
idxUpper=0L
idxCenter=0L 
iL=0L

idxgLower=0L
idxgUpper=0L

;print,size(name),size(gcentroid),size(gmax),size(gsigma),size(idxStart),size(ddidxWdith)



epsFileName='dumpRealProfile'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
probinfo=' '
for iL = 0L,maxLine do begin
	tmpIndex = iL / 32
	
	;print,"# (",iL,"/",maxLine,"); epsIndex tmpIndex : ",epsIndex, tmpIndex

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

	
    restore,filename[iL]
	;print,filename[iL]," *** iL = ",iL," size idxWidth = ",size(ddidxWidth)
	;print,"idxStart[iL] ddidxWidth[iL] = ",idxStart[iL]," ",ddidxWidth[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
	;print,filename[iL]," idxCenter = ", idxCenter
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower lt 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4

	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor





epsFileName='dumpPseudoProfile'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
for iL = 0L,maxLine do begin
	tmpIndex = iL / 32
	
	;print,"# ",iL," eps tmp : ",epsIndex, tmpIndex

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

    restore,filename[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower le 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
	minY = min(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
;	maxY = max(originalSignalBaseline[idxLower:idxUpper])
;	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", ",f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytile='Pseudo-profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1, xtics=1,yticks=3, $
			ytickv=[minY,0,maxYY,maxY],color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='black',thick=2

;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
;			title=strupcase(name[iL]),xmargin=[6,2], xstyle=1, xtics=1,yticks=3, color='grey',charsize=0.7, yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.5, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4

	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor
;print, filename, idxStart, idxWidth
;print, maxLine
ps_close







!p.multi=[0,3,7]
epsFileName='dumpBothProfile'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
for iL = 0L,maxLine do begin
	tmpIndex = iL / 7
	
	;print,"# ",iL," eps tmp : ",epsIndex, tmpIndex

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

    restore,filename[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower le 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	;------------------ first
	maxY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
	minY = min(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
;	maxY = max(originalSignalBaseline[idxLower:idxUpper])
;	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", ",f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytile='Pseudo-profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1, xtics=1,yticks=3, $
			ytickv=[minY,0,maxYY,maxY],color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='black',thick=2

;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
;			title=strupcase(name[iL]),xmargin=[6,2], xstyle=1, xtics=1,yticks=3, color='grey',charsize=0.7, yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.5, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4



	;---------------- second figure
	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',ddprob[iL],dtmp)
	djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4


	;---------------- third figure
	djs_plot,lam,originalSignalBaseline,xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7,xrange=[4000,5000], yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,4050,(0.7*maxY),probinfo,charsize=0.4




	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor
;print, filename, idxStart, idxWidth
;print, maxLine
ps_close
!p.multi=0




;---------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------
; plot all dumps with D100, D110, and D111 
;

readcol,'finaldump.summary',f='a,x,x,f,f,f,i,i,f,x,x,x,x,f,x,x,x,x,f,f,x,x,x,x,a',name,gcentroid,gmax,gsigma,idxStart,ddidxWidth,scent,tcent,ddprob,dfprob,dtcode
;readcol,'dump.detection',f='a,x,x,f,f,f,i,i,f,f',name,gcentroid,gmax,gsigma,idxStart,ddidxWidth,ddprob,dfprob
;readcol,'LAE1p5First.cat',f='a,x,i,i,x,x,x,f,f,f',name,idxStart,idxWidth,gcentroid,gmax,gsigma

filename = name+'Profile.sav'
idxtff = where(dtcode eq 'D100')
;maxLine = n_elements(name)
maxLine = n_elements(idxtff)
maxLine--
idxLower=0L
idxUpper=0L
idxCenter=0L 
iL=0L

idxgLower=0L
idxgUpper=0L

;print,size(name),size(gcentroid),size(gmax),size(gsigma),size(idxStart),size(ddidxWdith)

!p.multi=[0,3,7]
epsFileName='dumpAllBothProfileD100N'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
probinfo=' '

for jL = 0L,maxLine do begin
	iL = idxtff[jL]
	tmpIndex = jL / 7
	
	;print,"# ",iL," eps tmp : ",epsIndex, tmpIndex
	print,"Detection code : ",dtcode[idxtff[jL]]

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

    restore,filename[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower le 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	;------------------ first
	maxY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
	minY = min(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
;	maxY = max(originalSignalBaseline[idxLower:idxUpper])
;	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", ",f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytile='Pseudo-profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1, xtics=1,yticks=3, $
			ytickv=[minY,0,maxYY,maxY],color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='black',thick=2

;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
;			title=strupcase(name[iL]),xmargin=[6,2], xstyle=1, xtics=1,yticks=3, color='grey',charsize=0.7, yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.5, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4



	;---------------- second figure
	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',ddprob[iL],dtmp)
	djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4


	;---------------- third figure
	djs_plot,lam,originalSignalBaseline,xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7,xrange=[4000,5000], yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,4050,(0.7*maxY),probinfo,charsize=0.4




	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor
;print, filename, idxStart, idxWidth
;print, maxLine
ps_close
!p.multi=0

;---------------------------------------------------------------------------------------
filename = name+'Profile.sav'
idxtff = where(dtcode eq 'D110')
;maxLine = n_elements(name)
maxLine = n_elements(idxtff)
maxLine--
idxLower=0L
idxUpper=0L
idxCenter=0L 
iL=0L

idxgLower=0L
idxgUpper=0L

;print,size(name),size(gcentroid),size(gmax),size(gsigma),size(idxStart),size(ddidxWdith)

!p.multi=[0,3,7]
epsFileName='dumpAllBothProfileD110N'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
probinfo=' '

for jL = 0L,maxLine do begin
	iL = idxtff[jL]
	tmpIndex = jL / 7
	
	;print,"# ",iL," eps tmp : ",epsIndex, tmpIndex
	print,"Detection code : ",dtcode[idxtff[jL]]

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

    restore,filename[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower le 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	;------------------ first
	maxY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
	minY = min(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
;	maxY = max(originalSignalBaseline[idxLower:idxUpper])
;	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", ",f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytile='Pseudo-profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1, xtics=1,yticks=3, $
			ytickv=[minY,0,maxYY,maxY],color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='black',thick=2

;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
;			title=strupcase(name[iL]),xmargin=[6,2], xstyle=1, xtics=1,yticks=3, color='grey',charsize=0.7, yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.5, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4



	;---------------- second figure
	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',ddprob[iL],dtmp)
	djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4


	;---------------- third figure
	djs_plot,lam,originalSignalBaseline,xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7,xrange=[4000,5000], yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,4050,(0.7*maxY),probinfo,charsize=0.4
	djs_oplot,[scent[iL],scent[iL]],[-1000,1000],linestyle=1,color='grey', thick=5




	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor
;print, filename, idxStart, idxWidth
;print, maxLine
ps_close
!p.multi=0


;---------------------------------------------------------------------------------------
filename = name+'Profile.sav'
idxtff = where(dtcode eq 'D111')
;maxLine = n_elements(name)
maxLine = n_elements(idxtff)
maxLine--
idxLower=0L
idxUpper=0L
idxCenter=0L 
iL=0L

idxgLower=0L
idxgUpper=0L

;print,size(name),size(gcentroid),size(gmax),size(gsigma),size(idxStart),size(ddidxWdith)

!p.multi=[0,3,7]
epsFileName='dumpAllBothProfileD111N'
tmpName=''

epsIndex=-1L
tmpIndex=0L
isOpenNewEps = 0L
probinfo=' '

for jL = 0L,maxLine do begin
	iL = idxtff[jL]
	tmpIndex = jL / 7
	
	;print,"# ",iL," eps tmp : ",epsIndex, tmpIndex
	print,"Detection code : ",dtcode[idxtff[jL]]

	if tmpIndex ne epsIndex	then begin 
		isOpenNewEps = 1
	endif else begin
		isOpenNewEps = 0
	endelse

	if isOpenNewEps eq 1 then begin 
		ps_close
		print, "Closing the file : ",tmpName

		epsIndex = tmpIndex
		tmpName=epsFilename+string(format='(i03)',(epsIndex+1))
		ps_open,tmpName,/ps,/encap,/color
		device,/times,xsiz=5.5,ysiz=6.5,/inch
		print, "Opening the file : ",tmpName
	endif

    restore,filename[iL]
    idxCenter = long(idxStart[iL] + 0.5 * ddidxWidth[iL])
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower le 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	;------------------ first
	maxY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
	minY = min(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])
;	maxY = max(originalSignalBaseline[idxLower:idxUpper])
;	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(normalizedSignalSmooth[(idxLower+10):(idxUpper-10)])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", ",f5.3,")")',ddprob[iL],dtmp)
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytile='Pseudo-profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1, xtics=1,yticks=3, $
			ytickv=[minY,0,maxYY,maxY],color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='black',thick=2

;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
;			title=strupcase(name[iL]),xmargin=[6,2], xstyle=1, xtics=1,yticks=3, color='grey',charsize=0.7, yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.5, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4



	;---------------- second figure
	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]-3.0*gsigma[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (gcentroid[iL]+3.0*gsigma[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 

	dtmp=(1.0d -dfprob[iL])
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',ddprob[iL],dtmp)
	djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7, yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,lam[idxLower+5],(0.7*maxY),probinfo,charsize=0.4


	;---------------- third figure
	djs_plot,lam,originalSignalBaseline,xtitle=lambdatitle,ytile='Original Profile',$
			title=strupcase(strjoin(strsplit(name[iL],'_',/extract))),xmargin=[6,2], xstyle=1, ystyle=1,xtics=1,yticks=3,$
			ytickv=[minY,0,maxYY,maxY], color='grey',charsize=0.7,xrange=[4000,5000], yrange=[minY, maxY]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_xyouts,4050,(0.7*maxY),probinfo,charsize=0.4
	djs_oplot,[scent[iL],scent[iL]],[-1000,1000],linestyle=1,color='grey', thick=5
	djs_oplot,[tcent[iL],tcent[iL]],[-1000,1000],linestyle=1,color='grey', thick=2




	print,string(format='("**Name: ",a-20,"#",i5,"/",i5," appended to ",a-20)',filename[iL],iL,maxLine,tmpName)

endfor
;print, filename, idxStart, idxWidth
;print, maxLine
ps_close
!p.multi=0










;---------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------













;================================================================================================
; write each LAE to png file with all information
;


	

	;=============================== Read data
	readcol,'finaldump.summary',f='a,f,f,f,f,f,i,i,f,f,f,i,i,f,f,f,x,x,f,f',objname,ra,dec,fcent,ffmax,ffwhm,$
			fidx,fwd,scent,sfmax,sfwhm,sidx,swd,tcent,tfmax,tfwhm,afdp,affp
	readcol,'dump.detection',f='a,f,f,f,f,f,i,i,f,f,f,f,f,x,x,f,f,f',dobjname,dra,ddec,dfcent,dffmax,dffwhm,$
			dfidx,dfwd,dfdp,dffp,d1cent,d1fmax,d1fwhm,d2cent,d2fmax,d2fwhm
	infilename = objname+'Profile.sav'





	;== write data : directory setting
	if direxist('dumpsummarypictures') eq 0L then spawn, 'mkdir dumpsummarypictures'

	;-- for finaldump.summary
	idxNoneDetect = where(scent lt 3000.0)
	scent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
	idxNoneDetect = where(tcent lt 3000.0)
	tcent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	idxNoneDetect = where(d2cent lt 3000.0)
	d2cent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	;-- for dump.summary
;	idxNoneDetect = where(scent lt 0.0)
;	scent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
;	idxNoneDetect = where(tcent lt 0.0)
;	tcent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	numDump = n_elements(infilename)
	numDump = numDump-1

	iL=0L



	restore,'/Users/shong/work/newdata/monte/nosmoothfiterror/DFcontour.sav' ;; hard-wired
;	restore,'/Users/shong/work/newdata/monte/nosmoothresults/monte_results_refine_fluxscale_nosmooth.sav' ;; hard-wired
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;***********
	;window,0,retain=2,xsize=900,ysize=1200
	charangstrom = '('+STRING(197B)+')'
	!p.charsize=1.3
   	!p.multi=[0,2,3]

	for iL=0L,numDump do begin
	
	;===== prepare for eps output
	ps_open,'./dumpsummarypictures/'+objname[iL]+'INFO',/ps,/encap,/color
    device,/times,xsiz=5.5,ysiz=6.5,/inch

	;write_png,'./dumpsummarypictures/'+objname[iL]+'INFO.png',tvrd(true=1)

	;=============================== first figure
	djs_plot,fcent,scent,psym=7,yr=[3820,5000],xr=[4000,5000],xstyle=1,ystyle=1,symsize=0.5,$
		xtitle=TexToIDL("\lambda : Primary Detection"),ytitle=TexToIDL("\lambda : Secondary Detection"), title='1st & 2nd detections'
;	djs_oplot,[4255,4255],[-1000,10000],color='blue'
;	djs_oplot,[4620,4620],[-1000,10000],color='blue'
	djs_oplot,[4358.0-dSL,4358.0-dSL],[-1000,10000],color='blue' ;; 4358 HgI line
	djs_oplot,[4358.0+dSL,4358.0+dSL],[-1000,10000],color='blue' ;; 4358 HgI line
	djs_oplot,[4047.0-dSL,4047.0-dSL],[-1000,10000],color='blue' ;; 4047 HgI line
	djs_oplot,[4047.0+dSL,4047.0+dSL],[-1000,10000],color='blue' ;; 4047 HgI line

		;;horizontal lines
	djs_oplot,[-1000,10000],[4047.0,4047.0],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[-1000,10000],[4358.0,4358.0],color='blue',linestyle=2 ;; 4358 HgI line
    djs_oplot,[-1000,10000],[4000.0,4000.0],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[4000.0,4000.0],[-1000,4000],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[5000.0,5000.0],[-1000,4000],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[-1000,10000],[3820.0,3820.0],color='black',linestyle=0,thick=4 ;; d100 line




	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4358.0-dSL,4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL],color='red',thick=2
	djs_oplot,[4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL,4047.0-dSL],[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],color='red',thick=2
	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],color='red',thick=2
	djs_oplot,[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],color='red',thick=2
	;;;;;; overplot the current object
;    if scent[iL] gt 3950 then begin 
        djs_oplot,[fcent[iL]],[scent[iL]],psym=sym(13),color='blue',symsize=2.0,thick=3
;    endif else begin 
;        djs_oplot,[fcent[iL]],[3950],psym=sym(3),color='blue',symsize=2.0,thick=3
;    endelse


;	djs_oplot,[-1000,10000],[4620,4620],color='blue'
;	djs_oplot,[-1000,10000],[4358.0-dSL,4358.0-dSL],color='red' ;; 4358 HgI line
;	djs_oplot,[-1000,10000],[4358.0+dSL,4358.0+dSL],color='red' ;; 4358 HgI line
;	djs_oplot,[-1000,10000],[4047.0-dSL,4047.0-dSL],color='red' ;; 4047 HgI line
;	djs_oplot,[-1000,10000],[4047.0+dSL,4047.0+dSL],color='red' ;; 4047 HgI line

	;third detection
;	djs_plot,fcent,tcent,psym=7,xr=[4000,5000],yr=[3900,5000],xstyle=1,ystyle=1
;	djs_oplot,[4255,4255],[-1000,10000],color='blue'
;	djs_oplot,[4620,4620],[-1000,10000],color='blue'
;	djs_oplot,[4358.0-dSL,4358.0-dSL],[-1000,10000],color='red' ;; 4358 HgI line
;	djs_oplot,[4358.0+dSL,4358.0+dSL],[-1000,10000],color='red' ;; 4358 HgI line
;	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4255,4255,4620,4620,4255],color='blue',thick=2
;	djs_oplot,[4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL,4047.0-dSL],[4255,4255,4620,4620,4255],color='blue',thick=2

	;============================= second figure: plot histogram of dump.detect

	lambdatitle = TeXtoIDL("\lambda ")+charangstrom

	cgHistoplot,fcent,binsize=5.0,/fill, polycolor='magenta',datacolorname='magenta',$
        xtitle=lambdatitle, ytitle='Counts',yrange=[0,histymax],xrange=[4000,5000], title='Detection Histogram',xstyle=1
	djs_oplot,[fcent[iL],fcent[iL]],[-1,100],color='navy',linestyle=3,thick=2.0
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	;============================ third figure: plot the detection probability

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	labelLevels = [0.95,0.99]
	labelNames = ['0.95','0.99']
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)


	fratemoments = moment(fprob[iThres,iWds,*,*])
	fratestr = string(format='(a4,f6.4,a2,f9.6)',"F = ",fratemoments[0],string(177b),sqrt(fratemoments[1]))
	stitle = string(format='("D: (T,W) = (",f4.1,"," ,f4.1,")")',(dThreshold+iThres*0.1),(dWidth+iWds))
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
;    color=cgcolor('blue'), c_labels=replicate(1,2), $
;    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
;    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
;		xtitle=xtstr,levels=[0.01,0.05,0.5], $
;        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
;        c_labels=replicate(0,3), $
;        c_linestyle=0,c_charsize=0.5,/follow

    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=[0.5,0.95,0.99], $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	


	djs_oplot,ffmax*ffwhm*1.06447/dldp,ffwhm/dldp,psym=7,symsize=0.5,color='magenta'


	;;;
	tmpflux = ffmax[iL]*ffwhm[iL]*1.06447/dldp
	tmpfwhm =ffwhm[iL]/dldp
	if tmpflux gt 80.0d then tmpflux=80
	if tmpfwhm gt 15.0d then tmpfwhm=15.0
	;print,tmpflux,tmpfwhm
	djs_oplot,[tmpflux],[tmpfwhm],psym=sym(13),symsize=2.0,color='blue',thick=3.0

	;========================= fourth figure: plot RA-DEC position

	djs_plot,ra,dec,color='magenta',psym=1,xr=[max(dra),min(dra)],yr=[min(ddec),max(ddec)],$
			xtitle='RA', ytitle='DEC', title='Spatial Distribution',symsize=1.0
	djs_oplot,[ra[iL]],[dec[iL]],psym=sym(13),symsize=2.0,color='blue',thick=3.0


	;========================= fifth figure
    restore,infilename[iL]
	;print,infilename[iL]," *** iL = ",iL," size idxWidth = ",size(ddidxWidth)
	;print,"idxStart[iL] ddidxWidth[iL] = ",idxStart[iL]," ",ddidxWidth[iL]
    idxCenter = long(fidx[iL] + 0.5 * fwd[iL])
	;print,infilename[iL]," idxCenter = ", idxCenter
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower lt 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY

	maxYN = max(normalizedSignalSmooth[idxLower:idxUpper])
	minYN = min(normalizedSignalSmooth[idxLower:idxUpper])
	diffYN = abs(maxYN-minYN)
	maxYN = maxYN + 0.5*diffYN
	minYN = minYN - 0.2*diffYN
	maxYNN = max(normalizedSignalSmooth[idxLower:idxUpper])


	junky = min(abs(lam[idxLower:idxUpper] - (fcent[iL]-3.0*ffwhm[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (fcent[iL]+3.0*ffwhm[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 
	dtmp = 1.0d - affp[iL]
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',afdp[iL],dtmp)

;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],$
			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))), xstyle=1,$
			yticks=3,ystyle=1,ytickv=[minYN,0,maxYNN,maxYN],$
			xtitle=lambdatitle, xtics=1, color='grey', yrange=[minYN, maxYN],ytitle='Pseudoflux'
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='blue',thick=2

    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[ffmax[iL],fcent[iL],ffwhm[iL]/2.3548]),color='green',thick=0.5, linestyle=0
		;;horizontal lines
	djs_oplot,[4047.0,4047.0],[-1000,10000],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[4358.0,4358.0],[-1000,10000],color='blue',linestyle=2 ;; 4358 HgI line

 	djs_xyouts,lam[idxLower+5],(0.8*maxYN),probinfo,charsize=0.8

	;-- sixth figure
;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
    djs_plot,lam,originalSignalBaseline,$
			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,5000],ytitle='Original Profile'

    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
		;;horizontal lines
	djs_oplot,[4047.0,4047.0],[-1000,10000],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[4358.0,4358.0],[-1000,10000],color='blue',linestyle=2 ;; 4358 HgI line


	print,string(format='("DUMPS = ",i5,"/",i5)',iL,numDump)
	ps_close
	;write_png,'./dumpsummarypictures/'+objname[iL]+'INFO.png',tvrd(true=1)



	;*********************
	endfor
	;*********************





	;================================================================

	;== write data : directory setting for  "dump.detection"
	if direxist('dumpdetectionpictures') eq 0L then spawn, 'mkdir dumpdetectionpictures'

	infilename = dobjname+'Profile.sav'


	;-- for dump.summary	
;	idxNoneDetect = where(scent lt 0.0)
;	scent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
;	idxNoneDetect = where(tcent lt 0.0)
;	tcent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	numDump = n_elements(infilename)
	numDump = numDump-1


    ;***********
    ;window,0,retain=2,xsize=900,ysize=1200
	charangstrom = '('+STRING(197B)+')'
    !p.charsize=1.3
    !p.multi=[0,2,3]


	iL=0L
	;***********
	for iL=0L,numDump do begin

	;===== prepare for eps output
    ps_open,'./dumpdetectionpictures/'+dobjname[iL]+'INFO',/ps,/encap,/color
    device,/times,xsiz=5.5,ysiz=6.5,/inch

	;write_png,'./dumpdetectionpictures/'+dobjname[iL]+'INFO.png',tvrd(true=1)


	;=============================== first figure
	djs_plot,d1cent,d2cent,psym=7,yr=[3820,5000],xr=[4000,5000],xstyle=1,ystyle=1,symsize=0.5,$
		xtitle=TexToIDL("\lambda : Primary Detection"),ytitle=TexToIDL("\lambda : Secondary Detection"), title='1st & 2nd detections'
;	djs_oplot,[4255,4255],[-1000,10000],color='blue'
;	djs_oplot,[4620,4620],[-1000,10000],color='blue'
	djs_oplot,[4358.0-dSL,4358.0-dSL],[-1000,10000],color='blue' ;; 4358 HgI line
	djs_oplot,[4358.0+dSL,4358.0+dSL],[-1000,10000],color='blue' ;; 4358 HgI line
	djs_oplot,[4047.0-dSL,4047.0-dSL],[-1000,10000],color='blue' ;; 4047 HgI line
	djs_oplot,[4047.0+dSL,4047.0+dSL],[-1000,10000],color='blue' ;; 4047 HgI line

		;;horizontal lines
	djs_oplot,[-1000,10000],[4047.0,4047.0],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[-1000,10000],[4358.0,4358.0],color='blue',linestyle=2 ;; 4358 HgI line
    djs_oplot,[-1000,10000],[4000.0,4000.0],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[4000.0,4000.0],[-1000,4000],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[5000.0,5000.0],[-1000,4000],color='black',linestyle=0,thick=4 ;; d100 line
    djs_oplot,[-1000,10000],[3820.0,3820.0],color='black',linestyle=0,thick=4 ;; d100 line


	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4358.0-dSL,4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL],color='red',thick=2
	djs_oplot,[4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL,4047.0-dSL],[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],color='red',thick=2
	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],color='red',thick=2
	djs_oplot,[4047.0-dSL,4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL],[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],color='red',thick=2

	;;;;;; overplot the current object
;    if d2cent[iL] gt 3950 then begin 
        djs_oplot,[d1cent[iL]],[d2cent[iL]],psym=sym(13),color='blue',symsize=2.0,thick=3
;    endif else begin 
;        djs_oplot,[d1cent[iL]],[3950],psym=sym(3),color='blue',symsize=2.0,thick=3
;    endelse

;	djs_oplot,[-1000,10000],[4620,4620],color='blue'
;	djs_oplot,[-1000,10000],[4358.0-dSL,4358.0-dSL],color='red' ;; 4358 HgI line
;	djs_oplot,[-1000,10000],[4358.0+dSL,4358.0+dSL],color='red' ;; 4358 HgI line
;	djs_oplot,[-1000,10000],[4047.0-dSL,4047.0-dSL],color='red' ;; 4047 HgI line
;	djs_oplot,[-1000,10000],[4047.0+dSL,4047.0+dSL],color='red' ;; 4047 HgI line

	;third detection
;	djs_plot,fcent,tcent,psym=7,xr=[4000,5000],yr=[3900,5000],xstyle=1,ystyle=1
;	djs_oplot,[4255,4255],[-1000,10000],color='blue'
;	djs_oplot,[4620,4620],[-1000,10000],color='blue'
;	djs_oplot,[4358.0-dSL,4358.0-dSL],[-1000,10000],color='red' ;; 4358 HgI line
;	djs_oplot,[4358.0+dSL,4358.0+dSL],[-1000,10000],color='red' ;; 4358 HgI line
;	djs_oplot,[4358.0-dSL,4358.0+dSL,4358.0+dSL,4358.0-dSL,4358.0-dSL],[4255,4255,4620,4620,4255],color='blue',thick=2
;	djs_oplot,[4047.0-dSL,4047.0+dSL,4047.0+dSL,4047.0-dSL,4047.0-dSL],[4255,4255,4620,4620,4255],color='blue',thick=2

	;============================= second figure: plot histogram of dump.detect

	lambdatitle = TeXtoIDL("\lambda ")+charangstrom

	cgHistoplot,dfcent,binsize=5.0,/fill, polycolor='magenta',datacolorname='magenta',$
        xtitle=lambdatitle, ytitle='Counts',yrange=[0,histymax],xrange=[4000,5000], title='Detection Histogram',xstyle=1
	djs_oplot,[dfcent[iL],dfcent[iL]],[-1,100],color='navy',linestyle=3, thick=2
	djs_oplot,nflam,nfthru,color='green',thick=2.0

	;============================ third figure: plot the detection probability
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Contour in IDL is crazy
	;; Follow the below correctly to avoid IDL stupidity
	maxZ = 1.0 
	minz = 0.0 
	nlevels=100
	cgloadct,3, ncolors=nlevels, bottom=1
	zSteps = (maxZ - minZ)/(double(nlevels))
	zLevels = indgen(nlevels)*zSteps + minZ
	SetDecomposedState, 0, CurrentState=currentState

	labelLevels = [0.95,0.99]
	labelNames = ['0.95','0.99']
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	iThres = 0 
	iWds = 0 

	;print, n_elements(fwhm)

	;threesigmaFlux = dindgen(n_elements(fwhm))
	;fivesigmaFlux = dindgen(n_elements(fwhm))
	threesigmaFlux = 3.0*sqrt(fwhm*2.54799)
	fivesigmaFlux = 5.0*sqrt(fwhm*2.54799)


	fratemoments = moment(fprob[iThres,iWds,*,*])
	fratestr = string(format='(a4,f6.4,a2,f9.6)',"F = ",fratemoments[0],string(177b),sqrt(fratemoments[1]))
	stitle = string(format='("D: (T,W) = (",f4.1,"," ,f4.1,")")',(dThreshold+iThres*0.1),(dWidth+iWds))
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/fill,xtitle=xtstr,levels=zLevels, $
;    	c_colors=indgen(nlevels)+1, nlevels=nlevels, $
;   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1
;	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=labelLevels, $
;    color=cgcolor('blue'), c_labels=replicate(1,2), $
;    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
;    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
;		xtitle=xtstr,levels=[0.01,0.05,0.5], $
;        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
;         c_labels=replicate(0,3), $
;        c_linestyle=0,c_charsize=0.5,/follow
    contour,falseDetectionOneSig,xdensfunc,ydensfunc   ,$
;        /overplot,$
		xtitle=xtstr,levels=[0.01,0.05,0.5], $
        c_color=cgcolor('aqua'),c_thick=[1,1,4],$
         c_labels=replicate(0,3), $
        c_linestyle=0,c_charsize=0.5,/follow,$
   		ytitle=ytstr,title=stitle,xrange=[0,80],yrange=[1,15],ystyle=1,xstyle=1	
	contour, dProb[iThres,iWds,*,*],trueFlux,fwhm,/overplot,xtitle=xtstr,levels=[0.5,0.95,0.99], $
    color=cgcolor('blue'), c_labels=replicate(0,3), c_thick=[4,1,1], $
    ytitle=ytstr,title=stitle, c_linestyle=1,c_charsize=0.5,/follow
	djs_oplot,dffmax*dffwhm*1.06447/dldp,dffwhm/dldp,psym=7,symsize=0.5,color='magenta'


	;;;
	tmpflux = dffmax[iL]*dffwhm[iL]*1.06447/dldp
	tmpfwhm =dffwhm[iL]/dldp
    tmpsymnum=9
    if (tmpflux lt 80.0d) and (tmpfwhm lt 15.0d) then begin
        djs_oplot,[tmpflux],[tmpfwhm],psym=sym(13),symsize=2.0,color='blue',thick=3.0
    endif else begin
        if tmpfwhm gt 15.0d then begin
            tmpfwhm = 14.5d
            tmpsymnum=2
        endif
        if tmpflux gt 80.0d then begin
            tmpflux = 79.0d
            tmpsymnum=14
        endif
        djs_oplot,[tmpflux],[tmpfwhm],psym=sym(tmpsymnum),symsize=2.0,color='blue',thick=5.0
    endelse


	;========================= fourth figure: plot RA-DEC position

	djs_plot,dra,ddec,color='magenta',psym=1,xr=[max(dra),min(dra)],yr=[min(ddec),max(ddec)],$
		symsize=1.0,xtitle='RA', ytitle='DEC', title='Spatial Distribution'
	djs_oplot,[dra[iL]],[ddec[iL]],psym=sym(13),symsize=2.0,color='blue',thick=3.0


	;========================= fifth figure
    restore,infilename[iL]
	;print,infilename[iL]," *** iL = ",iL," size idxWidth = ",size(ddidxWidth)
	;print,"idxStart[iL] ddidxWidth[iL] = ",idxStart[iL]," ",ddidxWidth[iL]
    idxCenter = long(dfidx[iL] + 0.5 * dfwd[iL])
	;print,infilename[iL]," idxCenter = ", idxCenter
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower lt 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 0.5*diffY
	minY = minY - 0.2*diffY
	maxYY = max(originalSignalBaseline[idxLower:idxUpper])

	maxYN = max(normalizedSignalSmooth[idxLower:idxUpper])
	minYN = min(normalizedSignalSmooth[idxLower:idxUpper])
	diffYN = abs(maxYN-minYN)
	maxYN = maxYN + 0.5*diffYN
	minYN = minYN - 0.2*diffYN
	maxYNN = max(normalizedSignalSmooth[idxLower:idxUpper])


	junky = min(abs(lam[idxLower:idxUpper] - (dfcent[iL]-3.0*dffwhm[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (dfcent[iL]+3.0*dffwhm[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 
	dtmp = 1.0d - dffp[iL]
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',dfdp[iL],dtmp)

;    djs_oplot,lam[idxgLower:idxgUpper],gaussian(lam[idxgLower:idxgUpper],[gmax[iL],gcentroid[iL],gsigma[iL]/2.3548]),color='black',thick=0.7, linestyle=0
    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],xtitle=lambdatitle,ytitle='Pseudoprofile',$
			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))), xstyle=1,$
			yticks=3,ystyle=1,ytickv=[minYN,0,maxYNN,maxYN],$
			xtics=1, color='grey', yrange=[minYN, maxYN]
    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='blue',thick=2

     djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[dffmax[iL],dfcent[iL],dffwhm[iL]/2.3548]),color='green',thick=0.5, linestyle=0
		;;horizontal lines
	djs_oplot,[4047.0,4047.0],[-1000,10000],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[4358.0,4358.0],[-1000,10000],color='blue',linestyle=2 ;; 4358 HgI line

    djs_xyouts,lam[idxLower+5],(0.8*maxYN),probinfo,charsize=0.8


	;;------------- sixth figure
;    djs_plot,lam[idxLower:idxUpper],originalSignalBaseline[idxLower:idxUpper],$
    djs_plot,lam,originalSignalBaseline,xtitle=lambdatitle,ytitle='Original Profile',$
			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))), xstyle=1,$
			 xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,5000]
    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2
		;;horizontal lines
	djs_oplot,[4047.0,4047.0],[-1000,10000],color='blue',linestyle=2 ;; 4047 HgI line
	djs_oplot,[4358.0,4358.0],[-1000,10000],color='blue',linestyle=2 ;; 4358 HgI line

	print,string(format='("DETECTIONS = ",i5,"/",i5)',iL,numDump)
	ps_close



	;interrupt dummy
	;dummystring =''
	;read,dummystring , prompt='e for exit : '
	;if strupcase(dummystring) eq 'E' then stop 



	;write_png,'./dumpdetectionpictures/'+dobjname[iL]+'INFO.png',tvrd(true=1)
	;*********************
	endfor
	;*********************








;============
; making plots for matching emission lines 

	;== write data : directory setting
	if direxist('dumpmatchpictures') eq 0L then spawn, 'mkdir dumpmatchpictures'

	;-- for finaldump.summary
	idxNoneDetect = where(scent lt 3000.0)
	scent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
	idxNoneDetect = where(tcent lt 3000.0)
	tcent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	idxNoneDetect = where(d2cent lt 3000.0)
	d2cent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	;-- for dump.summary
;	idxNoneDetect = where(scent lt 0.0)
;	scent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
;	idxNoneDetect = where(tcent lt 0.0)
;	tcent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0

	infilename = objname+'Profile.sav'
	inwholefilename = objname+'WholeProfile.sav'
	inwholedetectionfilename = objname+'D.cat'
	numDump = n_elements(infilename)
	numDump = numDump-1

	iL=0L



;	restore,'/Users/shong/work/newdata/monte/nosmoothfiterror/DFcontour.sav' ;; hard-wired
;	restore,'/Users/shong/work/newdata/monte/nosmoothresults/monte_results_refine_fluxscale_nosmooth.sav' ;; hard-wired
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;***********
	;window,0,retain=2,xsize=900,ysize=1200
	charangstrom = '('+STRING(197B)+')'
	!p.charsize=1.3
   	!p.multi=[0,1,4]

	for iL=0L,numDump do begin
	
	;===== prepare for eps output
	ps_open,'./dumpmatchpictures/'+objname[iL]+'MATCH',/ps,/encap,/color
    device,/times,xsiz=7.0,ysiz=6.6,/inch

	;write_png,'./dumpsummarypictures/'+objname[iL]+'INFO.png',tvrd(true=1)


	;========================= fifth figure
    restore,infilename[iL]
    restore,inwholefilename[iL]
	readcol,inwholedetectionfilename[iL],f='x,x,x,x,x,x,x,x,x,f',allcentroids
	;print,infilename[iL]," *** iL = ",iL," size idxWidth = ",size(ddidxWidth)
	;print,"idxStart[iL] ddidxWidth[iL] = ",idxStart[iL]," ",ddidxWidth[iL]
    idxCenter = long(fidx[iL] + 0.5 * fwd[iL])
	;print,infilename[iL]," idxCenter = ", idxCenter
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower lt 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 1.0*diffY
	minY = minY - 0.5*diffY

	maxYN = max(normalizedSignalSmooth[idxLower:idxUpper])
	minYN = min(normalizedSignalSmooth[idxLower:idxUpper])
	diffYN = abs(maxYN-minYN)
	maxYN = maxYN + 0.5*diffYN
	minYN = minYN - 0.2*diffYN
	maxYNN = max(normalizedSignalSmooth[idxLower:idxUpper])


	junky = min(abs(lam[idxLower:idxUpper] - (fcent[iL]-3.0*ffwhm[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (fcent[iL]+3.0*ffwhm[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 
	dtmp = 1.0d - affp[iL]
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',afdp[iL],dtmp)

;    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))), xstyle=1,$
;			yticks=3,ystyle=1,ytickv=[minYN,0,maxYNN,maxYN],$
;			xtitle=lambdatitle, xtics=1, color='grey', yrange=[minYN, maxYN]
;    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='blue',thick=2

;    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[ffmax[iL],fcent[iL],ffwhm[iL]/2.3548]),color='green',thick=0.5, linestyle=0
;    djs_xyouts,lam[idxLower+5],(0.8*maxYN),probinfo,charsize=0.8


;    djs_plot,lam,originalSignalBaseline,$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
;			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2


	;;;;first	
	idxgLower += wholeXminIdx
	idxgUpper += wholeXminIdx

    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,6500],ystyle=1
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2
    djs_xyouts,4200,(0.8*maxY),probinfo,charsize=0.8
    djs_xyouts,5600,(0.8*maxY),strupcase(strjoin(strsplit(objname[iL],'_',/extract))),charsize=0.8

	;;;;second
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,6500],ystyle=1
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1,thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (fcent[iL] - 1215.67)/1215.67 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick=3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 
	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor


	;-------------------------------



	;;;; third
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,6500],ystyle=1 
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1, thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (fcent[iL] - 3727.08)/3727.08 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick =3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 

	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor
	;-------------------------------


	;;;; fourth
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,6500],ystyle=1 
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1, thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (fcent[iL] - 1549.053)/1549.053 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick =3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 

	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor
	;-------------------------------










	print,string(format='("DUMPS = ",i5,"/",i5)',iL,numDump)
	ps_close
	;write_png,'./dumpsummarypictures/'+objname[iL]+'INFO.png',tvrd(true=1)



	;*********************
	endfor
	;*********************





;============
; making plots for matching emission lines 

	;== write data : directory setting
	if direxist('dumpdetectmatchpictures') eq 0L then spawn, 'mkdir dumpdetectmatchpictures'

	;-- for d1 d2 nondetect

	idxNoneDetect = where(d1cent lt 3000.0)
	d1cent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0
	idxNoneDetect = where(d2cent lt 3000.0)
	d2cent[idxNoneDetect] = 3905.0; for the case of no secondary detection, we set the centroid lambda 3999.0


	infilename = dobjname+'Profile.sav'
	inwholefilename = dobjname+'WholeProfile.sav'
	inwholedetectionfilename = dobjname+'D.cat'
	numDump = n_elements(infilename)
	numDump = numDump-1

	iL=0L



;	restore,'/Users/shong/work/newdata/monte/nosmoothfiterror/DFcontour.sav' ;; hard-wired
;	restore,'/Users/shong/work/newdata/monte/nosmoothresults/monte_results_refine_fluxscale_nosmooth.sav' ;; hard-wired
	xtstr = TeXtoIDL("Pseudoflux")
	ytstr = TeXtoIDL("FWHM")
	stitle = ''
	fratestr = ''

	;***********
	;window,0,retain=2,xsize=900,ysize=1200
	charangstrom = '('+STRING(197B)+')'
	!p.charsize=1.3
   	!p.multi=[0,1,4]

	for iL=0L,numDump do begin
	
	;===== prepare for eps output
	ps_open,'./dumpdetectmatchpictures/'+dobjname[iL]+'MATCH',/ps,/encap,/color
    device,/times,xsiz=7.0,ysiz=6.6,/inch

	;write_png,'./dumpsummarypictures/'+objname[iL]+'INFO.png',tvrd(true=1)


	;========================= fifth figure
    restore,infilename[iL]
    restore,inwholefilename[iL]
	readcol,inwholedetectionfilename[iL],f='x,x,x,x,x,x,x,x,x,f',allcentroids
	;print,infilename[iL]," *** iL = ",iL," size idxWidth = ",size(ddidxWidth)
	;print,"idxStart[iL] ddidxWidth[iL] = ",idxStart[iL]," ",ddidxWidth[iL]
    idxCenter = long(dfidx[iL] + 0.5 * dfwd[iL])
	;print,infilename[iL]," idxCenter = ", idxCenter
    idxLower = idxCenter - 50
    idxUpper = idxCenter + 50

    if idxLower lt 0 then begin
        idxLower = 0
    endif

    if idxUpper ge n_elements(lam) then begin
        idxUpper = n_elements(lam)-1
    endif

	maxY = max(originalSignalBaseline[idxLower:idxUpper])
	minY = min(originalSignalBaseline[idxLower:idxUpper])
	diffY = abs(maxY-minY)
	maxY = maxY + 1.0*diffY
	minY = minY - 0.5*diffY

	maxYN = max(normalizedSignalSmooth[idxLower:idxUpper])
	minYN = min(normalizedSignalSmooth[idxLower:idxUpper])
	diffYN = abs(maxYN-minYN)
	maxYN = maxYN + 0.5*diffYN
	minYN = minYN - 0.2*diffYN
	maxYNN = max(normalizedSignalSmooth[idxLower:idxUpper])


	junky = min(abs(lam[idxLower:idxUpper] - (dfcent[iL]-3.0*dffwhm[iL]/2.3548)),idxgLower)
	junky = min(abs(lam[idxLower:idxUpper] - (dfcent[iL]+3.0*dffwhm[iL]/2.3548)),idxgUpper)
	idxgLower = idxgLower + idxLower
	idxgUpper = idxgUpper + idxLower
	
	;print,iL," : ", tmpIndex, epsIndex, " Indexes: ", idxLower, idxCenter, idxUpper," Gaussian: ", gcentroid[iL], gsigma[iL],idxgLower,idxgUpper, n_elements(lam) 
	dtmp = 1.0d - dffp[iL]
	probinfo=string(format='("(C,R) = (",f5.3,", " ,f5.3,")")',dfdp[iL],dtmp)

;    djs_plot,lam[idxLower:idxUpper],normalizedSignalSmooth[idxLower:idxUpper],$
;			title=strupcase(strjoin(strsplit(objname[iL],'_',/extract))), xstyle=1,$
;			yticks=3,ystyle=1,ytickv=[minYN,0,maxYNN,maxYN],$
;			xtitle=lambdatitle, xtics=1, color='grey', yrange=[minYN, maxYN]
;    djs_oplot,lam[idxgLower:idxgUpper],normalizedSignalSmooth[idxgLower:idxgUpper],color='blue',thick=2

;    djs_oplot,lam[idxLower:idxUpper],gaussian(lam[idxLower:idxUpper],[ffmax[iL],fcent[iL],ffwhm[iL]/2.3548]),color='green',thick=0.5, linestyle=0
;    djs_xyouts,lam[idxLower+5],(0.8*maxYN),probinfo,charsize=0.8


;    djs_plot,lam,originalSignalBaseline,$
;			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),$
;			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY]
;    djs_oplot,lam[idxgLower:idxgUpper],originalSignalBaseline[idxgLower:idxgUpper],color='black',thick=2


	;;;;first	
	idxgLower += wholeXminIdx
	idxgUpper += wholeXminIdx

    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,7500],ystyle=1
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2
    djs_xyouts,4200,(0.8*maxY),probinfo,charsize=0.8
    djs_xyouts,5600,(0.8*maxY),strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),charsize=0.8

	;;;;second
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,7500],ystyle=1
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1,thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (dfcent[iL] - 1215.67)/1215.67 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick=3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 
	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor


	;-------------------------------



	;;;; third
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,7500],ystyle=1 
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1, thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (dfcent[iL] - 3727.08)/3727.08 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick =3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 

	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor
	;-------------------------------


	;;;; fourth
    djs_plot,wholelam,wholeprofile,$
;			title=strupcase(strjoin(strsplit(dobjname[iL],'_',/extract))),$
			xtitle=lambdatitle, xstyle=1, xtics=1,yticks=3, color='grey', yrange=[minY, maxY],xrange=[4000,7500],ystyle=1 
    djs_oplot,wholelam[idxgLower:idxgUpper],wholeprofile[idxgLower:idxgUpper],color='black',thick=2


	;----------plot all detected centroids
           for j = 1, n_elements(allcentroids) -1 do begin
              djs_oplot, [allcentroids(j), allcentroids(j)], [-1e4,1e5], $
                 color='grey', linestyle=1, thick=3
           endfor 
	;----------plot emission lines
		   diff = (maxY - minY)
		   z = (dfcent[iL] - 1549.053)/1549.053 ;; change it into redshift scale
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si IV', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='blue', linestyle=2, thick =3
              djs_xyouts, (1+Z)*lines(j), minY+diff*(0.1+0.06*(-1)^j), $
                 labels(j), charsize=0.7, color='blue'
           endfor 

	;-----------plot skyem

           skylines= [4046.56,4165,4168,4358.34,4420,4423,$
                4465,4469, 4827, 4832,4983,5199,5461  ,5577.34,5685, 5889.95,5895.92,6300.64,6363.78]
           skylabels = ['HgI(4047)','NaI(4165+68)',' ','HgI(4358)','NaI(4420+23)',' ',$
                'NaI(4465+69)',' ','HgI(4827+32)',' ','NaI(4983)','NI(5199)','HgI(5461)',$
                'OI(5577)', 'NaI(5683+88)', 'NaI(5890)',$
                'NaI(5896)','[OI](6301)','[OI](6363)']
           for j = 0, n_elements(skylines) -1 do begin
              djs_oplot, (1+0.0)* [skylines(j), skylines(j)], [-1e4,1e5], $
                 color='green', linestyle=1,thick=3
              ;djs_xyouts, (1+0.0002)*skylines(j), maxY-diff*(0.3+0.04*(j mod 5)), $
              ;   skylabels(j), charsize=0.7, color='green'
           endfor
	;-------------------------------










	print,string(format='("DUMPS = ",i5,"/",i5)',iL,numDump)
	ps_close
	;write_png,'./dumpsummarypictures/'+dobjname[iL]+'INFO.png',tvrd(true=1)



	;*********************
	endfor
	;*********************

















endif ; skipdump ends
;*********

end

