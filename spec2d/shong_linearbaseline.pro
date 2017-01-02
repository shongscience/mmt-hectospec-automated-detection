function shong_linearbaseline, lamArr, signalArr, nSigmaArr
;-----------------------------
;Sungryong Hong (1/*/2012)
;
; return the baseline-subtracted signal
;
;Need:
;	mpfit packages
;
;-----------------------------

	;---------------------------------------------------------------
	; EXIT if the dimenstions of lamArr and nSigmaArr are different
	if (n_elements(lamArr) eq n_elements(nSigmaArr)) then begin 

	;-------------------------
	; Initialize arrays
	;-------------------------
	arrSize = n_elements(lamArr) 
	lam = lamArr ;; assign from input lambda
	signal = signalArr
	noiseSigma = nSigmaArr ;; assign from input noise sigma	
	mask = indgen(arrSize)                                                                                                
	baseline = dindgen(arrSize)                                                                                                
	mask[*] = 0

	;------------------------------------------------
	;Temporal input, instead of inputs from arguments
	;lam = dindgen(1000)                                                                                                
	;signal = dindgen(1000)                                                                                             
	;noiseSigma = 3.0D + randomn(seed, 1000)                                                                                  
	;signal = 10.0D*exp(-1.0*(lam-150.)*(lam-150.)/(2.0*4.0*4.0))+ 15.0D*exp(-1.0*(lam-450.)*(lam-450.)/(2.0*6.0*6.0)) 
	;signal = signal + randomn(seed,1000) * noiseSigma    
	

	;------------------------------------------------
	;Calculate the mean and stdev for the signal array
	tmpStdev = 0.D
	tmpMean = 0.D
	upsig = 0.D
	downsig = 0.D
	tmpStdev = stdev(signal,tmpMean)
	
	;print, tmpStdev,tmpMean
	;print, size(tmpStdev)

	
	;------------------------------------------------
	;Mask out +- 2 sigma outliers
	upsig = tmpMean + 2.0 * tmpStdev
	downsig = tmpMean - 2.0 * tmpStdev

	tmpidx = where( signal[*] gt upsig )
	mask[tmpidx] = 1
	tmpidx = where( signal[*] lt downsig )
	mask[tmpidx] = 1 
	
	maskIdx = where( mask[*] eq 0)


	;---------------------------------------------------------------
	; MPFIT on the masked profile : reject outliers above +- 2 sigma
	;	1) fit using linear profile
	;	2) remove the base line from signal
	fitfunc = 'P[0] + X * P[1]'
	initValues = [0.0D,0.0D]
	print, "<<<Base Line Fit>>>"
	para = mpfitexpr(fitfunc,lam[maskIdx], signal[maskIdx],noiseSigma[maskIdx],initValues)
	baseline = para(0) + para(1)*lam



	;-------------------------
	; Plot the related files 
	;-------------------------
   	goplot=0
    if (goplot) then begin
		window,15,retain=2,xsize=800,ysize=1100
	    !p.charsize=2.0
	    !p.multi=[0,1,3]

		djs_plot, lam, signal, xtitle='lambda', ytitle='Signal and Sigma'
		djs_oplot, lam, noiseSigma, color='red'
		djs_xyouts, (arrSize/20.),(0.9*max(signal)),'Input Signal + Sigma', color='green'


		djs_plot, lam, signal, xtitle='lambda', ytitle='Signal and Sigma'
		djs_oplot, lam[maskIdx], signal[maskIdx], color='blue'
		djs_oplot, lam, para(0) + para(1)*lam, color='green', thick=3
		djs_xyouts, (arrSize/20.),(0.9*max(signal)),'Base Noise + Function Fit', color='green'

		djs_plot, lam, signal-baseline, xtitle='lambda', ytitle='Signal and Sigma'
		djs_oplot, lam, noiseSigma, color='red'
		djs_xyouts, (arrSize/20.),(0.9*max(signal-baseline)),'Baseline subtracted + Sigma', color='green'

        ;;;;interupt
        voidstring = ''
        read, voidstring, prompt='<<<<enter to exit baseline fit>>>>>'

		wdelete, 15
		!p.multi=0
	endif

	
	; EXIT if the dimenstions of lamArr and nSigmaArr are different
	endif else begin 
		print, "Error: Input arrays are not compatible in dimensions. Check size(inputArray)"
	endelse
	;-----------------------------------------------------------------------------------------

	;baseline = signal - baseline ;; return baseline instead of baseline-subtracted signal
	return, baseline

end 

