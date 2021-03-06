function mygauss,x,p
    reval = 0.0d
    reval = gaussian(x,p[0:2]) ;; gaussian.pro p[0]=max,p[1]=centroid,p[2]=sigma 
    return,reval
end

pro shong_dumpNoDetectionNoSmoothGoodFit,plugmap,dTh,dWd, lamArr, signalArr, nSigmaArr, goplot=goplot, interactive=interactive
;-----------------------------
;Sungryong Hong (1/*/2012)
;	Input: 
;		1) plugmap: object logs, dTh: detection threshold, dWd: detection width, 
;					lamArr: lambda info, signalArr: signal info, nSigmaArr: sigma info 
;	output:
;		1) filename + each detection: idl binary => sav format
;		2) lamArr, signalArr, baseline subtracted signal Arr, detection array, 
;			detected width, start idx, end idx
;	Comments:
;	 	For the given dTh and dWd, save as ``sav'' file for further analysis
;		Automatically create ./nodetection/ directory
;	
;		** this is a modified version of dumpDetection with smoothing length=0.2 . This is zero smoothing. or deltafunction smooth
;
;	Need:
;		type '.run compileAllShong.pro' for necessary pro files
;-----------------------------

	if not keyword_set(goplot) then goplot=0 ;; default is visualizing the combining processes
	if not keyword_set(interactive) then interactive=0 ;; default is visualizing the combining processes

	;---------------------------------------------------------------
	; EXIT if the dimenstions of lamArr and nSigmaArr are different
	if (n_elements(lamArr) eq n_elements(nSigmaArr)) then begin 

	;-------------------------
	; Initialize arrays
	;-------------------------
	arrSize = n_elements(lamArr) 
	lam= lamArr ;; assign from input lambda
	signal= signalArr
	signalOriginal = signalArr
	signalSmooth= signalArr
	normalizedSignalSmooth= signalArr ;; signal / nsig ;; normalized S/N
	noise=dindgen(arrSize)	
	noiseSmooth=dindgen(arrSize)	
	noiseSigma=nSigmaArr ;; assign from input noise sigma	
	noiseSigmaOriginal=nSigmaArr ;; assign from input noise sigma	
	noiseSigmaSmooth=nSigmaArr ;; assign from input noise sigma	
	baseLine=dindgen(arrSize)	
	mixedSignal=dindgen(arrSize)	
	mixedSignalSmooth=dindgen(arrSize)	
	detectionFlag=dindgen(arrSize)	
	thresholdLine=dindgen(arrSize)	
	zeroLine=dindgen(arrSize)	
	dtmparr=dindgen(arrSize)	
	detectionFlag[*] = 0.0D ;; initialize the detection array
	dThreshold  = dTh
	dWidth = dWd

	if direxist('nodetection') eq 0L then spawn, 'mkdir nodetection'

	;Temporal input, not from arguments
	;lam = dindgen(1000)                                                                                                
	;nsig = dindgen(1000)                                                                                               
	;signal = dindgen(1000)                                                                                             
	;noiseSigma = 3.0D + randomn(seed, 1000)                                                                                  
	;signal = 5.0D*exp(-1.0*(lam-150.)*(lam-150.)/(2.0*4.0*4.0))+ 7.0D*exp(-1.0*(lam-450.)*(lam-450.)/(2.0*6.0*6.0)) 
	;signal = signal + randomn(seed,1000) * noiseSigma    
	;signalSmooth= signal
	;normalizedSignalSmooth= signal
	;noiseSigmaSmooth=noiseSigma ;; assign from input noise sigma	
	



	;; Assign a normal noise (1) assign sigma value
	;noiseSigma = 5.0D ;; assign from input argument... So, an obsolete line
	;; -----------------------------------------------------




	;-------------------------------------------------------------
	; Subtract the baseline 
	baseLine = shong_linearbaseline(lam,signal,noiseSigma)
	originalSignal = signal
	signal = signal - baseLine 
	originalSignalBaseline = signal

	signal = signal/noiseSigma ;; normalize the signal to S/N
	noiseSigma[*] = 1.0D ;; now sigma is normalized to 1.0D 
	signalSmooth = signal ;;signal will be smoothed

	;; set iteration and related quantities
	iterL = 0L
	iDetect = 0L
	iDetectMax = 0L
	dtmp = 0.0D
	dtmpNoiseSigmaMax = max(noiseSigma)



	;goplot=0
	;-------------------------
	; Plot the related files 
	;-------------------------
	if goplot eq 1 then begin
   	window,10,retain=2,xsize=1800,ysize=1000
    !p.charsize=2.0
    !p.multi=[10,2,3]

	wset,10
	endif
	;--------------------------------
	; Create smooth kernels
	;; SMOOTH 
	smoothSigma = 0.2D ;; this make a delta function kernel
	gaussian_ker = psf_gaussian([1.0,floor(2.5*smoothSigma),smoothSigma],npixel=floor(5.0*smoothSigma)) ;5sigma Gaussian profile
    gaussian_ker = gaussian_ker/total(gaussian_ker) ; normalize 
	gaussian_ker2 = gaussian_ker * gaussian_ker 	
	noiseSigmaSmooth = noiseSigma * noiseSigma ;; variance
	noiseSigmaSmooth = convol(noiseSigmaSmooth,gaussian_ker2,/center) ;; smoothed variance
	noiseSigmaSmooth = sqrt(noiseSigmaSmooth) ;; change variance to std dev 
	signalSmooth = convol(signalSmooth,gaussian_ker,/center) ;; smoothed signal
	normalizedSignalSmooth = signalSmooth/noiseSigmaSmooth
	;--------------------------------




	;-------------------------------------------------------------------------------------------------------------
	; Show results and fitting the detected regions
	;
	



		        detectionFlag[*] = 0.0D
		        thresholdLine[*] = dThreshold
		        zeroLine[*] = 0.0D
		        detectedPixels = where(normalizedSignalSmooth ge dThreshold,count)
		        undetectedPixels = where(normalizedSignalSmooth lt dThreshold,count)
		        ;detectionFlag[detectedPixels] = 1.0D
		        ;numDetected = n_elements(detectedPixels)
		        detections = findAllClusterWithWidthCut(detectedPixels, dWidth)

				print, "<<<<<<<< Detections >>>>>>>>"
				print, "Object Name: ", plugmap.objtype
				;print, "size :", size(detections)
				;print, "contents :"
				;print,  detections


				;; Trime detections above dWidth
				


				
				iDetect = 0
				iDetectMax = 0

			if (long((size(detections))[0]) eq 1) then begin ;; single detection, so 1D array
				if (detections[0] eq -1) then begin
					iDetectMax = 0
					print,"**NON detection***"
				endif else begin
					iDetectMax = 1
				endelse
			endif else begin                                 ;; multiple detection, so 2D array
				iDetectMax = long((size(detections))[2])
			endelse
			
			;------Define a new detection array------
			newDetections = lindgen(2,1000) ;; Make a new detection array which has a[0,0] = a[1,0] = Number of detection
											;; a[0,i] = detectionIdx, a[1,i] = detectionWidth 
											;; a[0,0] = a[1,0] = 1 for single detection
											;; a[0,0] = a[1,0] = 0 for non-detection

			newDetections[0,0] = iDetectMax 
			newDetections[1,0] = iDetectMax 

			iDetectMax -= 1 
			if (iDetectMax gt 0) then begin ;; if multiple detection, use 2D data
				for iDetect=0L, iDetectMax do begin
					newDetections[0,iDetect+1] = detections[0,iDetect]
					newDetections[1,iDetect+1] = detections[1,iDetect]
				endfor
			endif else begin 
				if (detections[0] ne -1) then begin
					newDetections[0,1] = detections[0]
					newDetections[1,1] = detections[1]
				endif
			endelse



			
	;;;;------------ when no detection
	if iDetectMax lt 0 then begin

				savfilename=strcompress(plugmap.objtype,/remove_all)+'specProfile'+'.sav'
				savfilename='./nodetection/'+savfilename
				save,lam,originalSignalBaseline,normalizedSignalSmooth,noiseSigmaOriginal,noiseSigmaSmooth,filename=savfilename
				print,">> Save profiles to : "+savfilename
	endif ;------------------not nondetect




	endif else begin ;-----------------------------------------------------------------------------------------------------
	;EXIT Control
		print, "ERROR: The lambda and nSigma arrays are not compatible in dimension: Check size(arr)" 
	endelse


	;if strupcase(isWrite) eq 'Y' then begin ;; if file writing flag, 
	;endif
	

	;----------------------------------------------------------------------------------------------------------------------
end

;-----------------------------------------------------------------------------------
; from the index data by where(DetectionArray), find detected clumps
function findAllClusterWithWidthCut, arr, widthCut
	tempidxCurrent = 0L
	tempidxPrevious = 0L
	tempidx = 0L
	tempwidth = 0L
	chunk = lindgen(2,1000) ;; the maximum detections is 1000
	iDetection=0L
	chunk[*,*] = 0L
	dWidth = long(widthCut)

	n_arr = long(n_elements(arr))
	iL=0L
	
	if (n_arr gt 0) then begin 
		iDetection=0L
		chunk[0,iDetection] = arr[0]
		chunk[1,iDetection] = 1L
		tempwidth = 1L
		tempidx = arr[0]
		tempidxPrevious = arr[0]		
		iL=1L
		while (iL lt n_arr) do begin
			tempidxCurrent = arr(iL)
			if (tempidxCurrent eq (tempidxPrevious+1)) then begin
				tempwidth = tempwidth +1
			endif else begin ;; update cluster info if this is a new larger cluster
				if (tempwidth ge dWidth) then begin 				
					;print, "tempwidth, dWidth :", tempwidth, dWidth
					chunk[0,iDetection] = tempidx
					chunk[1,iDetection] = tempwidth
					iDetection=iDetection + 1
				endif
				tempidx = arr[iL] ;; initialize the info with the new cluster
				tempwidth = 1L 
			endelse 
			tempidxPrevious = tempidxCurrent
			iL = iL + 1
		endwhile
		
		;finally 
		if (tempwidth ge dWidth) then begin 				
			;print, "finally: tempwidth, dWidth :", tempwidth, dWidth
			chunk[0,iDetection] = tempidx
			chunk[1,iDetection] = tempwidth
		endif else begin
			iDetection = iDetection -1 ;; removed `+1'ed iDetection 
		endelse 		
	endif else begin
		iDetection = -1 ;; when null input array, then set iDetection negative
	endelse

	;print,"In FindFunction, iDetection = ", iDetection	

	if (iDetection eq -1) then begin
		chunk[0,0] = -1
		chunk[1,0] = -1
		return, chunk[*,0]
	endif else begin
		if (iDetection eq 0) then begin
			return, chunk[*,0]
		endif else begin	
			return, chunk[*,0:iDetection]
		endelse
	endelse
end 





;-----------------------------------------------------------------------------------
; from the index data by where(DetectionArray), find detected clumps
function findAllCluster, arr 
	tempidxCurrent = 0L
	tempidxPrevious = 0L
	tempidx = 0L
	tempwidth = 0L
	chunk = lindgen(2,1000) ;; the maximum detections is 1000
	iDetection=0L
	chunk[*,*] = 0L

	n_arr = long(n_elements(arr))
	iL=0L
	
	if (n_arr gt 0) then begin 
		iDetection=0L
		chunk[0,iDetection] = arr[0]
		chunk[1,iDetection] = 1L
		tempwidth = 1L
		tempidx = arr[0]
		tempidxPrevious = arr[0]		
		iL=1L
		while (iL lt n_arr) do begin
			tempidxCurrent = arr(iL)
			if (tempidxCurrent eq (tempidxPrevious+1)) then begin
				tempwidth = tempwidth +1
			endif else begin ;; update cluster info if this is a new larger cluster
				chunk[0,iDetection] = tempidx
				chunk[1,iDetection] = tempwidth
				iDetection=iDetection + 1
				tempidx = arr[iL] ;; initialize the info with the new cluster
				tempwidth = 1L 
			endelse 
			tempidxPrevious = tempidxCurrent
			iL = iL + 1
		endwhile
		
		;finally 
		chunk[0,iDetection] = tempidx
		chunk[1,iDetection] = tempwidth
	endif else begin
		iDetection = -1 ;; when null input array, then set iDetection negative
	endelse

	;print,"In FindFunction, iDetection = ", iDetection	

	if (iDetection eq -1) then begin
		chunk[0,0] = -1
		chunk[1,0] = -1
		return, chunk[*,0]
	endif else begin
		if (iDetection eq 0) then begin
			return, chunk[*,0]
		endif else begin	
			return, chunk[*,0:iDetection]
		endelse
	endelse
end 



;-----------------------------------------------------------------------------------
; from the index data by where(DetectionArray), find the largest chunk of detection
function findLargestCluster, arr 
	tempidxCurrent = 0L
	tempidxPrevious = 0L
	tempidx = 0L
	tempwidth = 0L
	chunk = lindgen(2)
	chunk(0) = 0L
	chunk(1) = 0L

	n_arr = long(n_elements(arr))
	iL=0L
	
	if (n_arr gt 0) then begin 
		chunk(0) = arr[0]
		chunk(1) = 1L
		tempwidth = 1L
		tempidx = arr[0]
		tempidxPrevious = arr[0]		
		iL=1L
		while (iL lt n_arr) do begin
			tempidxCurrent = arr(iL)
			if (tempidxCurrent eq (tempidxPrevious+1)) then begin
				tempwidth = tempwidth +1
			endif else begin ;; update cluster info if this is a new larger cluster
				if (tempwidth gt chunk(1)) then begin 
					chunk(0) = tempidx
					chunk(1) = tempwidth
				endif 
				tempidx = arr[iL] ;; initialize the info with the new cluster
				tempwidth = 1L 
			endelse 
			tempidxPrevious = tempidxCurrent
			iL = iL + 1
		endwhile
		
		;finally 
		if (tempwidth gt chunk(1)) then begin 
			chunk(0) = tempidx
			chunk(1) = tempwidth
		endif 
	endif

	return, chunk
end 

function shong_histogramRead, outHisto, inFraction
;--------------------------------------------------
; Read a fraction value from the histogram output by sungryong hong in 2012
; Input: 
;		1) the output histogram array 2x1000 from shong_hisgtogram.pro, 
;		2) the fraction such as 99% (0.99) or 99.9%(0.999) 
;			At which, we will read the flux value for that fraction
;
; Output:
; 		the integer(long) value between 0 - 999; normalized to 1000 for maximum
; Comments:
;		if max = min, histogram is trivial. In this case, this returns 0.0
; Example:
;		inData = randomn(seed,10000)
;		histoData = shong_histogram(inData)
;
;		print, shong_histogramRead(histoData, 0.95) ;; return the normalized flux out of 1000
;--------------------------------------------------

	reIdx = 0L

	if (outHisto[0,0] gt -0.5) then begin ;; outHisto[0,0] = -1.0 if max = min.. this case return 0.0
		junk = min(abs(outHisto[1,*] - inFraction),reIdx)
	endif
	return, double(reIdx)/999.0D

end

function shong_histogram, inData
;--------------------------------------------------
; Making a histogram by sungryong hong in 2012
; Input: one dimensional data
; Output:
; 		out[0][1000] = histogram for the data from min to max with 1000 bins (0.1% bins)
;		out[1][1000] = cumulative version of out[0][1000]
; Example:
;		inData = randomn(seed,10000)
;		histoData = shong_histogram(inData)
;
;		x = dindgen(1000)/10.0D ;; percent unit
;		plot, x, histoData[0,*]
;		plot, x, histoData[1,*]
; Comments:
;		If min = max, the histogram is trivial. Return negative value
;--------------------------------------------------

	outData = dindgen(2,1000)


	numInData = long(n_elements(inData))
	dnumInData = double(numInData) ;; total for normalization
	numInData = numInData - 1
	maxData = double(max(inData))
	minData = double(min(inData))
	lengData = maxData - minData

;	print,"In shong_histogram.pro: Input data N =", numInData, " :: Max, Min, leng = ", maxData, minData, lengData

	if (lengData ne 0.0) then begin

		idx = 0L
		iL = 0L
		outData[*,*] = 0.0D
		for iL=0L, numInData do begin   ;; tricky with inData = max because min = 0 and max = 1 in cumulative dist.
			idx = long(1000.0*(inData[iL] - minData)/lengData)
;			print,"In shong_histogram.pro: returned idx =", idx, " :: Max, inData = ", maxData, inData[iL]
			if (idx lt 1000) then begin
				outData[0,idx] = outData[0,idx] + 1.0D
			endif else begin
				outData[0,999] = outData[0,999] + 1.0D
			endelse
		endfor

		tmpCumulCout = 0.0D
		for iL=0L, 999 do begin
			tmpCumulCout = tmpCumulCout + outData[0,iL]
			outData[1,iL] = tmpCumulCout
		endfor

		outData = outData/dnumIndata ;; normalize
	endif else begin
		outData[*,*] = -1.0D ;; negative value when max = min ... histogram is trivial
	endelse

	return, outData

end

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

