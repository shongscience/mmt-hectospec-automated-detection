pro shong_runDetection,thMin,thMax,thN,wdMin,wdMax,wdN,Ntries, lamArr, signalArr, nSigmaArr
;-----------------------------
;Sungryong Hong (1/*/2012)
;	Input: 
;		1) thresholds, widths, number of random tries, three arrays for lambda, flux, flux error
;	output:
;		2) plot the detection probabilities
;	Comments:
;	 	parameter space: (Thresholds, Widths) for given fluxSignal and fluxSigma 
;
;	Need:
;		shong_histogram.pro, shong_histogramRead.pro
;-----------------------------

	;---------------------------------------------------------------
	; EXIT if the dimenstions of lamArr and nSigmaArr are different
	if (n_elements(lamArr) eq n_elements(nSigmaArr)) then begin 

	;-------------------------
	; Initialize arrays
	;-------------------------
	arrSize = n_elements(lamArr) 
	lam= lamArr ;; assign from input lambda
	signal= signalArr
	signalSmooth= signalArr
	normalizedSignalSmooth= signalArr ;; signal / nsig ;; normalized S/N
	noise=dindgen(arrSize)	
	noiseSmooth=dindgen(arrSize)	
	noiseSigma=nSigmaArr ;; assign from input noise sigma	
	noiseSigmaSmooth=nSigmaArr ;; assign from input noise sigma	
	baseLine=dindgen(arrSize)	
	mixedSignal=dindgen(arrSize)	
	mixedSignalSmooth=dindgen(arrSize)	
	detectionFlag=dindgen(arrSize)	
	thresholdLine=dindgen(arrSize)	
	dtmparr=dindgen(arrSize)	
	detectionFlag[*] = 0.0D ;; initialize the detection array


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
	


	fluxSignal = 0.0D  ;;Initialize:  signal Gaussian flux
	fluxSigma = 1.0D  ;; Initialize:  signal Gaussian sigma
	thresholdMin = double(thMin) ;;
	thresholdMax = double(thMax) ;;
	numThresholds = long(thN) ;;
	widthMin = double(wdMin)  ;;
	widthMax = double(wdMax) ;;
	numWidths = long(wdN) ;;
	numOfTrial = long(NTries) ;;Initialize:  number of trials for statistics
	dnumOfTrial = double(numOfTrial)


	;; Make a grid in parameter space (thershold, width)
	thresholds=thresholdMin + dindgen(numThresholds)/(double(numThresholds)) * (thresholdMax - thresholdMin)
	widths=widthMin + dindgen(numWidths)/(double(numWidths)) * (widthMax - widthMin)
	probResults = dindgen(numThresholds,numWidths)
	probResults[*,*] = 0.0D
	detectResults = dindgen(numThresholds,numWidths)
	detectResults[*,*] = 0.0D
	detectResultsLam = dindgen(arrSize,numThresholds,numWidths)
	detectResultsLam[*,*,*] = 0.0D
	detectLam = dindgen(arrSize)
	detectLamWeighted = dindgen(arrSize)

	iThres = 0L ;; looping parameter space
	iWidths = 0L
	iThresMax = numThresholds -1
	iWidthsMax = numWidths -1 



	;; Assign a simulated profile -------------------------
	fluxSignal = 0.0D ;; zero flux input
	fluxSigma = 5.0D
	
    ;dtmparr=-1.0*(lam-50.)*(lam-50.)/(2.0*fluxSigma*fluxSigma)
    ;signal=fluxSignal*exp(dtmparr)	


	;; Assign a normal noise (1) assign sigma value
	;noiseSigma = 5.0D ;; assign from input argument... So, an obsolete line
	;; -----------------------------------------------------




	;-------------------------------------------------------------
	; Subtract the baseline 
	baseLine = shong_linearbaseline(lam,signal,noiseSigma)
	signalSmooth = signalSmooth - baseLine 

	;; set iteration and related quantities
	iterL = 0L
	iterMax = numOfTrial - 1L	
	iDetect = 0L
	iDetectMax = 0L
	dtmp = 0.0D
	dtmpNoiseSigmaMax = max(noiseSigma)
	goplot=0

	;-------------------------
	; Plot the related files 
	;-------------------------
   	window,10,retain=2,xsize=1800,ysize=1000
    !p.charsize=2.0
    !p.multi=[10,2,3]

	wset,10

	;--------------------------------
	; Create smooth kernels
	;; SMOOTH 
	smoothSigma = 1.0D
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
	; Plot the results (1/3)
	;
				
				noise=randomn(seed, arrSize)*noiseSigma
				noiseSmooth = convol(noise,gaussian_ker,/center)

			    djs_plot, lam, signal, xtitle='lambda', ytitle='Signal + Sigma + BaseLine'
			    djs_oplot, lam, noiseSigma, color='red'
			    djs_oplot, lam, baseLine, color='blue', thick=3
				djs_xyouts, lam[arrSize/20],(0.8*max(signal,/nan)),'Input Signal + Sigma + BaseLine', color='green'

			    djs_plot, lam, signalSmooth, xtitle='lambda', ytitle='[Smoothed] BaseLine Subtracted Signal + Sigma'
			    djs_oplot, lam, noiseSigmaSmooth, color='red'
			    djs_oplot, lam, noiseSmooth, color='blue'
				djs_xyouts, lam[arrSize/20],(0.9*max(signalSmooth,/nan)),'Smooth + BaseLineSubtracted', color='green'

			    djs_plot, lam, normalizedSignalSmooth, xtitle='lambda', ytitle='Signal / Sigma'
				thresholdLine[*] = min(thresholds,/nan)
			    djs_oplot, lam, thresholdLine, color='yellow'
				thresholdLine[*] = max(thresholds,/nan)
			    djs_oplot, lam, thresholdLine, color='yellow'
				djs_xyouts, lam[arrSize/20],(0.9*max(normalizedSignalSmooth,/nan)),'Normalized', color='green'





	;---------------------------------
	; Get False Detection rates
	; Description:
	;	1) Make random noises (iterMax times) 
	;	2) for each noise, probe the detectin on the (thes, widths) space
	; Final output is probResults[*,*]
	for iterL = 0L, iterMax do begin ;;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		
		;; Assign a normal noise (2) make a random noise
		noise=randomn(seed, arrSize)*noiseSigma

		;; Mix the singal and noise
		;mixedSignal = signal + noise ;;; Only for theoretical gaussian model
		mixedSignal = noise
	 	mixedSignalSmooth = convol(mixedSignal,gaussian_ker,/center)
		;; normalize the signa to S/N
		mixedSignal = mixedSignal/noiseSigma
		mixedSignalSmooth = mixedSignalSmooth/noiseSigmaSmooth


		;; Signal is done!
		;; Now calculate the fraction for each (threshold, width)
		;; --- Looping the parameter spaces (threshold, widths)  
		iThres = 0L
		iWidths = 0L
		iLam = 0L
		tmpIdxLam = 0L
		iLamMax = 0L

		for iThres=0L, iThresMax do begin ;;*******************************************


		;; threshold cut and count the survived points' width 
			detectionFlag[*] = 0.0D
			thresholdLine[*] = thresholds[iThres]
			detectedPixels = where(mixedSignalSmooth ge thresholds[iThres],count)
			;detectionFlag[detectedPixels] = 1.0D
			;numDetected = n_elements(detectedPixels)
			detections = findAllCluster(detectedPixels)
	
			;Debugging logs
			;print, numDetected
			;print, size(detectedPixels)
			;print, detectedPixels


			;------------
			;Each detection
			;-------------

			if (long((size(detections))[0]) eq 1) then begin ;; single detection, so 1D array
				iDetectMax = 1
			endif else begin                                 ;; multiple detection, so 2D array
				iDetectMax = long((size(detections))[2])
			endelse
			iDetectMax = iDetectMax -1
			iDetect=0L			

			tmpMaxIdx = 0L  ;; largest detection cluster
			tmpMaxWidth = 0.0D ;; largest detection cluster

			if (iDetectMax gt 0) then begin ;; if multiple detection, use 2D data
				;print, "Detections Matrix = ", (size(detections))[1],  (size(detections))[2]
				;print, "size(detections) = ", size(detections)
				;print, detections,"iDetectMax = ",iDetectMax

				tmpMaxWidth = double(detections[1,iDetect]);; save the maximum cluster: width
				tmpMaxIdx = detections[0,iDetect]          ;; save the maximum cluster: index
				for iDetect=0L, iDetectMax do begin
	
					;print, "iDetect = ",iDetect,"/",iDetectMax
					dtmp = double(detections[1,iDetect])
					tmpIdxLam = detections[0,iDetect]
					if (dtmp gt tmpMaxWidth) then begin
						tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
						tmpMaxWidth = dtmp    ;; save the maximum cluster: Width
					endif
					iLamMax = long(dtmp)
					iLamMax = iLamMax - 1
					;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

					;for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%  ProbLam is removed
					;	if (dtmp ge widths[iWidths]) then begin
					;		for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
					;			probResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
					;		endfor
					;	endif
					;endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probResultsLam record all the detected clusters at lambda
				endfor


				for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%% probResults record only once for the largest cluster
						if (tmpMaxWidth ge widths[iWidths]) then begin
							probResults[iThres,iWidths] += 1.0D
						endif
				endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			endif else begin               ;; if single detection, use 1D data
				;print, "iDetect = ",iDetect,"/",iDetectMax
				dtmp = double(detections[1])
				tmpIdxLam = detections[0]

				tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
				tmpMaxWidth = dtmp    ;; save the maximum cluster: Width

				;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

				if (tmpIdxLam ne (-1)) then begin   ;; (-1) is a flag for null detection
					for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%
						if (dtmp ge widths[iWidths]) then begin
							probResults[iThres,iWidths] += 1.0D
							;for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
							;	probResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
							;endfor
						endif
					endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				endif else begin
					;print,"Null detction"
				endelse
			endelse


			;----------------
			;Each detection Ends
			;-------------------


		;;;;interupt
		;voidstring = ''
		;read, voidstring, prompt='<<<<<<<<<<<<<<<<enter for the next results>>>>>>>>>>>>>>>>>>>>'

		endfor ;; *********************************************************************

	endfor ;;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	probResults = probResults/dnumOfTrial
	probResults = 1.0D - probResults
	;End of Getting False Detection: output = probResults[][]
	;-------------------------





	;-------------------------------------------------------------------------------------------------------------
	; Plot the results (2/3)
	; 	False detection rates

				contour, probResults,thresholds,widths, nlevels=20,/fill, xtitle='Thresholds', ytitle='Widths'
				djs_xyouts,(0.4*max(thresholds,/nan)),(0.8*max(widths,/nan)),'Detection Probability', color='blue'
				;contour, detectResults,thresholds,widths, nlevels=20,/fill, xtitle='Thresholds', ytitle='Widths'







	;-------------------------------------------------------------------------------------------------------------
	; Get Detection Rate using the input Signal
	; Input: normalizedSignalSmooth
	; Procedure: 
	;	1) in the para space, (thres, widths), probe the detection
	
		iThres = 0L
		iWidths = 0L
		iLam = 0L
		tmpIdxLam = 0L
		iLamMax = 0L
		detectResults[*,*] = 0.0D
		detectResultsLam[*,*,*] = 0.0D
		detectLam[*] = 0.0D
		detectLamWeighted[*] = 0.0D

		goplot=0
		for iThres=0L, iThresMax do begin ;;*******************************************


		;; threshold cut and count the survived points' width 
			detectionFlag[*] = 0.0D
			thresholdLine[*] = thresholds[iThres]
			detectedPixels = where(normalizedSignalSmooth ge thresholds[iThres],count)
			;detectionFlag[detectedPixels] = 1.0D
			;numDetected = n_elements(detectedPixels)
			detections = findAllCluster(detectedPixels)
	
			;Debugging logs
			;print, numDetected
			;print, size(detectedPixels)
			;print, detectedPixels


			;------------
			;Each detection
			;-------------

			if (long((size(detections))[0]) eq 1) then begin ;; single detection, so 1D array
				iDetectMax = 1
			endif else begin                                 ;; multiple detection, so 2D array
				iDetectMax = long((size(detections))[2])
			endelse
			iDetectMax = iDetectMax -1
			iDetect=0L			

			tmpMaxIdx = 0L  ;; largest detection cluster
			tmpMaxWidth = 0.0D ;; largest detection cluster

			if (iDetectMax gt 0) then begin ;; if multiple detection, use 2D data
				;print, "Detections Matrix = ", (size(detections))[1],  (size(detections))[2]
				;print, "size(detections) = ", size(detections)
				;print, detections,"iDetectMax = ",iDetectMax

				tmpMaxWidth = double(detections[1,iDetect]);; save the maximum cluster: width
				tmpMaxIdx = detections[0,iDetect]          ;; save the maximum cluster: index
				for iDetect=0L, iDetectMax do begin
	
					;print, "iDetect = ",iDetect,"/",iDetectMax
					dtmp = double(detections[1,iDetect])
					tmpIdxLam = detections[0,iDetect]
					if (dtmp gt tmpMaxWidth) then begin
						tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
						tmpMaxWidth = dtmp    ;; save the maximum cluster: Width
					endif
					iLamMax = long(dtmp)
					iLamMax = iLamMax - 1
					;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

					for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%  ProbLam is removed
						if (dtmp ge widths[iWidths]) then begin
							for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
								detectResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
							endfor
						endif
					endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probResultsLam record all the detected clusters at lambda
				endfor


				for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%% probResults record only once for the largest cluster
						if (tmpMaxWidth ge widths[iWidths]) then begin
							detectResults[iThres,iWidths] += 1.0D
						endif
				endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			endif else begin               ;; if single detection, use 1D data
				;print, "iDetect = ",iDetect,"/",iDetectMax
				dtmp = double(detections[1])
				tmpIdxLam = detections[0]

				tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
				tmpMaxWidth = dtmp    ;; save the maximum cluster: Width

				;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

				if (tmpIdxLam ne (-1)) then begin   ;; (-1) is a flag for null detection
					for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%
						if (dtmp ge widths[iWidths]) then begin
							detectResults[iThres,iWidths] += 1.0D
							for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
								detectResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
							endfor
						endif
					endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				endif else begin
					;print,"Null detction"
				endelse
			endelse
			;----------------
			;Each detection Ends
			;-------------------

          if (goplot) then begin
				print, normalizedSignalSmooth
                print,"Detection Size = ", size(detections)
				print,"Maximums: signal, signalSmooth, norm.Sig.Smooth =", max(signal),max(signalSmooth),max(normalizedSignalSmooth)
                detectionFlag[detectedPixels] = 1.0D
				detectLam = total(total(detectResultsLam,2),2)

				noise=randomn(seed, arrSize)*noiseSigma
				noiseSmooth = convol(noise,gaussian_ker,/center)

			    djs_plot, lam, signal, xtitle='lambda', ytitle='Signal and Sigma'
			    djs_oplot, lam, noiseSigma, color='red'
			    djs_oplot, lam, noise, color='blue'
				djs_xyouts, (arrSize/20.),(0.9*max(signal)),'Input Signal + Sigma', color='green'

			    djs_plot, lam, signalSmooth, xtitle='lambda', ytitle='Signal and Sigma'
			    djs_oplot, lam, noiseSigmaSmooth, color='red'
			    djs_oplot, lam, noiseSmooth, color='blue'
				djs_xyouts, (arrSize/20.),(0.9*max(signalSmooth)),'SMOOTH', color='green'

			    djs_plot, lam, normalizedSignalSmooth, xtitle='lambda', ytitle='Signal and Sigma'
			    djs_oplot, lam, thresholdLine, color='red'
				djs_xyouts, (arrSize/20.),(0.7*max(normalizedSignalSmooth,/nan)),'Normalized', color='green'

                djs_plot, lam, detectionFlag, xtitle='x: PIXEL INDEX', ytitle='Flux', color='blue', yrange=[0.0,1.5]
                djs_xyouts, (arrSize/20.),1.3,'Width', color='green'
                for iDetect=0L, iDetectMax do begin
                    djs_xyouts, lam[detections[0,iDetect]],1.1,'V', color='blue'
                endfor
                djs_xyouts, lam[tmpMaxIdx],1.1,'V', color='red'

                djs_plot, lam, detectResultsLam[*,2,0], $
						xtitle='Thresholds', ytitle='Prob(th,wd)'
                for iWidthsTmp=1L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%
                    djs_oplot, lam, detectResultsLam[*,2,iWidthsTmp],color=iWidthsTmp
                endfor
                djs_plot, lam, detectLam/max(detectLam), xtitle='lam', ytitle='Prob(th,wd)', yrange=[0.0,1.1]
            endif


		;;;;interupt
		;voidstring = ''
		endfor ;; *********************************************************************
		


		detectLam = total(total(detectResultsLam,2),2)
		iL = 0L
		iMaxL = arrSize -1 
		for iL=0, iMaxL do begin 
			detectResultsLam[iL,*,*] = detectResultsLam[iL,*,*] * probResults[*,*]
		endfor 
		detectLamWeighted = total(total(detectResultsLam,2),2)

	; End of Getting Detection Rate
	;-------------------------------------------------------------------------------------------------------------











	;---------------------------------
	; Get Detection Limit
	; Description:
	;	1) Make random noises (iterMax times) 
	;	2) for each noise, probe the detectin on the (thes, widths) space
	;	3) Calculate the two detection rates for each random realization
	;		3-1) Prob(lam, N) without weighting
	;		3-2) Prob(lam, N) with weighting with ProbResults[][] 
	; Comments:
	;	- This is a redundant step of the first. 
	;		If we want to calculate all of the results in the first step, 
	;		the memory consumption is too high. 
	;		Hence, at the first step, we get the probResults[][]. 
	;		In this step, we get the two lamda-dependent probs below 
	; Final output is probResultsLamWeighted[][], probResults


	probResultsLam = dindgen(arrSize,numThresholds,numWidths)
	probResultsLamNoWeight = dindgen(arrSize,numOfTrial)
	probResultsLamWeight = dindgen(arrSize, numOfTrial)
	probResultsLam[*,*,*] = 0.0D
	probResultsLamNoWeight[*,*] = 0.0D
	probResultsLamWeight[*,*] = 0.0D

	for iterL = 0L, iterMax do begin ;;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		
		;; Assign a normal noise (2) make a random noise
		noise=randomn(seed, arrSize)*noiseSigma

		;; Mix the singal and noise
		;mixedSignal = signal + noise ;;; Only for theoretical gaussian model
		mixedSignal = noise
	 	mixedSignalSmooth = convol(mixedSignal,gaussian_ker,/center)
		;; normalize the signa to S/N
		mixedSignal = mixedSignal/noiseSigma
		mixedSignalSmooth = mixedSignalSmooth/noiseSigmaSmooth


		;; Signal is done!
		;; Now calculate the fraction for each (threshold, width)
		;; --- Looping the parameter spaces (threshold, widths)  
		iThres = 0L
		iWidths = 0L
		iLam = 0L
		tmpIdxLam = 0L
		iLamMax = 0L

		for iThres=0L, iThresMax do begin ;;*******************************************


		;; threshold cut and count the survived points' width 
			detectionFlag[*] = 0.0D
			thresholdLine[*] = thresholds[iThres]
			detectedPixels = where(mixedSignalSmooth ge thresholds[iThres],count)
			;detectionFlag[detectedPixels] = 1.0D
			;numDetected = n_elements(detectedPixels)
			detections = findAllCluster(detectedPixels)
	
			;Debugging logs
			;print, numDetected
			;print, size(detectedPixels)
			;print, detectedPixels


			;------------
			;Each detection
			;-------------

			if (long((size(detections))[0]) eq 1) then begin ;; single detection, so 1D array
				iDetectMax = 1
			endif else begin                                 ;; multiple detection, so 2D array
				iDetectMax = long((size(detections))[2])
			endelse
			iDetectMax = iDetectMax -1
			iDetect=0L			

			tmpMaxIdx = 0L  ;; largest detection cluster
			tmpMaxWidth = 0.0D ;; largest detection cluster

			if (iDetectMax gt 0) then begin ;; if multiple detection, use 2D data
				;print, "Detections Matrix = ", (size(detections))[1],  (size(detections))[2]
				;print, "size(detections) = ", size(detections)
				;print, detections,"iDetectMax = ",iDetectMax

				tmpMaxWidth = double(detections[1,iDetect]);; save the maximum cluster: width
				tmpMaxIdx = detections[0,iDetect]          ;; save the maximum cluster: index
				for iDetect=0L, iDetectMax do begin
	
					;print, "iDetect = ",iDetect,"/",iDetectMax
					dtmp = double(detections[1,iDetect])
					tmpIdxLam = detections[0,iDetect]
					if (dtmp gt tmpMaxWidth) then begin
						tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
						tmpMaxWidth = dtmp    ;; save the maximum cluster: Width
					endif
					iLamMax = long(dtmp)
					iLamMax = iLamMax - 1
					;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

					for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%  ProbLam is removed
						if (dtmp ge widths[iWidths]) then begin
							for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
								probResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
							endfor
						endif
					endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probResultsLam record all the detected clusters at lambda
				endfor


				;for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%% probResults record only once for the largest cluster
				;		if (tmpMaxWidth ge widths[iWidths]) then begin
				;			probResults[iThres,iWidths] += 1.0D
				;		endif
				;endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			endif else begin               ;; if single detection, use 1D data
				;print, "iDetect = ",iDetect,"/",iDetectMax
				dtmp = double(detections[1])
				tmpIdxLam = detections[0]

				tmpMaxIdx = tmpIdxLam ;; save the maximum cluster: index
				tmpMaxWidth = dtmp    ;; save the maximum cluster: Width

				;print, " Start Lambda = ",lam[tmpIdxLam]," Width = ", dtmp

				if (tmpIdxLam ne (-1)) then begin   ;; (-1) is a flag for null detection
					for iWidths=0L, iWidthsMax do begin ;;%%%%%%%%%%%%%%%%%
						if (dtmp ge widths[iWidths]) then begin
							;probResults[iThres,iWidths] += 1.0D
							for iLam=0L, iLamMax do begin ;;%%%%%%%%%%%%%%%%%
								probResultsLam[tmpIdxLam+iLam,iThres,iWidths] += 1.0D
								;print,"Prob(", thresholds[iThres], widths[iWidths], lam[tmpIdxLam+iLam],")", $ 
								;		iLam,"/", iLamMax," Delta Width = ", dtmp
							endfor
						endif
					endfor;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				endif else begin
					;print,"Null detction"
				endelse
			endelse


			;----------------
			;Each detection Ends
			;-------------------


		;;;;interupt
		;voidstring = ''
		;read, voidstring, prompt='<<<<<<<<<<<<<<<<enter for the next results>>>>>>>>>>>>>>>>>>>>'

		endfor ;; *********************************************************************



		probResultsLamNoWeight[*,iterL] = total(total(probResultsLam,2),2)
		iL = 0L
		iMaxL = arrSize -1 
		for iL=0, iMaxL do begin 
			probResultsLam[iL,*,*] = probResultsLam[iL,*,*] * probResults[*,*]
		endfor 
		probResultsLamWeight[*,iterL] = total(total(probResultsLam,2),2)

		probResultsLam[*,*,*] = 0.0D ;;initialize
	endfor ;;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	;summarize the results
	probResultsLamNoWeight950 = dindgen(arrSize)
	probResultsLamNoWeight990 = dindgen(arrSize)
	probResultsLamNoWeight999 = dindgen(arrSize)
	probResultsLamWeight950 = dindgen(arrSize)
	probResultsLamWeight990 = dindgen(arrSize)
	probResultsLamWeight999 = dindgen(arrSize)

	tmparr = dindgen(arrSize)
	iL = 0L
	iMaxL = arrSize -1 
	outHisto = dindgen(2,1000)
	tmpMin = 0.0D
	tmpMax = 0.0D

	for iL=0, iMaxL do begin 
		outHisto[*,*] = 0.0D
		tmparr = probResultsLamWeight[iL,*]
		tmpMin = min(tmparr)
		tmpMax = max(tmparr)
		outHisto = shong_histogram(tmparr)
		probResultsLamWeight990[iL] = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.99)
		probResultsLamWeight999[iL] = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.999)

		outHisto[*,*] = 0.0D
		outHisto = shong_histogram(probResultsLamNoWeight[iL,*])
		tmpMin = min(probResultsLamNoWeight[iL,*])
		tmpMax = max(probResultsLamNoWeight[iL,*])
		probResultsLamNoWeight990[iL] = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.99)
		probResultsLamNoWeight999[iL] = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.999)
	endfor 
	
	tmpTotalW = reform(probResultsLamWeight, n_elements(probResultsLamWeight))
	outHisto = shong_histogram(tmpTotalW)
	tmpMin = min(tmpTotalW)
	tmpMax = max(tmpTotalW)
	tmpW990 = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.99)
	tmpW999 = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.999)

	tmpTotalNoW = reform(probResultsLamNoWeight, n_elements(probResultsLamNoWeight))
	outHisto = shong_histogram(tmpTotalNoW)
	tmpMin = min(tmpTotalNoW)
	tmpMax = max(tmpTotalNoW)
	tmpNoW990 = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.99)
	tmpNoW999 = tmpMin+(tmpMax-tmpMin)*shong_histogramRead(outHisto,0.999)





	;
	;End of Getting False Detection: output = probResultsLamWeight and NoWeight
	;-------------------------








	;-------------------------------------------------------------------------------------------------------------
	; Plot the results (3/3)
	;
	;	Detections and their limits



                djs_plot, lam, detectLam, xtitle='lam',$
					ytitle='Detection Rates', yrange=[0.0,2.0*tmpNoW999], color='magenta'
				djs_oplot, lam, probResultsLamNoWeight990, color='green'
				djs_oplot, lam, probResultsLamNoWeight999, color='green'
				thresholdLine[*] = tmpNoW990
				djs_oplot,lam, thresholdLine, color='green', thick=3
				thresholdLine[*] = tmpNoW999
				djs_oplot,lam, thresholdLine, color='green', thick=3

                djs_plot, lam, detectLamWeighted, color='yellow', $
					xtitle='lam', ytitle='Detection Rates (Weighted)', yrange=[0.0,2.0*tmpW999]
				djs_oplot, lam, probResultsLamWeight990, color='green'
				djs_oplot, lam, probResultsLamWeight999, color='green'
				thresholdLine[*] = tmpW990
				djs_oplot,lam, thresholdLine, color='green', thick=3
				thresholdLine[*] = tmpW999
				djs_oplot,lam, thresholdLine, color='green', thick=3



;				djs_plot, probResultsLamNoWeight[100,*]
;				djs_plot, probResultsLamWeight[100,*]

;				outHisto = shong_histogram(probResultsLamNoWeight[100,*])
;				print, min(probResultsLamNoWeight[100,*])
;				print, max(probResultsLamNoWeight[100,*])
;				djs_plot, outHisto[0,*]
;				djs_plot, outHisto[1,*]

			;	print,"<<In shong_runDetection.pro>>"
			;	print,"maxes: signal, signalSmooth, normalizedSignalSmooth = ", $
			;			max(signal,/nan), max(signalSmooth,/nan), max(normalizedSignalSmooth,/nan)











	endif else begin ;-----------------------------------------------------------------------------------------------------
	;EXIT Control
		print, "ERROR: The lambda and nSigma arrays are not compatible in dimension: Check size(arr)" 
	endelse
	;----------------------------------------------------------------------------------------------------------------------
end

;-----------------------------------------------------------------------------------
; from the index data by where(DetectionArray), find the largest chunk of detection
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

