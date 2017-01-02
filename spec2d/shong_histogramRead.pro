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
