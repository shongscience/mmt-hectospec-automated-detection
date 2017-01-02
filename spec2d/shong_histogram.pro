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
