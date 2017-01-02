;+
; NAME:
;   shong_readcat
;
;	******* Modified version of hs_readcat.pro
;
; PURPOSE:
;   Read the catalog file in order to match them to the observation
;
; CALLING SEQUENCE:
;  hs_readcat, catname, mapfile, columnnames, colformat
;
; INPUTS:
;    catname - name of input catalog to be matched
;    mapfile - name of HS map file to be matched
;    columnnames - [ncol] - string array containing the names of column
;                           entries.  Note that should match the standard
;                           plugcat standards for anything you want to get
;                           included in the plugmap.  These are listed in the
;                           documentation.
;    colformat - [ncol] - Format string for each columns
;                         A - string
;                         D - Double
;                         F - Float
;                         I - integer
;                         L - Long
;
; OPTIONAL KEYWORDS:
;     comment  - comment string (default to '#')
;     skyfile  - file containing catalog of sky locations
;     skycolumnnames - same as column names for skies
;     skycolformat - same as column names for skies
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;   plugcat/mapfile_XXXX files which will server as the plugmaps
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   March 2005 - Written by R Cool - UofA
;-
;------------------------------------------------------------------------------
PRO shong_readcat, catname, mapfile, columnnames, colformat, comment=comment, $
                skyfile=skyfile, skycolumnnames=skycolumnnames, $
                skycolformat=skycolformat, nomatch=nomatch
      
   IF NOT keyword_set(comment) THEN comment = '#'
   colname = strupcase(columnnames)
   file = rc_readfile(catname, 0)
   
   ;;Find the lines that does not have a comment string in front of it
   k = where(strmid(file, 0, 1) NE comment AND strtrim(file) NE '', ct)
   file = file(k)
   
   delvarx, output
   
	;;error finding
	print,"LOG:: in hs_readcat.pro -- colname: ",colname  
	print,"LOG:: in hs_readcat.pro -- size(colname): ",size(colname)  
	;print,"LOG:: in hs_readcat.pro -- colname: ",colformat
	;print,"LOG:: in hs_readcat.pro -- size(colname): ",size(colformat)  
	;print,"LOG:: in hs_readcat.pro -- file: ",file  
	;print,"LOG:: in hs_readcat.pro -- size(file): ",size(file)  

   create_struct, output, 'file_input', colname, colformat, $
      dimen=ct
 

   FOR iline = 0l, ct -1l DO BEGIN
      
		;print,"LOG:: in hs_readcat.pro -- iline: ",iline
      split = strsplit(file(iline),  /extract)

		;;error finding
		;print,"LOG:: in hs_readcat.pro -- split: ",split
		;print,"LOG:: in hs_readcat.pro -- size(split): ",size(split)
		;print,"LOG:: in hs_readcat.pro -- n_elements(split) N_elements(colname): ",n_elements(split),n_elements(colname)

      IF n_elements(split) NE n_elements(colname) THEN BEGIN
         splog, 'WARNING - line ' + string(k(iline)) + ' does not match'
         splog, 'the input column names'
         
         return
      endIF
      
      FOR icol = 0, n_elements(colname) -1 DO BEGIN
         output[iline].(icol) = split(icol)
      ENDFOR
   ENDFOR
      
   
   ;;Now, do the sky
   IF keyword_set(skyfile) THEN BEGIN
      
      scolname = strupcase(skycolumnnames)
      file = rc_readfile(skyfile, 0)
      
      ;;Find the lines that does not have a comment string in front of it
      k = where(strmid(file, 0, 1) NE comment AND strtrim(file) NE '', ct)
      file = file(k)
      
      delvarx, skyoutput
  
      create_struct, skyoutput, 'sky_input', scolname, skycolformat, $
         dimen=ct
      
      FOR iline = 0l, ct -1l DO BEGIN
         
         split = strsplit(file(iline), /extract)
         IF n_elements(split) NE n_elements(scolname) THEN BEGIN
            splog, 'WARNING - line ' + string(k(iline)) + 'does not match'
            splog, 'the input scolumn names'
            
            return
         endIF
         
         FOR icol = 0, n_elements(scolname) -1 DO BEGIN
            skyoutput[iline].(icol) = split(icol)
         ENDFOR
         
      ENDFOR
   endIF
   
   ;;Now create the needed vectors
   tags = tag_names(output)
   ncol = n_elements(output)

	;; error finding
	print,"LOG:: in hs_readcat.pro -- TAGS: ",tags  
	print,"LOG:: in hs_readcat.pro -- size(TAGS): ",size(tags) 
	print,"LOG:: in hs_readcat.pro -- ncol: ",ncol
	 
  
   k = where(tags EQ 'RA', ct)
   IF ct NE 0 THEN ra = output.ra *15     ELSE BEGIN
      splog, 'RA not found - this is bad'
      ra = replicate(0.0, ncol)
   ENDELSE
	print,"LOG:: in hs_readcat.pro -- RA index: ",k

 
   k = where(tags EQ 'DEC', ct)
   IF ct NE 0 THEN dec = output.dec ELSE BEGIN
      splog, 'DEC not found - this is bad'
      dec = replicate(0.0, ncol)
   ENDELSE
	print,"LOG:: in hs_readcat.pro -- DEC index: ",k
   
   k = where(tags EQ 'TYPE', ct)
   IF ct NE 0 THEN type = output.type ELSE BEGIN
      splog, 'TYPE not found' 
      splog, 'Setting all types to TARGET'
      type = replicate('TARGET', ncol)
   ENDELSE
	print,"LOG:: in hs_readcat.pro -- TYPE index: ",k
   
   k = where(tags EQ 'RANK', ct)
   IF ct NE 0 THEN rank = output.rank ELSE rank = replicate(0, ncol)
	print,"LOG:: in hs_readcat.pro -- RANK index: ",k
   
   k = where(tags EQ 'RMAG', ct)
   IF ct NE 0 THEN rmag = output.rmag ELSE BEGIN
      splog, 'RMAG is not set - setting to zero'
      rmag = replicate(0, ncol)
   ENDELSE
   
   k = where(tags EQ 'RAPMAG', ct)
   IF ct NE 0 THEN rapmag = output.rapmag ELSE BEGIN
      splog, 'RAPMAG not set'
      k = where(rmag GT 0, ct)
      IF ct EQ 0 THEN BEGIN
         splog, 'Setting to zero'
         rapmag = replicate(0, ncol)
      endIF
      IF ct NE 0 THEN BEGIN
         splog, 'Setting to RMAG'
         rapmag = rmag
      endIF
      
      
   ENDELSE
   
   K  = WHERE(TAGS EQ 'BCODE', CT)
   IF CT NE 0 THEN BCODE = OUTPUT.BCODE ELSE BEGIN
      splog, 'bcode not set'
      
      BCODE = REPLICATE(0, NCOL)
   ENDELSE
   
   K = WHERE(TAGS EQ 'ICODE', CT)
   IF CT NE 0 THEN ICODE = OUTPUT.ICODE ELSE BEGIN
      splog, 'icode not set'
      ICODE = REPLICATE(0, NCOL)
   ENDELSE
   
   k = where(tags EQ 'RCODE', ct)
   IF ct NE 0 THEN rcode = output.rcode ELSE BEGIN
      rcode = replicate(0, ncol)
      splog, 'rcode not set'
   ENDELSE
   
   k = where(tags EQ 'NEWICODE', ct)
   IF ct NE 0 THEN newicode = output.newicode ELSE BEGIN
      newicode = replicate(0, ncol)
      splog, 'newicode not set'
   ENDELSE
   
   
   
   
   IF keyword_set(skyfile) THEN BEGIN
      
      tags = tag_names(skyoutput)
      
      k = where(tags EQ 'RA', ct)
      IF ct NE 0 THEN ras = skyoutput.ra *15 ELSE BEGIN
         splog, 'RA not set in Sky file - THIS IS BAD!'
      ENDELSE
      
      K = WHERE(TAGS EQ 'DEC', CT)
      IF CT NE 0 THEN DECS = skyOUTPUT.DEC ELSE BEGIN
         SPLOG, 'DEC NOT SET IN SKY FILE - THIS IS BAD!'
      ENDELSE
      
   ENDIF
   
   
   
   ;;Setup the output and format the names
      
   if direxist('plugdir') eq 0L then spawn, 'mkdir plugdir'
   junk=strsplit(mapfile, ' ', length=ct)
   outname=strmid(mapfile, 0,ct-9)
   
   
   ;;Read the map  
   readcol, mapfile, junk, junk1, target1, ra1, dec1, junk, fiber1, $
      junk,junk, format='L,L,A,A,A,A,A,A,A'
   k = where(strmatch(target1, 'broken*') EQ 1l, ct1)
   IF ct1 GT 0 THEN target1(k) = 'unused'
   
   ;;Now do the matching on the targets
   k = where(strmatch(target1,'sky*') eq 0l AND $
             strmatch(target1,'unused*') eq 0l)

   fiber1 = fiber1(k)
   target1=target1(k)
   ra1=ra1(k)
   dec1=dec1(k)
   
   for i= 0l, n_elements(ra1) -1l do begin
      ra1(i) = hms2dec(ra1(i))*15
      dec1(i) = hms2dec(dec1(i))
   endfor
   
    spherematch, ra, dec, ra1, dec1, 1.0/3600., match1, match2,dist
    
      
    ra2 = dindgen(300)*0.0-1.0
    dec2 = dindgen(300)*0.0-1.0
    target2 = strarr(300)
    dmin2= lindgen(300)*0.0
    rmag2 = dindgen(300)*0.0-1.0
    rapmag2 = dindgen(300)*0.0-1.0
    bcode2 = strarr(300)
    icode2 = lindgen(300)*0+1
    newicode2 = lindgen(300)*0-1
    rcode2 = lindgen(300)*0-1
    target = strarr(300)
    fiber = lindgen(300)*0-1
    ct = n_elements(match1)
    ra2(0:ct-1) = ra(match1)
    dec2(0:ct-1) = dec(match1)
    fiber(0:ct-1) = fiber1(match2)
    rmag2(0:ct-1) = rmag(match1)
    rapmag2(0:ct-1) = rapmag(match1)
    bcode2(0:ct-1) = bcode(match1)
    icode2(0:ct-1) = icode(match1)
    newicode2(0:ct-1) = newicode(match1)
    rcode2(0:ct-1) = rcode(match1)
    dmin2(0:ct-1) = 0.0
    target(0:ct-1) = 'target'
    
   IF keyword_set(ras) THEN begin
      ;;Now match the skies
      
      spherematch, ra1, dec1, ras, decs, 1.0/3600., match1, match2, dist
      
      if match1(0) ne -1 then begin
         
         cts = n_elements(match1)
         ra2(ct:cts+ct-1) = ra(match2)
         dec2(ct:cts+ct-1) = dec(match2)
         fiber(ct:cts+ct-1) = fiber1(match1)
         rmag2(ct:cts+ct-1) = 0
         rapmag2(ct:cts+ct-1) = 0
         bcode2(ct:cts+ct-1) = 0
         icode2(ct:cts+ct-1) = 0
         newicode2(ct:cts+ct-1) = 0
         rcode2(ct:cts+ct-1) = 0
         dmin2(ct:cts+ct-1) = 0
         target(ct:cts+ct-1) = 'sdsssky'
         
      endif else cts =0 
      
   ENDIF ELSE cts = 0

   
   ;;Now find the unused
   ;;Read the map
      
   readcol, mapfile, junk, junk1, target1, ra1, dec1, junk, fiber1, $
       junk,junk, format='L,A,A,A,A,A,A,A,A'
    
    k = where(strmatch(target1, 'broken*') EQ 1l, ct1)  
    IF ct1 GT 0 THEN target1(k) = 'unused'           
    

    k = where(strmatch(target1, 'unused*'))
    
    if k(0) ne -1 then begin
       ctu = n_elements(k)
       
       ra2(cts+ct:cts+ct+ctu-1) = 0.0
       dec2(cts+ct:cts+ct+ctu-1) = 0.0
       fiber(cts+ct:cts+ct+ctu-1) = fiber1(k)
       rmag2(cts+ct:cts+ct+ctu-1) = 0
       rapmag2(cts+ct:cts+ct+ctu-1) = 0
       bcode2(cts+ct:cts+ct+ctu-1) = 0
       icode2(cts+ct:cts+ct+ctu-1) = 0
       newicode2(cts+ct:cts+ct+ctu-1) = 0
       rcode2(cts+ct:cts+ct+ctu-1) = 0
       dmin2(cts+ct:cts+ct+ctu-1) = 0
       target(cts+ct:cts+ct+ctu-1) = 'unused'
       
    endif else ctu = 0


    ;;Now the normal skies
    ;;Read the map
    
   readcol, mapfile, junk, junk1, target1, ra1, dec1, junk, fiber1, $
       junk, format='L,A,A,A,A,A,A,A,A'
    
    k = where(strmatch(target1,'sky*'))
    fiber1 = fiber1(k)
    target1=target1(k)
    ra1=ra1(k)
    dec1=dec1(k)
    
    k = where(strmatch(target1,'sky*'))
    if k(0) ne -1 then begin
       cts1 = n_elements(k)
       ra2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0.0
       dec2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0.0
       fiber(ctu+cts+ct:ctu+ct+cts-1+cts1) = fiber1(k)
       rmag2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       rapmag2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       bcode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       icode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       newicode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       rcode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       dmin2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
       target(ctu+cts+ct:ctu+ct+cts-1+cts1) = 'sky'
    endif else ctu = 0
    

    sort1 = sort(fiber)
    fiber = fiber(sort1)
    dmin = dmin2(sort1)
    ra = ra2(sort1)/15.0
    dec = dec2(sort1)
    rmag = rmag2(sort1)
    rapmag = rapmag2(sort1)
    bcode = bcode2(sort1)
    icode = icode2(sort1)
    newicode = newicode2(sort1)
    rcode = rcode2(sort1)
    target = target(sort1)

    ;;This is a hack incase of a nonmatch
    if keyword_set(nomatch) then begin
   
        ;;Read the map  
   		readcol, mapfile, junk, junk1, target1, ra1, dec1, junk, fiber1, $
                 junk,junk, format='L,L,A,A,A,A,A,A,A'
        
        for i = 1, 300 do begin
            k = where(fiber eq i, ct)
         
            if ct eq 0 then begin
                l = where(fiber1 eq i)
                j = where(fiber eq -1)
                fiber(j(0)) = i
                dmin(j(0)) = -1.0
                ra(j(0)) = hms2dec(ra1(l(0)))
                dec(j(0)) = hms2dec(dec1(l(0)))
                target(j(0)) = 'target'
            endif
        endfor
        sort1 = sort(fiber)
        fiber = fiber(sort1)
        dmin = dmin(sort1)
        ra = ra(sort1)
        dec = dec(sort1)
        rmag = rmag(sort1)
        rapmag = rapmag(sort1)
        bcode = bcode(sort1)
        icode = icode(sort1)
        newicode = newicode(sort1)
        rcode = rcode(sort1)
        target = target(sort1)
        

    endif

    






    output = string(fiber) + ' ' + string(dmin) + ' '+ string(ra) + ' ' + $
       string(dec) + ' ' + string(rmag) + ' ' + string(rapmag) + ' ' + $
       string(bcode) + ' ' + string(icode) + ' ' +string(rcode) + ' ' + $
       string(newicode)
    

    openw, 1, 'plugdir/' + outname+'_all'
    for i = 0, n_elements(output) -1 do begin
       printf, 1, output(i)
    endfor
    close, 1

    k = where(strmatch(target,'target'))
    if k(0) ne -1 then begin $
       openw, 1, 'plugdir/' + outname+'_targets'
       for i = 0, n_elements(k) -1 do begin
          printf, 1, output(k(i))
       endfor
    endif

    output = string(fiber) +'  '+ target
    close, 1
    k = where(strmatch(target,'*sky*'))
    if k(0) ne -1 then begin $
       openw, 1, 'plugdir/' + outname+'_sky'
       for i = 0, n_elements(k) -1 do begin
          printf, 1, output(k(i))
       endfor
    endif
    close, 1
    
    k = where(strmatch(target,'unused'))
    if k(0) ne -1 then begin $
       openw, 1, 'plugdir/' + outname+'_unused'
       for i = 0, n_elements(k) -1 do begin
          printf, 1, output(k(i))
       endfor
    endif
    close, 1


END
