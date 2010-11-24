; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2010  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
function getsize_fits, filename, nmaps = nmaps, nside = nside, mlpol = mlpol, ordering = order, obs_npix =obs_npix, type = type, header = header, tags = tags, extension=extension_id
;+
;  result = getsize_fits(filename, [nmaps=, nside=, mlpol=, ordering=,
;  obs_npix=, type=, header=, tags=, extension=])
;
;     Get the number of pixels stored in a map FITS file.
;     Each pixel is counted only once 
;     (even if several information is stored on each of them, see nmaps).
;     Depending on the data storage format, result may be :
;       - equal or smaller to the number Npix of Healpix pixels available 
;          over the sky for the given resolution (Npix = 12*nside*nside)
;       - equal or larger to the number of non blank pixels 
;         (obs_npix)
;     
;  INPUT
;     filename = (IN) name of the (compressed) FITS file containing the map
;
;  OPTIONAL INPUT
;     extension = (IN) number of extension (0 based) to read data from
;     /header   = (IN) if set, filename is indeed as FITS header instead of a
;          FITS file
;
;  OPTIONAL OUTPUT
;     nmaps = (OPTIONAL, OUT) number of maps in the file
;     ordering = (OPTIONAL, OUT) Healpix ordering scheme used
;                  '' = unknown
;                  'RING'
;                  'NESTED'
;     obs_npix =  (OPTIONAL, OUT) number of non blanck pixel.    
;                It is set to -1 if it can not be determined from header
;                information alone
;     nside = (OPTIONAL, OUT) Healpix parameter Nside
;     mlpol = (OPTIONAL, OUT) maximum multipole
;     type = (OPTIONAL, OUT) Healpix/FITS file type
;                  <0 : file not found, or not valid
;                  0  : image only fits file, deprecated Healpix format
;                        (result = 12 * nside * nside)
;                  1  : ascii table, generally used for C(l) storage
;                  2  : binary table : with implicit pixel indexing (full sky)
;                        (result = 12 * nside * nside)
;                  3  : binary table : with explicit pixel indexing (generally cut sky)
;                        (result <= 12 * nside * nside)
;                999  : unable to determine the type
;
;      adapted from F90 getsize_fits
;      EH, 2000-11
;      2008-04-01: accepts compressed files
;       Jan 2009: calls init_astrolib
;-

routine = 'getsize_fits'
if n_params() eq 0 then begin
    print,'Syntax : size='+routine+'(filename, nmaps= , nside= , mlpol= , ordering= , obs_npix=, type= , header=, tags=, extension=)'
    message,'Abort'
    return,-1
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

fsize  = 0
type  = -1
nside =  0
nmaps =  0
mlpol =  0
obs_npix = -1
order = ''
fake = 0 ; <<<< bug correction Feb 13, 2003

from_header = 0
if keyword_set(header) then begin
    from_header = 1
    fits_header = filename
endif
from_file = 1- from_header

if (from_file) then begin
; get the primary header information
    hdr = headfits(filename, errmsg = errmsg)
    if (strtrim(errmsg) ne '') then begin
        print,errmsg          ; stop execution
        retall
    endif
    fits_info,filename, /silent, n_ext=n_ext
    file_err = filename
endif else begin
    hdr = fits_header
    fake = 0
    if (n_elements(hdr) eq 1) then begin
       fake = 1
    endif else begin
       junk = sxpar(hdr,'XTENSION', count=n_ext)
       junk = sxpar(hdr,'END', count=n_end,/silent)
       fake = (n_end eq 0)
    endelse
    file_err ='???'
endelse

if (fake) then begin
    if defined(tags) then begin
       select_pix = index_word(tags,'PIXEL',err=error2)
       if error2 eq 1 then begin
          type = 3
          return,fsize
       endif else begin
          type = 2
          return,fsize
       endelse
    endif
    return,fsize
endif
;no extension: backward compatibility with old format
if (n_ext eq 0) then begin
    print,'the file '+file_err+' is in a deprecated format'
    print,' you can update to current HEALPix format using WRITE_FITS_MAP'

    type = 0

    naxis  = sxpar(hdr,'NAXIS')
    naxes = replicate(1L,2)
    naxes[0] = sxpar(hdr,'NAXIS1')
    if naxis ge 2 then naxes[1]= sxpar(hdr,'NAXIS2')
    
    npix = naxes[0]*naxes[1]

    nmaps = min(naxis)

    obs_npix = -1

    nside_fits = npix2nside(npix,err=error)
    if (error eq 0) then nside = nside_fits

    junk = sxpar(hdr,'ORDERING', count = countfits)
    if (countfits ne 0) then begin
        order = junk
    endif else begin
        print,'No ordering information can be found'
    endelse

    return, npix
endif else begin
    type = 999
    if undefined(extension_id) then extension_id = 0
    if (from_file) then begin
        if (extension_id + 1) gt n_ext then begin
            print,' Requested extension ',extension_id,' (0 based) from ',filename
            print,' Only found ',n_ext,' extensions.'
            message,' Abort'
        endif 
                                ; open extension 1 and get header
        errmsg=''
        use_fxb = 1
        fxbopen, lun, filename, extension_id+1, xhdr, errmsg=errmsg
        if (errmsg ne '') then begin ; maybe is it compressed ?
            xhdr = headfits(filename, exten=extension_id+1, /silent)
            use_fxb = 0
        endif
    endif else begin
        xhdr = fits_header
    endelse

    xtension = strtrim(strupcase(sxpar(xhdr,'XTENSION')))
    if (xtension eq 'TABLE') then begin ; ascii table (generally for power spectrum)
        type = 1
    endif

    junk  = ROUND(FLOAT(sxpar(xhdr,'TFIELDS',count=countfits)))
    if (countfits ne 0) then nmaps = junk

    nrows  = ROUND(FLOAT(sxpar(xhdr,'NAXIS2',count=countfits)))
    tforms = strtrim(sxpar(xhdr,'TFORM*'),2)
    sl = strupcase(strmid(tforms,0,1)) ; first letter of 'TFORM*'
    if (sl[0] eq 'E' or sl[0] eq 'I') then begin
        fsize = nrows
    endif else begin
        fsize = nrows * max(long(tforms))
    endelse

    ;if the same keyword is present twice, SXPAR takes the last one
    junk  = STRTRIM(sxpar(xhdr,'ORDERING',count=countfits,/silent),2)
    if (countfits ne 0) then order = junk[0]

    junk  = LONG(sxpar(xhdr,'NSIDE',count=countfits,/silent))
    if (countfits ne 0) then nside = junk[0]

    junk  = LONG(sxpar(xhdr,'MAX-LPOL',count=countfits))
    if (countfits ne 0) then mlpol = junk[0]

    junk  = LONG(sxpar(xhdr,'OBS_NPIX',count=countfits))
    if (countfits ne 0) then obs_npix = junk[0]

        
    ; find the file type
    if (type eq 999) then begin
                                ; most stringent test
        junk = STRUPCASE(STRTRIM(sxpar(xhdr,'INDXSCHM',count = countfits,/silent),2))
        if (countfits ne 0) then begin
            type = 3
            if (junk eq 'IMPLICIT') then type = 2
            goto, found
        endif
                                ; 2nd most stringent test
        junk = LONG(sxpar(xhdr,'GRAIN',count = countfits,/silent))
        if (countfits ne 0) then begin
            type = 3
            if (junk eq 0) then type = 2
            goto, found
        endif
                                ; 3rd most stringent test
        junk = STRUPCASE(STRTRIM(sxpar(xhdr,'TTYPE*',count=countfits),2))
        if (countfits ne 0 and junk[0] ne '') then begin
            type = 2
            if junk[0] eq 'PIXEL' then type = 3
            goto, found
        endif
                                ; lousy test, for backward compatibility
        junk = STRUPCASE(STRTRIM(sxpar(xhdr,'OBJECT',count=countfits),2))
        if (countfits ne 0 and junk[0] ne '') then begin
            if (junk eq 'PARTIAL') then type = 3
            if (junk eq 'FULLSKY') then type = 2
            if (type ne 999) then goto, found
        endif
                                ; very lousy test
        if (nside gt 0) then begin
            npix = nside2npix(nside, error = error)
            if (error eq 0) then begin
                if (npix gt fsize) then type = 3
                if (npix eq fsize) then type = 2
                if (type ne 999) then goto, found
            endif
        endif
        
    endif

found:

endelse

if (from_file && use_fxb) then begin
; close the file
    fxbclose, lun
endif


return, fsize
end

