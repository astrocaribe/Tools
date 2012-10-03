;+
; :Routine Name: line_strength
;
; :Description:
;    Return the coefficients of a gaussian fit to a line,
;    and measure the FWHM of said line.
;
; :Parameters:
; 	:Params:
; 	 outfile  - output file for Gaussian plot
;    wave     - wavelength [arbitrary]
;    flux     - flux [arbitrary]
;
;
; :Output:
;    fwhm - The strength of the spectral line. The units are the
;           same as that input into the function
;           
;    outfile - A hardcopy plot of the spectral region, including the 
;              smoothed line, gaussian fit and the line location 
;              for input spectrum.
;
; :Versions History:
;   0.1 - This version assumes that wavelength and flux are 
;         already ordered (sorted), and the "window" for the
;         line has already been choosen (appropriate wavelength
;         range). This version will generate a plot of the line 
;         that includes a gaussian estimate, and the function 
;         will return the coefficients of the gaussian fit.
;
; :Author:
; 	leblanc
;
; :Date:
; 	Sep 27, 2012
;-
function line_strength, outfile, wave, flux

  ;Convert multi-dimension arrays into 1-dim
  wave = reform(wave)
  flux = reform(flux)
  
  ;Fit a Gaussian to the input line (using the normalized spectrum)
  fit = gaussfit(wave, flux, coeff, nterms = 3)
  
  print
  print, "Gaussian fit coefficients (A, B, C): "
  print, coeff[0], coeff[1], coeff[2]
  print

  ;Determine lower and upper limits of the wings for estimating the
  ;continuum, ± 4sigma
  llim = coeff[1]-4.*coeff[2]
  ulim = coeff[1]+4.*coeff[2]
  ind = where(wave lt llim OR wave gt ulim)
    

  ;Calculate the strength of the line (FWHM) for the input line
  ;Units are the same of the input spectrum
  fwhm = 2.*sqrt(2.*alog(2.))*coeff[2]
  
  
  ;Generate a plot of the line with the estimated Gaussian fit
  xlim = [min(wave), max(wave)]
  pos = idl_setplot(outfile, 16., 9.)   ;sets up and opens a PS device
  leg = ['Original', 'Gaussian Fit', 'Line Position']

  plot, wave, flux, xtit = 'Wavelength', ytit = 'Flux', xr = xlim, $
          xstyle = 1, position=pos, yr=[min(flux),max(flux)], /nodata
  oplot, wave, smooth(flux, 10), color = cgColor('dark green'), thick = 3.
  oplot, wave, fit, color = cgColor('red'), thick = 1.5
  oplot, [coeff[1], coeff[1]], [!Y.CRANGE[0], !Y.CRANGE[1]], $
         linestyle = 2, color = cgColor('black'), thick = 2.       

  al_legend, leg, psym=[-0, -0, -0], linestyle = [0,0,2], $
             color=['dark green', 'red', 'black']
        
  line_pos = strtrim(string(coeff[1], format = '(f10.3)') ,2)    
  line_str = strtrim(string(fwhm, format = '(f8.4)') ,2)
        
  xyouts, 0.7, 0.42, 'Position: ' + line_pos, /normal      
  xyouts, 0.7, 0.40, 'Strength: ' + line_str, /normal      
  
  ;Close the graphics plotting device
  device, /close
  
  
  ;Return the calculation on the strength of the line (FWHM)
  return, fwhm
    
end