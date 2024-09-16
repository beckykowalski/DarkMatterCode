# DarkMatterCode

To calculate event rate:


   RunScripts.py
   Input for CalculateEventRate:


   	 cross section, M_DM, List of elements, starting energy, ending energy, bins in calculation, resolution



   Function returns for the list of energies provided:
   	 list of energies, list of Sm, list of S0, list of Form Factor, list of x values
	 currently <list of x values> is empty, should be filled with v_min/v0 eventually for debugging purposes


   An intermediate function exists: RunEventRateCalculation.py, to make a plot and text file of calculated event rate




To calculate AM sensitivity:



   RunSensitivity.py
   Input for SensitivityCalculation:


   	 cross section, M_DM, exposure, CL value, list of elements, bkg counts/keV/day, starting energy, ending energy, bins in calculation, resolution



   Change name of txt (line 33) to save mass, cross section in a text file

