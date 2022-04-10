# Energy Calibration of Neutron Detectors
This code calibrates the neutron detectors used for the Musun experiment.</br>
Neutron detectors see a lot gamma rays due to muon decay. Thus a convolution of Compton scattering and Klein Nishina formula is suitable used for 
calibrating the detectors. 

## Compton scattering
Gamma ray sources undergo Compton scattering and so they were used for the energy calibration of the neutron detectors. The specific energy of the 
Compton edge of the gamma ray, due to Compton scattering was compared with the number of channels read by the neutron counter corresponding to that energy.
The Compton scattering wavelength ∆λ is given by, </br>
<img width="195" alt="Screen Shot 2022-04-10 at 8 57 28 AM" src="https://user-images.githubusercontent.com/27436642/162619302-5defbe1f-15e9-4138-b1db-e2f9ae57c2d5.png"></br>
where me is the mass of the electron and θ is the angle by which the incident gamma ray scatters. 
At the maximum scattering angle θ = 180°, the energy difference carried away by the scattered electron is its maximum kinetic energy ECompton given by the equation,</br>
<img width="224" alt="Screen Shot 2022-04-10 at 8 58 24 AM" src="https://user-images.githubusercontent.com/27436642/162619332-d0760055-92a9-4ea3-bd3d-76357e6769de.png"></br>
The radioactive sources used during the experiment for calibration runs were 60Co and 137Cs. 60Co is an unstable isotope of Cobalt with a half life of
about 5.27 years which undergoes a beta decay to an excited state of 60Ni, which in turn emits gamma rays with energies 1.17 and 1.33 MeV respectively. 
Since these spectral lines are in close proximity we took the average of these to be an equivalent spectral line emitted by this source, yielding an of 
energy 1.25 MeV. 137Cs undergoes a beta decay to 137Ba which emits a distinct photo peak of 662 keV.

## Klein Nishina Formula
Klein Nishina formula gives the expression for differential scattering cross section of photons scattered by single electron corresponding to the lowest 
order of quantum electrodynamics.</br> This also considered the details of relativistic quantum mechanical effect in the high energy regime and so we 
expected it to give accurate results. The energy distribution of the differential cross section is given by,</br>
<img width="555" alt="Screen Shot 2022-04-10 at 9 22 48 AM" src="https://user-images.githubusercontent.com/27436642/162620258-f78c0e8f-a82a-4a7f-8c84-860b493a9816.png"></br>
## Fit function
The convolution is done numerically using C++ and ROOT using this code. 
A Gaussian function of suitable width w, denoted by G(w) was superimposed, shifted by a bin over the Klein Nishina formula denoted by 
K (i.e. Eq.(5.12)) and finally integrated over the entire range of interest from E1 to E2. This procedure is illustrated by the equation,</br>
<img width="287" alt="Screen Shot 2022-04-10 at 9 24 10 AM" src="https://user-images.githubusercontent.com/27436642/162620323-8c244be9-75d8-4617-91aa-93ffab9de625.png"></br>

The channel corresponding to the average Compton edge returned by the fit function was used as the final calibration gain from the two source.
The example below show the energy distribution from one neutron counter fitted with the convolution function</br>
<img width="625" alt="Screen Shot 2022-04-10 at 9 31 35 AM" src="https://user-images.githubusercontent.com/27436642/162620607-6d1df287-48bb-4aa8-8e56-4c2e1a53a1c0.png">

This is a simple C code which just need ROOT int to produce the results.
