# Improved-TOA-Localization-Using-MPR

Existing communication systems are required to provide highly precise locations of users or terminals to support diverse applications. This study enhances the performance of time of arrival localization by utilizing modified polar representation, regardless of the range of users or terminals. We propose two solutions: 1) the two-stage semidefinite relaxation (2S-SDR), enhancing accuracy by refining the weighting matrix after initial SDR; 2) the maximum likelihood estimator (MLE), using Gauss-Newton iteration to minimize its cost function. Both solutions achieve the Cramér-Rao lower bound (CRLB) under moderate noise levels. Particularly, the MLE demonstrates CRLB accuracy and lower bias in the low noise region. Simulation results validate the effectiveness of the proposed solutions.

**If you use any of the following codes in your research, please cite the paper as a reference in your publication. Thank you!**

## Improved TOA Localization Using MPR (07/18/2024)

### <u>Reference</u>
>C. Si, R. Fan, Y. Yang and Y. Sun, “Improved TOA Localization Using Modified Polar Representation,” *IEEE Commun. Lett.*, vol. 28, no. 9, pp. 2051-2055, Sep. 2024.

### Code List:
- TOA_MPR_2SSDR.m: proposed 2S-SDR Solution.
- TOA_MPR_MLE.m: MLE implemented by Gauss-Newton iteration for TOA localization in MPR.
- CRLB_TOA.m: CRLB of TOA localization in MPR.
- TOA_2WLS_Cart.m: Cartesian TOA localization, from (Z. Ma and K. C. Ho, "TOA localization in the presence of random sensor position errors," in 2011 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2011, pp. 2468-2471.).
- mainTOAMPR.m: main function of generating Fig. 2 and Fig. 3.
