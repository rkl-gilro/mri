/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon                    */



/***************************************************************************
*              3D Adaptive Multiresolution Non-Local Means Filter          *
* Pierrick Coupe a, Jose V. Manjon, Montserrat Robles , D. Louis Collins   *
***************************************************************************/



/*                          Details on ONLM filter                        */
/***************************************************************************
 *  The ONLM filter is described in:                                       *
 *                                                                         *
 *  P. Coupé, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, * 
 *  Avril 2008                                                             *
 ***************************************************************************/

/*                          Details on Wavelet mixing                     */
/***************************************************************************
 *  The hard wavelet subbands mixing is described in:                      *
 *                                                                         *
 *  P. Coupé, S. Prima, P. Hellier, C. Kervrann, C. Barillot.              *
 *  3D Wavelet Sub-Bands Mixing for Image Denoising                        *
 *  International Journal of Biomedical Imaging, 2008                      * 
 ***************************************************************************/

/*                      Details on Rician adaptation                      */
/***************************************************************************
 *  The adaptation to Rician noise is described in:                        *
 *                                                                         *
 *  N. Wiest-Daesslé, S. Prima, P. Coupé, S.P. Morrissey, C. Barillot.     *
 *  Rician noise removal by non-local means filtering for low              *
 *  signal-to-noise ratio MRI: Applications to DT-MRI. In 11th             *
 *  International Conference on Medical Image Computing and                *
 *  Computer-Assisted Intervention, MICCAI'2008,                           *
 *  Pages 171-179, New York, États-Unis, Septembre 2008                    *
 ***************************************************************************/



Matlab code for 3D Adaptive Multiresolution Non-Local Means Filter

1) Select the directory where the files have been extracted 
as work directory in matlab.

2) compile the mex codes of the filters
In the matlab prompt:

	> mex onlm.c

for the Gaussian version.

and 
	> mex ornlm.c

for the Rician version.

3) Open a script : GaussianNoise.m or RicianNoise.m 
(scripts used in experiment section of the paper)

4) Choose your modality by commenting the two others

	name ='data/t1_icbm_normal_1mm_pn0_rf0.rawb';
	%  name ='data/t2_icbm_normal_1mm_pn0_rf0.rawb';
	% name ='data/pd_icbm_normal_1mm_pn0_rf0.rawb';

5) Chosse the size of the 3D volume

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Comment this part to filter the entire image.
	ima=ima(:,:,70:100);
	Label=Label(:,:,70:100);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Comment this to work on the entire image 
(it takes several hours for all the level of noise with all the filter).

6) Run the script and go to take a coffee

